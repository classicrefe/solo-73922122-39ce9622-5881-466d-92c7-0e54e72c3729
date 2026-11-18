"""Reference FASTA parsing and alternate-haplotype construction.

The ``/haplotype`` endpoint overlays the records that intersect one interval
of a pinned VCF version onto an uploaded single-chromosome reference FASTA.
Only simple biallelic records are supported: single-base substitutions plus
VCF-anchored insertions/deletions (one allele is one base and both alleles
share their first base).  All records are treated as one haplotype, so a
duplicate, an overlapping reference span, a REF mismatch or an unsupported
record rejects the whole request.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import pysam

MAX_INTERVAL = 100_000
# Reference bases must be uppercase ACGT; soft-masked / ambiguity codes are
# deliberately not accepted because they cannot be checked against VCF alleles.
_SEQ_RE = re.compile(r"[ACGT]+")

SNV = "snv"
INSERTION = "insertion"
DELETION = "deletion"


class FastaError(ValueError):
    """The uploaded reference FASTA violates the acceptance rules."""


class HaplotypeError(ValueError):
    """The requested interval or record set cannot form one haplotype."""


@dataclass(frozen=True)
class Fasta:
    name: str
    sequence: str


@dataclass(frozen=True)
class Variant:
    line: str
    pos: int
    ref: str
    alt: str
    kind: str | None

    @property
    def span_end(self) -> int:
        """1-based inclusive end of the record's reference span."""
        return self.pos + len(self.ref) - 1


def parse_fasta(data: bytes) -> Fasta:
    """Parse a single-record, uppercase ACGT FASTA held in memory."""
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise FastaError(f"upload is not valid UTF-8 text: {exc}")

    lines = text.lstrip("\ufeff").splitlines()
    first = 0
    while first < len(lines) and not lines[first].strip():
        first += 1
    if first == len(lines):
        raise FastaError("upload is empty")

    header = lines[first].strip()
    if not header.startswith(">"):
        raise FastaError("missing '>' FASTA header line")
    name_parts = header[1:].split()
    if not name_parts:
        raise FastaError("FASTA header does not name a chromosome")
    name = name_parts[0]

    seq_parts = []
    for raw in lines[first + 1:]:
        if raw.startswith(">"):
            raise FastaError("FASTA must contain exactly one chromosome record")
        line = raw.strip()
        if line:
            seq_parts.append(line)
    sequence = "".join(seq_parts)
    if not sequence:
        raise FastaError("FASTA contains no sequence")
    if not _SEQ_RE.fullmatch(sequence):
        bad = next(ch for ch in sequence if ch not in "ACGT")
        raise FastaError(f"sequence must contain only A/C/G/T (found {bad!r})")
    return Fasta(name=name, sequence=sequence)


def check_interval(start: int, end: int, seq_len: int) -> None:
    """Validate the requested 1-based closed interval against the reference."""
    if end < start:
        raise HaplotypeError("end must be >= start")
    if end - start + 1 > MAX_INTERVAL:
        raise HaplotypeError(f"interval may cover at most {MAX_INTERVAL} bases")
    if end > seq_len:
        raise HaplotypeError(
            f"interval [{start}, {end}] extends past the {seq_len}-base reference sequence"
        )


def resolve_chrom(requested: str | None, fasta_name: str, dataset_chrom: str) -> str:
    """Require the requested, FASTA and published chromosomes to agree."""
    target = requested if requested else fasta_name
    if target != fasta_name:
        raise HaplotypeError(
            f"requested chromosome '{target}' does not match FASTA record '{fasta_name}'"
        )
    if target != dataset_chrom:
        raise HaplotypeError(
            f"FASTA record '{target}' does not match the published dataset "
            f"chromosome '{dataset_chrom}'"
        )
    return target


def fetch_intersecting(gz_path: str, chrom: str, start: int, end: int) -> list[str]:
    """Read records whose reference span intersects [start, end] via Tabix."""
    tb = pysam.TabixFile(gz_path)
    try:
        # 1-based closed [start, end] -> 0-based half-open fetch.  htslib's
        # VCF index stores the full REF span, so deletions starting left of
        # the interval are returned too.
        raw = list(tb.fetch(chrom, start - 1, end))
    finally:
        tb.close()

    records = []
    for line in raw:
        fields = line.split("\t", 4)
        if len(fields) < 5:
            # Let record validation report the malformed line.
            records.append(line)
            continue
        try:
            pos = int(fields[1])
        except ValueError:
            records.append(line)
            continue
        ref = fields[3]
        if pos <= end and pos + len(ref) - 1 >= start:
            records.append(line)
    return records


def _classify(pos: int, ref: str, alt: str, issues: list[str]) -> str | None:
    if not ref or not alt:
        issues.append(f"malformed variant at POS {pos}: empty REF or ALT")
        return None
    if len(ref) == 1 and len(alt) == 1:
        if ref == alt:
            issues.append(f"unsupported variant at POS {pos}: REF and ALT are identical")
            return None
        return SNV
    if len(ref) == 1:
        if alt[0] != ref[0]:
            issues.append(
                f"unsupported variant at POS {pos}: insertion ALT '{alt}' does not "
                f"start with the REF anchor base '{ref}'"
            )
            return None
        return INSERTION
    if len(alt) == 1:
        if ref[0] != alt[0]:
            issues.append(
                f"unsupported variant at POS {pos}: deletion REF '{ref}' does not "
                f"start with the ALT anchor base '{alt}'"
            )
            return None
        return DELETION
    issues.append(
        f"unsupported variant at POS {pos}: only single-base substitutions and "
        "anchored insertions/deletions are supported"
    )
    return None


def validate_records(lines: list[str], reference: str) -> list[Variant]:
    """Parse, check and order the records that will form the haplotype.

    Collects every problem and raises one error so the caller can reject the
    whole request with an explanation.
    """
    variants: list[Variant] = []
    issues: list[str] = []
    seen: set[tuple[int, str, str]] = set()

    for line in lines:
        fields = line.split("\t")
        if len(fields) != 8:
            issues.append(f"malformed VCF record (expected 8 fields): {line[:80]!r}")
            continue
        try:
            pos = int(fields[1])
        except ValueError:
            issues.append(f"malformed VCF record (bad POS): {line[:80]!r}")
            continue
        ref, alt = fields[3], fields[4]
        kind = _classify(pos, ref, alt, issues)

        key = (pos, ref, alt)
        if key in seen:
            issues.append(f"duplicate record at POS {pos} ({ref}>{alt})")
        else:
            seen.add(key)

        variant = Variant(line=line, pos=pos, ref=ref, alt=alt, kind=kind)
        variants.append(variant)

        if pos < 1 or variant.span_end > len(reference):
            issues.append(
                f"reference sequence does not cover POS {pos}-{variant.span_end} "
                f"required by record ({ref}>{alt})"
            )
        else:
            actual = reference[pos - 1:variant.span_end]
            if actual != ref:
                issues.append(
                    f"REF mismatch at POS {pos}: reference has '{actual}', "
                    f"VCF has '{ref}'"
                )

    ordered = sorted(variants, key=lambda v: (v.pos, v.span_end))
    furthest: Variant | None = None
    for variant in ordered:
        if furthest is not None and variant.pos <= furthest.span_end:
            issues.append(
                "overlapping reference spans: "
                f"POS {furthest.pos} ({furthest.pos}-{furthest.span_end}) and "
                f"POS {variant.pos} ({variant.pos}-{variant.span_end})"
            )
        if furthest is None or variant.span_end > furthest.span_end:
            furthest = variant

    if issues:
        # Keep the first occurrence of each explanation, in a stable order.
        raise HaplotypeError("; ".join(dict.fromkeys(issues)))
    return ordered


def apply_variants(reference: str, start: int, end: int,
                   variants: list[Variant]) -> str:
    """Build the alternate sequence for [start, end] from ordered variants.

    Coordinates are applied against the reference.  An insertion is assigned
    to its left anchor base (the VCF POS), so an insertion anchored inside the
    interval is included after that base.  Deletions clip to the interval, so
    a deletion crossing either boundary affects the overlapping bases.
    """
    substitutions: dict[int, str] = {}
    insertions: dict[int, str] = {}
    deleted: list[tuple[int, int]] = []

    for variant in variants:
        if variant.kind == SNV:
            substitutions[variant.pos] = variant.alt
        elif variant.kind == INSERTION:
            insertions[variant.pos] = variant.alt[1:]
        elif variant.kind == DELETION:
            lo = max(start, variant.pos + 1)   # ALT keeps the anchor base
            hi = min(end, variant.span_end)
            if lo <= hi:
                deleted.append((lo, hi))

    deleted.sort()
    delete_i = 0
    out: list[str] = []
    for pos in range(start, end + 1):
        while delete_i < len(deleted) and deleted[delete_i][1] < pos:
            delete_i += 1
        if delete_i < len(deleted) and deleted[delete_i][0] <= pos:
            continue
        out.append(substitutions.get(pos, reference[pos - 1]))
        if pos in insertions:
            out.append(insertions[pos])
    return "".join(out)
