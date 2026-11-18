# Versioned VCF Variant Store

FastAPI + pysam + SQLite service for uploading single-chromosome VCF files and
querying them by 1-based closed chromosome intervals, with atomic versioned
replacement and safe concurrent access from two app instances. It also builds
alternate DNA sequences for a reference interval by overlaying the published
variants as one haplotype (`POST /haplotype`).

## Run

```bash
docker compose up --build
```

Two instances come up sharing one data volume:

- instance 1: http://localhost:8001/docs
- instance 2: http://localhost:8002/docs

Both serve the same published data; uploads may go to either one. `/docs` is
the interactive entry point for every endpoint (including `/haplotype`).

## API

### `POST /vcf?base_version=<n>` (multipart file field `file`)

Accepts a VCF that is:

- at most 10 MiB,
- single-chromosome (all records share one CHROM),
- exactly 8 columns (`#CHROM POS ID REF ALT QUAL FILTER INFO`, no FORMAT/sample
  columns),
- biallelic with `REF`/`ALT` containing only `A/C/G/T` (uppercase; no `N`, no
  symbolic alleles, no `*`, no multi-allelic `A,G`).

Records are sorted by position, BGZF-compressed and Tabix-indexed. The upload
carries the version the client last saw (`base_version`, `0` when the store is
empty). If two replacements race — even across the two instances — exactly one
wins; the loser gets `409` and must re-read `/version` and retry.
Returns `201` with `{"version", "chrom", "n_records"}`.

### `GET /variants?chrom=<c>&start=<s>&end=<e>`

1-based closed interval `[start, end]`. Returns every raw VCF record whose
reference span `[POS, POS+len(REF)-1]` intersects the interval — including
deletions that start before the interval. Fetching goes through the Tabix
index (`TabixFile.fetch`), never a full scan. The response carries the
`version` that was actually read; that version's data and index are pinned for
the duration of the query, so version, data and index always match.

### `POST /haplotype?start=<s>&end=<e>[&chrom=<c>][&version=<v>]`

Multipart file field `file`: a reference FASTA (max 10 MiB) containing exactly
one record with an uppercase `A/C/G/T` sequence. The header's first word names
the chromosome and must match the published dataset; `chrom` may be given
explicitly (it must agree with both) or omitted.

`start`/`end` are 1-based, inclusive, and cover at most 100,000 bases. Every
published record whose reference span intersects the interval is read through
the Tabix index and treated as one haplotype. Supported records:

- single-base substitution (REF and ALT both length 1),
- insertion (REF length 1; ALT starts with the REF base) — assigned to its left
  anchor base, i.e. included when `POS` is inside the interval, not when only
  `POS+1` is,
- deletion (ALT length 1; REF starts with the ALT base) — the anchor base is
  kept, the remaining REF bases are removed; a deletion crossing either
  interval boundary deletes only the bases inside the interval.

Any duplicate record, overlapping reference spans (touching spans are fine),
REF mismatch against the uploaded FASTA, reference sequence too short for a
record, or unsupported/complex record rejects the whole request with `400` and
an explanation listing every problem found. Nothing is published or modified.

The response is:

```json
{
  "version": 1,
  "chrom": "chr1",
  "start": 100,
  "end": 110,
  "reference": "ACGTAGAAAAC",
  "sequence": "AATAAAACT",
  "count": 3,
  "records": [
    "chr1\t100\t.\tACGT\tA\t.\t.\t.",
    "chr1\t105\t.\tG\tT\t.\t.\t.",
    "chr1\t110\t.\tC\tCT\t.\t.\t."
  ]
}
```

Here the deletion at 100 keeps its anchor `A` and removes `CGT`, the SNP at 105
turns `G` into `T`, and the insertion at 110 appends `T` after its anchor `C`.

`version` is the VCF version the sequence was built from and `records` lists
the intersecting VCF records that participated (sorted by position).

The version is pinned when the request starts, before the FASTA upload is read.
By default the version current at that moment is used for the whole request, so
another instance publishing a correction mid-request cannot mix versions.
Passing `version=<v>` explicitly selects a fixed published version (e.g. to
reproduce a build while a correction is being prepared); a version that was
never published or has already been reclaimed answers `404`. Pins keep that
version's files alive until the request finishes. A failed request never
changes published data, and the endpoint never publishes anything.

### `GET /version`

Current published version, chromosome and record count.

## Consistency design

- **Atomic publish**: a new version is fully built (sorted, compressed,
  indexed) in a per-instance staging directory, moved into `versions/`, and
  only then the current-version pointer in SQLite is flipped inside a single
  `BEGIN IMMEDIATE` transaction. Readers never see a half-built file set.
- **Optimistic concurrency**: the pointer flip checks `base_version` inside
  the same transaction, so concurrent replacements from both instances are
  serialized and only one succeeds.
- **Old queries finish**: a query pins its version in a `refs` table before
  opening files. A background reaper deletes a superseded version's data and
  index only when no pin references it. Instances heartbeat; pins of crashed
  instances expire and never block reclamation.
- **Crash recovery**: if the process dies mid-build or mid-publish, the
  staging/orphan directories left behind are swept (with a grace period) and
  the previously published version keeps serving. The pointer flip is the last
  step, so a restart never exposes a partial file set.
- **Fixed versions for haplotype builds**: `/haplotype` pins the version it
  starts with (current by default, or an explicit `version`) before opening
  the index. The reaper keeps a superseded version on disk while such a
  request is in flight, so instance 2 can publish a correction while instance
  1 still finishes a build from version N. `/haplotype` is read-only and never
  publishes, so a failed build leaves the published data untouched.

## Local development

```bash
pip install -r requirements.txt
python -m pytest tests/          # 84 tests
uvicorn --factory app.main:create_app   # single dev instance (uses /data or $VCF_DATA_DIR)
```

Configuration via environment: `VCF_DATA_DIR`, `VCF_INSTANCE_ID`,
`VCF_MAX_UPLOAD_BYTES`, `VCF_REAPER_INTERVAL`, `VCF_HEARTBEAT_INTERVAL`,
`VCF_DEAD_AFTER`, `VCF_ORPHAN_GRACE`.
