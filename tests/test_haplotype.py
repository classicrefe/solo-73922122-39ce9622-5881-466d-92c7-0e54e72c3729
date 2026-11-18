"""Alternate-haplotype endpoint: FASTA upload, record checks, sequence build."""
import os

import pytest

from conftest import SAMPLE, upload, vcf

from app.haplotype import MAX_INTERVAL


def ref(length: int = 300, overrides: dict[int, str] | None = None) -> str:
    """Reference of `length` A's with 1-based segments replaced."""
    bases = bytearray(b"A" * length)
    for pos, text in (overrides or {}).items():
        bases[pos - 1:pos - 1 + len(text)] = text.encode()
    return bases.decode()


def fasta(sequence: str, header: str = "chr1") -> bytes:
    body = "\n".join(sequence[i:i + 60] for i in range(0, len(sequence), 60))
    return f">{header}\n{body}\n".encode()


def sample_reference() -> str:
    """Reference matching SAMPLE: del 100-103, SNP G105T, ins C110>CT, del 200-204."""
    return ref(overrides={100: "ACGT", 105: "G", 110: "C", 200: "TAAAA"})


def haplotype(client, sequence: str, start: int, end: int, **params):
    data = params.pop("body", None)
    if data is None:
        data = fasta(sequence, params.pop("header", "chr1"))
    query = {"start": start, "end": end, **params}
    return client.post(
        "/haplotype",
        params=query,
        files={"file": ("reference.fa", data)},
    )


def test_snv_insertion_and_deletion_are_applied(client):
    upload(client, SAMPLE, base_version=0)
    sequence = sample_reference()

    r = haplotype(client, sequence, 95, 120)
    assert r.status_code == 200, r.text
    body = r.json()
    # 95-100 = AAAAAA, 101-103 (CGT) deleted, 104 A, 105 G>T, 106-109 AAAA,
    # 110 C + inserted T, 111-120 A * 10.
    expected = "AAAAAA" + "A" + "T" + "AAAA" + "CT" + "A" * 10
    assert body["sequence"] == expected
    assert body["reference"] == sequence[94:120]
    assert body["version"] == 1
    assert body["chrom"] == "chr1"
    assert body["start"] == 95 and body["end"] == 120
    assert body["count"] == 3
    assert [rec.split("\t")[1] for rec in body["records"]] == ["100", "105", "110"]


def test_records_outside_interval_are_ignored(client):
    upload(client, SAMPLE, base_version=0)
    sequence = sample_reference()
    body = haplotype(client, sequence, 95, 120).json()
    assert "200" not in [rec.split("\t")[1] for rec in body["records"]]
    # The other deletion lies wholly outside and must not touch the sequence.
    assert len(body["sequence"]) == 26 - 3 + 1


def test_deletion_crossing_left_boundary(client):
    upload(client, SAMPLE, base_version=0)
    sequence = sample_reference()
    body = haplotype(client, sequence, 102, 108).json()
    # 100-103 is deleted; only 102-103 fall inside [102, 108]:
    # remaining reference is 104-108 = AGAAA, SNP 105 G>T -> ATAAA.
    assert body["sequence"] == "ATAAA"
    assert body["records"][0].split("\t")[1] == "100"


def test_deletion_crossing_right_boundary(client):
    upload(client, SAMPLE, base_version=0)
    sequence = sample_reference()
    body = haplotype(client, sequence, 197, 202).json()
    # deletion 200 (TAAAA>T) removes 201-204; clip at the interval end.
    assert body["sequence"] == "AAAT"
    assert [rec.split("\t")[1] for rec in body["records"]] == ["200"]


def test_deletion_entirely_inside_interval(client):
    upload(client, SAMPLE, base_version=0)
    sequence = sample_reference()
    body = haplotype(client, sequence, 198, 206).json()
    assert body["sequence"] == "AA" + "T" + "AA"
    assert body["count"] == 1


def test_insertion_belongs_to_left_anchor(client):
    upload(client, SAMPLE, base_version=0)
    sequence = sample_reference()
    # Anchor at 110: query starting at 110 keeps it ...
    assert haplotype(client, sequence, 110, 112).json()["sequence"] == "CTAA"
    # ... a query starting at 111 must not see the insertion.
    assert haplotype(client, sequence, 111, 112).json()["sequence"] == "AA"


def test_insertion_anchored_at_interval_end_is_included(client):
    upload(client, SAMPLE, base_version=0)
    sequence = sample_reference()
    # 110 is the last base of the interval and the left anchor of C>CT.
    body = haplotype(client, sequence, 106, 110).json()
    assert body["sequence"] == "AAAACT"


def test_deletion_covering_whole_interval_gives_empty_sequence(client):
    upload(client, vcf([(100, "ACGTACGTA", "A")]), base_version=0)
    sequence = ref(200, {100: "ACGTACGTA"})
    # REF spans 100-108; 101-108 are deleted by ALT "A".
    assert haplotype(client, sequence, 101, 108).json()["sequence"] == ""
    # the anchor base itself is kept when it is inside the interval
    assert haplotype(client, sequence, 100, 108).json()["sequence"] == "A"


def test_interval_without_records_returns_reference(client):
    upload(client, SAMPLE, base_version=0)
    sequence = sample_reference()
    body = haplotype(client, sequence, 130, 140).json()
    assert body["sequence"] == sequence[129:140]
    assert body["count"] == 0
    assert body["records"] == []


def test_insertion_at_interval_start_and_left_of_it(client):
    body_vcf = vcf([(100, "A", "AG"), (110, "C", "CT")])
    upload(client, body_vcf, base_version=0)
    sequence = ref(200, {100: "A", 110: "C"})
    # insertion anchored exactly at the first base of the interval is kept
    result = haplotype(client, sequence, 100, 101).json()
    assert result["sequence"] == "AGA"
    # an insertion anchored one base before the interval is not included
    result = haplotype(client, sequence, 111, 112).json()
    assert result["sequence"] == "AA"


def test_multi_base_insertion(client):
    upload(client, vcf([(105, "G", "GATTAC")]), base_version=0)
    sequence = ref(200, {105: "G"})
    body = haplotype(client, sequence, 104, 106).json()
    assert body["sequence"] == "AGATTACA"


def test_duplicate_records_rejected(client):
    body = vcf([(105, "G", "T"), (105, "G", "T")])
    upload(client, body, base_version=0)
    r = haplotype(client, sample_reference(), 100, 115)
    assert r.status_code == 400
    assert "duplicate" in r.json()["detail"]
    assert "POS 105" in r.json()["detail"]


def test_overlapping_reference_spans_rejected(client):
    upload(client, vcf([(100, "ACGT", "A"), (102, "GT", "G")]), base_version=0)
    r = haplotype(client, sample_reference(), 100, 115)
    assert r.status_code == 400
    assert "overlapping reference spans" in r.json()["detail"]


def test_three_way_overlap_reports_each_collision(client):
    # 100-103 vs 102-104 overlap; 102-104 vs 104-107 share base 104.
    upload(
        client,
        vcf([(100, "ACGT", "A"), (102, "GTA", "G"), (104, "AGAA", "A")]),
        base_version=0,
    )
    r = haplotype(client, sample_reference(), 99, 120)
    assert r.status_code == 400
    detail = r.json()["detail"]
    assert detail.count("overlapping reference spans") == 2
    assert "POS 100 (100-103)" in detail
    assert "POS 102 (102-104)" in detail
    assert "POS 104 (104-107)" in detail


def test_adjacent_records_are_allowed(client):
    upload(client, vcf([(100, "ACGT", "A"), (104, "A", "C")]), base_version=0)
    r = haplotype(client, sample_reference(), 100, 110)
    assert r.status_code == 200, r.text
    # 101-103 deleted, then the SNV at 104; 105-110 unchanged (GAAAAC).
    assert r.json()["sequence"] == "ACGAAAAC"


def test_ref_mismatch_rejected(client):
    upload(client, vcf([(105, "C", "T")]), base_version=0)
    r = haplotype(client, sample_reference(), 100, 115)
    assert r.status_code == 400
    assert "REF mismatch at POS 105" in r.json()["detail"]
    assert "'G'" in r.json()["detail"] and "'C'" in r.json()["detail"]


def test_reference_not_covering_record_rejected(client):
    # REF spans 298-302 but the reference ends at 300.
    upload(client, vcf([(298, "AAAAA", "A")]), base_version=0)
    r = haplotype(client, ref(300), 290, 300)
    assert r.status_code == 400
    assert "does not cover" in r.json()["detail"]


UNSUPPORTED = {
    "mnv": [(105, "GA", "CT")],
    "insertion_not_anchored": [(105, "G", "AT")],
    "deletion_not_anchored": [(105, "GA", "T")],
    "ref_equals_alt": [(105, "G", "G")],
    "complex": [(105, "GA", "CTTT")],
}


@pytest.mark.parametrize("name", sorted(UNSUPPORTED))
def test_unsupported_variants_rejected(client, name):
    upload(client, vcf(UNSUPPORTED[name]), base_version=0)
    r = haplotype(client, sample_reference(), 100, 115)
    assert r.status_code == 400
    assert "unsupported variant at POS 105" in r.json()["detail"]


def test_all_problems_reported_in_one_rejection(client):
    upload(client, vcf([(105, "GA", "CT"), (110, "T", "A")]), base_version=0)
    r = haplotype(client, sample_reference(), 100, 115)
    assert r.status_code == 400
    detail = r.json()["detail"]
    assert "unsupported variant at POS 105" in detail
    assert "REF mismatch at POS 110" in detail


def test_chromosome_must_match_fasta_and_dataset(client):
    upload(client, SAMPLE, base_version=0)
    r = haplotype(client, sample_reference(), 100, 110, chrom="chr2")
    assert r.status_code == 400
    assert "chr2" in r.json()["detail"]

    r = haplotype(client, ref(300), 100, 110, header="chrX")
    assert r.status_code == 400
    assert "FASTA record 'chrX'" in r.json()["detail"]
    assert "'chr1'" in r.json()["detail"]

    # matching chrom explicitly is fine, as is omitting it
    assert haplotype(client, sample_reference(), 100, 110, chrom="chr1").status_code == 200
    assert haplotype(client, sample_reference(), 100, 110).status_code == 200


def test_interval_validation(client):
    upload(client, SAMPLE, base_version=0)
    sequence = sample_reference()

    r = haplotype(client, sequence, 110, 100)
    assert r.status_code == 400 and "end must be >= start" in r.json()["detail"]

    r = haplotype(client, sequence, 1, 0)
    assert r.status_code == 422  # 1-based query constraint

    r = haplotype(client, sequence, 1, 301)
    assert r.status_code == 400 and "extends past" in r.json()["detail"]

    long_seq = ref(MAX_INTERVAL + 100)
    r = haplotype(client, long_seq, 1, MAX_INTERVAL + 1)
    assert r.status_code == 400 and "at most" in r.json()["detail"]


def test_maximum_interval_accepted(client):
    upload(client, vcf([(1, "A", "C")]), base_version=0)
    r = haplotype(client, ref(MAX_INTERVAL + 100), 1, MAX_INTERVAL)
    assert r.status_code == 200, r.text
    body = r.json()
    assert len(body["sequence"]) == MAX_INTERVAL
    assert body["sequence"][0] == "C"


BAD_FASTAS = {
    "empty": b"",
    "no_header": b"ACGT\n",
    "empty_sequence": b">chr1\n",
    "two_records": b">chr1\nACGT\n>chr2\nACGT\n",
    "lowercase": b">chr1\nACGTacgt\n",
    "ambiguity_code": b">chr1\nACGTN\n",
    "not_utf8": b">chr1\n\xff\xfeACGT\n",
    "header_without_name": b">\nACGT\n",
}


@pytest.mark.parametrize("name", sorted(BAD_FASTAS))
def test_invalid_fasta_rejected(client, name):
    upload(client, SAMPLE, base_version=0)
    r = haplotype(client, "", 1, 10, body=BAD_FASTAS[name])
    assert r.status_code == 400, f"{name}: {r.status_code} {r.text}"
    assert "invalid FASTA" in r.json()["detail"]


def test_oversize_fasta_rejected(client, app):
    upload(client, SAMPLE, base_version=0)
    limit = app.state.cfg.max_upload_bytes
    body = b">chr1\n" + b"A" * limit
    assert len(body) > limit
    r = haplotype(client, "", 1, 10, body=body)
    assert r.status_code == 413


def test_exactly_at_limit_accepted(client, app):
    upload(client, SAMPLE, base_version=0)
    limit = app.state.cfg.max_upload_bytes
    header = b">chr1\n"
    seq = b"A" * (limit - len(header))
    r = haplotype(client, "", 1, 10, body=header + seq)
    assert r.status_code == 200, r.text


def test_haplotype_without_published_data_is_404(client):
    r = haplotype(client, ref(100), 1, 10)
    assert r.status_code == 404


def test_unknown_version_is_404(client):
    upload(client, SAMPLE, base_version=0)
    r = haplotype(client, sample_reference(), 100, 110, version=7)
    assert r.status_code == 404
    assert "version 7" in r.json()["detail"]


def test_pins_current_version_by_default(client):
    upload(client, SAMPLE, base_version=0)
    upload(client, vcf([(105, "G", "C")]), base_version=1)
    body = haplotype(client, sample_reference(), 104, 106).json()
    assert body["version"] == 2
    assert body["sequence"] == "ACA"


def test_explicit_version_is_used_while_a_newer_one_is_current(client, app):
    upload(client, SAMPLE, base_version=0)
    upload(client, vcf([(105, "G", "C")]), base_version=1)

    seq = sample_reference()
    body = haplotype(client, seq, 104, 106, version=1).json()
    assert body["version"] == 1
    assert body["sequence"] == "ATA"  # G>T from version 1, not G>C

    body = haplotype(client, seq, 104, 106).json()
    assert body["version"] == 2 and body["sequence"] == "ACA"


def test_fixed_version_survives_concurrent_correction_and_reaper(tmp_path):
    """Instance 2 publishes a correction while instance 1 reads version 1."""
    from fastapi.testclient import TestClient

    from conftest import make_config
    from app.main import create_app

    app1 = create_app(make_config(tmp_path, instance_id="inst-1"))
    app2 = create_app(make_config(tmp_path, instance_id="inst-2"))
    with TestClient(app1) as c1, TestClient(app2) as c2:
        assert upload(c1, SAMPLE, base_version=0).status_code == 201

        # Instance 1 starts a fixed-version build: the pin holds version 1.
        store1 = app1.state.store
        pinned = store1.pin("inst-1", version=1)
        assert pinned[0] == 1
        ref_id = pinned[3]
        try:
            # The other instance corrects the list (v2) and reaps.
            assert upload(c2, vcf([(105, "G", "C")]), base_version=1).status_code == 201
            c2.app.state.reap_once()
            v1_dir = os.path.join(
                app1.state.cfg.data_dir, store1._connect().execute(
                    "SELECT path FROM versions WHERE version = 1"
                ).fetchone()["path"]
            )
            assert os.path.isdir(v1_dir), "pinned version must survive the reaper"

            body = haplotype(c1, sample_reference(), 104, 106, version=1).json()
            assert body["version"] == 1
            assert body["sequence"] == "ATA"
        finally:
            store1.unpin(ref_id)

        # Once unpinned, the reaper may reclaim v1 and the explicit version is gone.
        app1.state.reap_once()
        assert haplotype(c1, sample_reference(), 104, 106, version=1).status_code == 404
        # Current version still works.
        body = haplotype(c1, sample_reference(), 104, 106).json()
        assert body["version"] == 2 and body["sequence"] == "ACA"


def test_current_version_pinned_before_a_midrequest_publish(client, app, monkeypatch):
    """A publish that lands mid-request must not change the version used."""
    upload(client, SAMPLE, base_version=0)

    import app.main as main

    real_fetch = main.fetch_intersecting

    def publish_then_fetch(gz_path, chrom, start, end):
        # Another instance replaces the list after this request pinned v1.
        assert upload(client, vcf([(105, "G", "C")]), base_version=1).status_code == 201
        return real_fetch(gz_path, chrom, start, end)

    monkeypatch.setattr(main, "fetch_intersecting", publish_then_fetch)
    body = haplotype(client, sample_reference(), 104, 106).json()
    assert body["version"] == 1
    assert body["sequence"] == "ATA"
    assert client.get("/version").json()["current_version"] == 2
    # the next request sees v2
    monkeypatch.setattr(main, "fetch_intersecting", real_fetch)
    body = haplotype(client, sample_reference(), 104, 106).json()
    assert body["version"] == 2 and body["sequence"] == "ACA"


def test_pins_are_released_after_success_and_failure(client, app):
    upload(client, SAMPLE, base_version=0)
    assert haplotype(client, sample_reference(), 100, 110).status_code == 200
    assert haplotype(client, sample_reference(), 100, 400).status_code == 400
    assert haplotype(client, sample_reference(), 100, 110, version=9).status_code == 404
    with app.state.store._connect() as con:
        assert con.execute("SELECT COUNT(*) FROM refs").fetchone()[0] == 0


def test_failed_build_does_not_change_published_data(client):
    upload(client, SAMPLE, base_version=0)
    before = client.get("/version").json()
    variants_before = client.get(
        "/variants", params={"chrom": "chr1", "start": 100, "end": 110}
    ).json()

    for body, params in [
        (fasta(ref(300)), {}),                       # REF mismatch against SAMPLE
        (fasta(ref(300)), {"start": 1, "end": 400}),  # bad interval
        (b"not a fasta", {}),                        # bad FASTA
    ]:
        r = client.post(
            "/haplotype",
            params={"start": 100, "end": 110, **params},
            files={"file": ("reference.fa", body)},
        )
        assert r.status_code >= 400

    assert client.get("/version").json() == before
    assert client.get(
        "/variants", params={"chrom": "chr1", "start": 100, "end": 110}
    ).json() == variants_before


def test_whole_interval_from_position_one(client):
    upload(client, vcf([(1, "A", "C"), (2, "A", "G")]), base_version=0)
    sequence = ref(10)
    body = haplotype(client, sequence, 1, 3).json()
    assert body["sequence"] == "CGA"


def test_docs_and_openapi_expose_haplotype(client):
    assert client.get("/docs").status_code == 200
    schema = client.get("/openapi.json").json()
    assert "/haplotype" in schema["paths"]
    assert "post" in schema["paths"]["/haplotype"]
    params = {p["name"] for p in schema["paths"]["/haplotype"]["post"]["parameters"]}
    assert {"start", "end", "chrom", "version"} <= params


def test_original_endpoints_still_work(client):
    assert upload(client, SAMPLE, base_version=0).status_code == 201
    q = client.get("/variants", params={"chrom": "chr1", "start": 102, "end": 103})
    assert q.status_code == 200 and q.json()["count"] == 1
    assert client.get("/version").json()["current_version"] == 1
    assert client.get("/", follow_redirects=False).status_code in (307, 308)
