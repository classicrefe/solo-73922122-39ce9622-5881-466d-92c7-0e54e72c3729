"""FastAPI application: versioned VCF upload + indexed interval queries.

Layout of the shared data directory:

    <data_dir>/meta.sqlite3          metadata, current-version pointer, ref pins
    <data_dir>/versions/<uuid>/      published file sets (variants.vcf.gz + .tbi)
    <data_dir>/staging/<instance>/<uuid>/   in-progress builds

Publishing is atomic: the file set is completed in staging, renamed into
versions/, and only then the current-version pointer is flipped inside one
SQLite transaction. Readers pin a version before opening its files, so an
in-flight query keeps its (version, data, index) trio consistent while a
replacement is published; the reaper frees old files once no pin remains.
"""
from __future__ import annotations

import asyncio
import os
import shutil
import time
import uuid
from contextlib import asynccontextmanager
from dataclasses import dataclass

import pysam
from fastapi import FastAPI, File, HTTPException, Query, UploadFile
from fastapi.responses import RedirectResponse
from pydantic import BaseModel, Field

from .db import Conflict, Store
from .haplotype import (
    FastaError,
    HaplotypeError,
    apply_variants,
    check_interval,
    fetch_intersecting,
    parse_fasta,
    resolve_chrom,
    validate_records,
)
from .vcfbuild import GZ_NAME, VCFError, validate_and_build

DEFAULT_MAX_UPLOAD = 10 * 1024 * 1024  # 10 MiB


class HaplotypeResponse(BaseModel):
    """Alternate sequence for one interval of one pinned VCF version."""

    version: int = Field(description="VCF version the sequence was built from")
    chrom: str
    start: int = Field(description="1-based interval start (inclusive)")
    end: int = Field(description="1-based interval end (inclusive)")
    reference: str = Field(description="Reference bases of the interval")
    sequence: str = Field(
        description="Alternate sequence after applying the participating records"
    )
    count: int = Field(description="Number of participating records")
    records: list[str] = Field(
        description="VCF records that intersect the interval, sorted by position"
    )


@dataclass
class Config:
    data_dir: str
    instance_id: str
    max_upload_bytes: int = DEFAULT_MAX_UPLOAD
    reaper_interval: float = 2.0      # seconds between reaper runs (0 disables)
    heartbeat_interval: float = 2.0   # seconds between heartbeats (0 disables)
    dead_after: float = 15.0          # instance considered dead after this silence
    orphan_grace: float = 60.0        # min age before orphan dirs are swept


def _config_from_env() -> Config:
    e = os.environ
    return Config(
        data_dir=e.get("VCF_DATA_DIR", "/data"),
        instance_id=e.get("VCF_INSTANCE_ID") or uuid.uuid4().hex,
        max_upload_bytes=int(e.get("VCF_MAX_UPLOAD_BYTES", DEFAULT_MAX_UPLOAD)),
        reaper_interval=float(e.get("VCF_REAPER_INTERVAL", "2")),
        heartbeat_interval=float(e.get("VCF_HEARTBEAT_INTERVAL", "2")),
        dead_after=float(e.get("VCF_DEAD_AFTER", "15")),
        orphan_grace=float(e.get("VCF_ORPHAN_GRACE", "60")),
    )


def create_app(config: Config | None = None) -> FastAPI:
    cfg = config or _config_from_env()
    store = Store(os.path.join(cfg.data_dir, "meta.sqlite3"))
    versions_dir = os.path.join(cfg.data_dir, "versions")
    staging_base = os.path.join(cfg.data_dir, "staging")
    staging_root = os.path.join(staging_base, cfg.instance_id)
    os.makedirs(versions_dir, exist_ok=True)
    # Leftover builds of a previous incarnation of this instance can never be
    # published any more; drop them before serving.
    shutil.rmtree(staging_root, ignore_errors=True)
    os.makedirs(staging_root, exist_ok=True)

    def reap_once() -> None:
        now = time.time()
        # 1) DB: drop dead instances' pins, delete reclaimable version rows.
        for _version, rel in store.reap_metadata(now, cfg.dead_after):
            shutil.rmtree(os.path.join(cfg.data_dir, rel), ignore_errors=True)
        # 2) FS: sweep version dirs with no DB row (crashed publish/reaper).
        cutoff = now - cfg.orphan_grace
        known = store.all_version_paths()
        for name in os.listdir(versions_dir):
            path = os.path.join(versions_dir, name)
            rel = os.path.relpath(path, cfg.data_dir)
            try:
                old = os.path.getmtime(path) < cutoff
            except OSError:
                continue
            if rel not in known and old:
                shutil.rmtree(path, ignore_errors=True)
        # 3) FS: sweep staging dirs of dead instances.
        live = store.live_instances(now, cfg.dead_after)
        for inst in os.listdir(staging_base):
            path = os.path.join(staging_base, inst)
            try:
                old = os.path.getmtime(path) < cutoff
            except OSError:
                continue
            if inst not in live and old:
                shutil.rmtree(path, ignore_errors=True)

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        store.heartbeat(cfg.instance_id, time.time())
        stopped = asyncio.Event()

        async def heartbeat_loop():
            while not stopped.is_set():
                await asyncio.sleep(cfg.heartbeat_interval)
                try:
                    store.heartbeat(cfg.instance_id, time.time())
                except Exception:
                    pass

        async def reaper_loop():
            while not stopped.is_set():
                await asyncio.sleep(cfg.reaper_interval)
                try:
                    await asyncio.to_thread(reap_once)
                except Exception:
                    pass

        tasks = []
        if cfg.heartbeat_interval > 0:
            tasks.append(asyncio.create_task(heartbeat_loop()))
        if cfg.reaper_interval > 0:
            tasks.append(asyncio.create_task(reaper_loop()))
        try:
            yield
        finally:
            stopped.set()
            for t in tasks:
                t.cancel()
            await asyncio.gather(*tasks, return_exceptions=True)
            try:
                store.unregister(cfg.instance_id)
            except Exception:
                pass

    app = FastAPI(
        title="Versioned VCF variant store",
        version="1.0.0",
        lifespan=lifespan,
        description=(
            "Upload single-chromosome biallelic VCFs (max 10 MiB, no sample "
            "columns, REF/ALT in ACGT); records are sorted, BGZF-compressed and "
            "Tabix-indexed. Query 1-based closed intervals; every record whose "
            "reference span [POS, POS+len(REF)-1] intersects the interval is "
            "returned via the index. Uploads are optimistic-concurrency "
            "controlled through `base_version`. POST /haplotype overlays the "
            "records intersecting one interval (max 100,000 bases) onto an "
            "uploaded reference FASTA and returns the alternate sequence."
        ),
    )
    app.state.cfg = cfg
    app.state.store = store
    app.state.reap_once = reap_once

    @app.get("/", include_in_schema=False)
    def root():
        return RedirectResponse("/docs")

    @app.post("/vcf", status_code=201)
    async def upload_vcf(
        file: UploadFile = File(...),
        base_version: int = Query(..., ge=0, description="Version the client last saw (0 when empty)"),
    ):
        staging = os.path.join(staging_root, uuid.uuid4().hex)
        os.makedirs(staging)
        try:
            upload_path = os.path.join(staging, "upload.vcf")
            size = 0
            with open(upload_path, "wb") as out:
                while True:
                    chunk = await file.read(1 << 20)
                    if not chunk:
                        break
                    size += len(chunk)
                    if size > cfg.max_upload_bytes:
                        raise HTTPException(
                            413, f"upload exceeds the {cfg.max_upload_bytes}-byte limit"
                        )
                    out.write(chunk)
            try:
                chrom, n_records = await asyncio.to_thread(
                    validate_and_build, upload_path, staging
                )
            except VCFError as exc:
                raise HTTPException(400, f"invalid VCF: {exc}")
            os.remove(upload_path)

            # File set is complete; move it into place, then flip the pointer.
            final_rel = os.path.join("versions", uuid.uuid4().hex)
            os.rename(staging, os.path.join(cfg.data_dir, final_rel))
            try:
                new_version = await asyncio.to_thread(
                    store.publish, base_version, chrom, final_rel, n_records, time.time()
                )
            except Conflict as c:
                shutil.rmtree(os.path.join(cfg.data_dir, final_rel), ignore_errors=True)
                raise HTTPException(
                    409,
                    f"base_version {base_version} does not match current version "
                    f"{c.current}; re-read and retry",
                )
        finally:
            if os.path.isdir(staging):
                shutil.rmtree(staging, ignore_errors=True)
        return {"version": new_version, "chrom": chrom, "n_records": n_records}

    @app.get("/variants")
    def get_variants(
        chrom: str = Query(..., description="Chromosome name as in the VCF"),
        start: int = Query(..., ge=1, description="1-based interval start (inclusive)"),
        end: int = Query(..., ge=1, description="1-based interval end (inclusive)"),
    ):
        if end < start:
            raise HTTPException(400, "end must be >= start")
        pinned = store.pin(cfg.instance_id, time.time())
        if pinned is None:
            raise HTTPException(404, "no VCF has been published yet")
        version, vchrom, rel, ref_id = pinned
        try:
            if chrom != vchrom:
                raise HTTPException(
                    400, f"the published dataset contains chromosome '{vchrom}' only"
                )
            gz_path = os.path.join(cfg.data_dir, rel, GZ_NAME)
            try:
                tb = pysam.TabixFile(gz_path)
                try:
                    # 1-based closed [start, end] -> 0-based half-open fetch.
                    lines = tb.fetch(vchrom, start - 1, end)
                    records = []
                    for ln in lines:
                        f = ln.split("\t", 4)
                        pos = int(f[1])
                        # Exact span check on top of the index iterator.
                        if pos <= end and pos + len(f[3]) - 1 >= start:
                            records.append(ln)
                finally:
                    tb.close()
            except HTTPException:
                raise
            except Exception as exc:
                raise HTTPException(500, f"indexed fetch failed: {exc}")
            return {
                "version": version,
                "chrom": vchrom,
                "start": start,
                "end": end,
                "count": len(records),
                "records": records,
            }
        finally:
            store.unpin(ref_id)

    @app.post("/haplotype", response_model=HaplotypeResponse)
    async def build_haplotype(
        file: UploadFile = File(..., description="Single-chromosome ACGT reference FASTA"),
        chrom: str | None = Query(
            None,
            description="Chromosome name; must match both the FASTA record and the dataset "
            "(may be omitted when the FASTA header names it)",
        ),
        start: int = Query(..., ge=1, description="1-based interval start (inclusive)"),
        end: int = Query(..., ge=1, description="1-based interval end (inclusive)"),
        version: int | None = Query(
            None,
            ge=1,
            description="VCF version to apply; defaults to the version current when the "
            "request starts",
        ),
    ):
        """Build the alternate DNA sequence for one reference interval.

        The whole reference FASTA is uploaded (at most 10 MiB) and validated,
        then every published record whose reference span intersects the
        interval is overlaid as one haplotype.  The VCF version is pinned
        before the upload is read and used for the whole request: a
        replacement published by another instance mid-request cannot change
        which list this build applies.
        """
        # Pin first: whichever version this request starts on is the one it
        # uses end-to-end, even if the other instance publishes meanwhile.
        pinned = store.pin(cfg.instance_id, time.time(), version)
        if pinned is None:
            if version is not None:
                raise HTTPException(
                    404,
                    f"version {version} is not available (never published, or "
                    "already reclaimed)",
                )
            raise HTTPException(404, "no VCF has been published yet")
        version_id, vchrom, rel, ref_id = pinned
        try:
            data = bytearray()
            while True:
                chunk = await file.read(1 << 20)
                if not chunk:
                    break
                data.extend(chunk)
                if len(data) > cfg.max_upload_bytes:
                    raise HTTPException(
                        413, f"upload exceeds the {cfg.max_upload_bytes}-byte limit"
                    )
            try:
                fasta = await asyncio.to_thread(parse_fasta, bytes(data))
            except FastaError as exc:
                raise HTTPException(400, f"invalid FASTA: {exc}")
            try:
                check_interval(start, end, len(fasta.sequence))
            except HaplotypeError as exc:
                raise HTTPException(400, str(exc))

            try:
                target = resolve_chrom(chrom, fasta.name, vchrom)
                gz_path = os.path.join(cfg.data_dir, rel, GZ_NAME)
                lines = await asyncio.to_thread(
                    fetch_intersecting, gz_path, target, start, end
                )
                variants = await asyncio.to_thread(
                    validate_records, lines, fasta.sequence
                )
                sequence = apply_variants(fasta.sequence, start, end, variants)
            except HaplotypeError as exc:
                raise HTTPException(400, str(exc))
            except HTTPException:
                raise
            except Exception as exc:
                raise HTTPException(500, f"indexed fetch failed: {exc}")
        finally:
            # Synchronous on purpose: unpinning must also happen when the
            # request task is cancelled (client disconnect).  It is one tiny
            # indexed DELETE, so blocking the loop here is harmless.
            store.unpin(ref_id)

        return HaplotypeResponse(
            version=version_id,
            chrom=target,
            start=start,
            end=end,
            reference=fasta.sequence[start - 1:end],
            sequence=sequence,
            count=len(variants),
            records=[v.line for v in variants],
        )

    @app.get("/version")
    def get_version():
        info = store.version_info()
        if info is None:
            return {"published": False, "current_version": 0}
        return {
            "published": True,
            "current_version": info["version"],
            "chrom": info["chrom"],
            "n_records": info["n_records"],
            "created_at": info["created_at"],
        }

    return app


# Run with: uvicorn --factory app.main:create_app
