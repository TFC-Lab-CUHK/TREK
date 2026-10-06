#!/usr/bin/env python3
"""
annotate_apa.py — Standalone PA-cluster annotation script.

For each unique PA cluster (cluster_key) in a unified APA sites file:
  1. Determines the cluster genomic range: min/max of original_site_position
  2. Extracts ±50 bp sequence around mode_site_position (5'→3' gene-strand)
  3. Applies the lineage-specific PAS motif model assigned to the species
  4. Assigns PAS-L1..PAS-L3 from cluster-sample recurrence and motif support
  5. Writes one annotation row per cluster to an output TSV

Usage:
    python annotate_apa.py \\
        --unified-file data/Homo_sapiens/Homo_sapiens_unified_apa_sites.txt \\
        --fasta /path/to/GCF_000001405.40_GRCh38.p14_genomic.fna \\
        --output data/Homo_sapiens/Homo_sapiens_unified_apa.anno.txt

Output TSV columns:
    cluster_key         – unique gene-level cluster identifier
    unified_ID          – coordinate label from the unifier
    chromosome          – RefSeq chromosome (e.g. NC_000001.11)
    strand              – '+' or '-'
    mode_site_position  – representative cleavage position (1-based)
    cluster_start       – min(original_site_position) across all rows
    cluster_end         – max(original_site_position) across all rows
    sequence            – 101 bp window (±50 bp) centred on mode_site_position, 5'→3' on gene strand
    pas_motif           – best motif found (hexamer or 4-mer), or '' if none
    pas_position        – distance from cleavage site to motif start, or ''
    pas_type            – lineage-specific motif category or 'not_assessed'
    search_level        – motif element selected by the active model
    motif_model         – active lineage model or 'not_assessed'
    cleavage_context    – oriented dinucleotide at positions -1 and 0
    u_rich_fraction     – plant T fraction across -10..+15, otherwise ''
    n_samples           – distinct sample_id values (normalized `sample` names if sample_id is absent)
    motif_support       – 'primary' | 'auxiliary' | 'none' | 'not_assessed'
    pa_level            – PAS-L1..PAS-L3 evidence-convergence level
    pa_level_basis      – evidence combination underlying pa_level

Notes:
  - Motif tables are bundled under resources/ next to this script.
  - A reusable FASTA index is stored next to the output, or at --fasta-index.
  - Each run regenerates the annotation and atomically replaces the output.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import os
from pathlib import Path
import sys
import tempfile
import time
from typing import Dict, List, Optional

if __package__:
    from . import pas_motif_models as models
else:
    import pas_motif_models as models

CANONICAL_HEXAMERS = models.CANONICAL_HEXAMERS
VARIANT_HEXAMERS = models.VARIANT_HEXAMERS
FISSION_YEAST = models.FISSION_YEAST
MAMMALIAN = models.MAMMALIAN
annotate_context = models.annotate_context
load_spombe_auxiliary_motifs = models.load_spombe_auxiliary_motifs
resolve_species_model = models.resolve_species_model

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

_INDEX_CACHE: Dict[str, dict] = {}
_FASTA_HANDLE_CACHE: Dict[str, object] = {}
RESOURCE_DIR = Path(__file__).resolve().with_name("resources")
DEFAULT_SPECIES_MODEL_TABLE = RESOURCE_DIR / "species_motif_groups.tsv"
DEFAULT_SPOMBE_MOTIF_TABLE = RESOURCE_DIR / "spombe_auxiliary_motifs.tsv"
VERSION = "1.0.0"


def _published_file_mode() -> int:
    """Regular-file permissions under the current umask (tempfile defaults to 0600)."""
    current = os.umask(0)
    os.umask(current)
    return 0o666 & ~current


def _close_reference(fasta_path: str) -> None:
    handle = _FASTA_HANDLE_CACHE.pop(fasta_path, None)
    if handle is not None:
        handle.close()
    _INDEX_CACHE.pop(fasta_path, None)


def _build_fasta_index(fasta_path: str, index_path: Optional[str] = None) -> None:
    idx_path = Path(index_path or (fasta_path + ".fidx"))
    stat = Path(fasta_path).stat()
    identity = dict(path=str(Path(fasta_path).resolve()), size=stat.st_size, mtime_ns=stat.st_mtime_ns)
    _close_reference(fasta_path)
    if idx_path.exists():
        try:
            with idx_path.open() as handle:
                cached = json.load(handle)
            if cached.get("format_version") == 1 and cached.get("source") == identity:
                sequences = cached["sequences"]
                if not isinstance(sequences, dict) or not sequences:
                    raise ValueError("Empty FASTA index")
                for entry in sequences.values():
                    if not all(isinstance(entry.get(k), int) for k in ("offset", "length", "line_len", "line_bytes")):
                        raise ValueError("Invalid FASTA index")
                    if not (entry["offset"] >= 0 and entry["length"] > 0
                            and entry["line_bytes"] >= entry["line_len"] > 0):
                        raise ValueError("Invalid FASTA index bounds")
                _INDEX_CACHE[fasta_path] = sequences
                return
        except (ValueError, TypeError, KeyError, AttributeError):
            log.warning("Rebuilding incompatible FASTA index: %s", idx_path)
    log.info(f"Building FASTA index for {fasta_path} ...")
    t0 = time.time()
    index: Dict[str, dict] = {}
    current_chrom: Optional[str] = None
    first_base_offset: Optional[int] = None
    line_len: Optional[int] = None
    line_bytes: Optional[int] = None
    accumulated_len = 0
    previous_line_full = True

    def finish_sequence():
        if current_chrom is None:
            return
        if not accumulated_len:
            raise ValueError(f"Empty FASTA sequence: {current_chrom}")
        index[current_chrom] = {
            "offset": first_base_offset, "length": accumulated_len,
            "line_len": line_len, "line_bytes": line_bytes,
        }

    with open(fasta_path, "rb") as fh:
        if fh.read(2) == b"\x1f\x8b":
            raise ValueError("Provide an uncompressed genome FASTA")
        fh.seek(0)
        while True:
            pos = fh.tell()
            raw = fh.readline()
            if not raw:
                break
            line = raw.decode("ascii")
            if line.startswith(">"):
                finish_sequence()
                identifiers = line[1:].split()
                if not identifiers:
                    raise ValueError("FASTA header has no sequence identifier")
                current_chrom = identifiers[0]
                if current_chrom in index:
                    raise ValueError(f"Duplicate FASTA identifier: {current_chrom}")
                first_base_offset = None
                line_len = None
                line_bytes = None
                accumulated_len = 0
                previous_line_full = True
            else:
                bases = raw.rstrip(b"\r\n")
                if current_chrom is None:
                    raise ValueError("FASTA sequence must begin with a >header")
                if not bases or bases.translate(None, b"ACGTRYSWKMBDHVNacgtryswkmbdhvn"):
                    raise ValueError(f"Invalid or blank sequence line in FASTA: {current_chrom}")
                if not previous_line_full:
                    raise ValueError(f"FASTA has inconsistent line wrapping: {current_chrom}")
                if first_base_offset is None:
                    first_base_offset = pos
                    line_len = len(bases)
                    line_bytes = len(raw)
                if len(bases) > line_len:
                    raise ValueError(f"FASTA has inconsistent line wrapping: {current_chrom}")
                previous_line_full = len(bases) == line_len and len(raw) == line_bytes
                accumulated_len += len(bases)
    finish_sequence()
    if not index:
        raise ValueError("Genome FASTA contains no sequences")
    now = Path(fasta_path).stat()
    if (now.st_size, now.st_mtime_ns) != (stat.st_size, stat.st_mtime_ns):
        raise ValueError("Genome FASTA changed during indexing")
    idx_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", dir=idx_path.parent, prefix=".fasta-index-", delete=False) as fh:
            temporary = Path(fh.name)
            json.dump(dict(format_version=1, source=identity, sequences=index), fh)
        os.chmod(temporary, _published_file_mode())
        os.replace(temporary, idx_path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    _INDEX_CACHE[fasta_path] = index
    log.info(
        f"Index built: {len(index):,} sequences in {time.time() - t0:.1f}s → {idx_path}"
    )


def _load_fasta_index(fasta_path: str) -> dict:
    if fasta_path not in _INDEX_CACHE:
        _build_fasta_index(fasta_path)
        log.info(f"Loaded FASTA index ({len(_INDEX_CACHE[fasta_path]):,} chroms)")
    return _INDEX_CACHE[fasta_path]


def _fasta_entry(fasta_path: str, chrom: str) -> Optional[dict]:
    index = _load_fasta_index(fasta_path)
    alias = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
    return index.get(chrom) or index.get(alias)


def _fasta_handle(fasta_path: str):
    handle = _FASTA_HANDLE_CACHE.get(fasta_path)
    if handle is None or handle.closed:
        handle = open(fasta_path, "rb")
        _FASTA_HANDLE_CACHE[fasta_path] = handle
    return handle


def _fetch_seq(fasta_path: str, chrom: str, start: int, end: int) -> str:
    entry = _fasta_entry(fasta_path, chrom)
    if not entry:
        raise ValueError(f"Chromosome {chrom!r} is absent from the genome FASTA")
    chrom_len = int(entry.get("length") or 0)
    start = max(0, int(start))
    end = min(chrom_len, int(end))
    if end <= start:
        return ""
    start = max(0, start)
    offset, line_len, line_bytes = (
        entry["offset"],
        entry["line_len"],
        entry["line_bytes"],
    )
    if not line_len or not line_bytes:
        return ""
    byte_start = offset + (start // line_len) * line_bytes + (start % line_len)
    byte_end = (
        offset + ((end - 1) // line_len) * line_bytes + ((end - 1) % line_len) + 1
    )
    fh = _fasta_handle(fasta_path)
    fh.seek(byte_start)
    raw = fh.read(byte_end - byte_start)
    sequence = (
        raw.replace(b"\n", b"")
        .replace(b"\r", b"")
        .decode("ascii")
        .upper()[: end - start]
    )
    if len(sequence) != end - start:
        raise ValueError(f"Incomplete FASTA sequence read for {chrom}:{start}-{end}")
    return sequence


def _rev_comp(seq: str) -> str:
    return seq.translate(
        str.maketrans(
            "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
            "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
        )
    )[::-1]


FLANK = 50


SCAN_START = -150
SCAN_END = 61


def _oriented_context_bounds(
    pos_0: int,
    strand: str,
    chrom_length: int,
) -> tuple[int, int, int]:
    if strand == "+":
        genomic_start = max(0, pos_0 + SCAN_START)
        genomic_end = min(chrom_length, pos_0 + SCAN_END)
        actual_relative_start = genomic_start - pos_0
    elif strand == "-":
        genomic_start = max(0, pos_0 - SCAN_END + 1)
        genomic_end = min(chrom_length, pos_0 - SCAN_START + 1)
        actual_relative_start = pos_0 - (genomic_end - 1)
    else:
        raise ValueError(f"Invalid strand: {strand!r}")
    return genomic_start, genomic_end, actual_relative_start


def _extract_oriented_context(
    chrom_seq: str, pos_1: int, strand: str
) -> tuple[str, int]:
    if not chrom_seq:
        return "", SCAN_START
    if not 1 <= pos_1 <= len(chrom_seq):
        raise ValueError("PAS coordinate lies outside the reference sequence")
    pos_0 = pos_1 - 1
    start, end, relative_start = _oriented_context_bounds(
        pos_0, strand, len(chrom_seq)
    )
    sequence = chrom_seq[start:end]
    if strand == "-":
        sequence = _rev_comp(sequence)
    return sequence, relative_start


def _fetch_oriented_context(
    fasta_path: str, chrom: str, pos_1: int, strand: str
) -> tuple[str, int]:
    entry = _fasta_entry(fasta_path, chrom)
    if not entry:
        raise ValueError(f"Chromosome {chrom!r} is absent from the genome FASTA; check the genome build and sequence names")
    if not 1 <= pos_1 <= int(entry["length"]):
        raise ValueError(f"PAS coordinate {chrom}:{pos_1} lies outside the reference sequence")
    pos_0 = pos_1 - 1
    start, end, relative_start = _oriented_context_bounds(
        pos_0, strand, int(entry["length"])
    )
    sequence = _fetch_seq(fasta_path, chrom, start, end)
    if strand == "-":
        sequence = _rev_comp(sequence)
    return sequence, relative_start


def _display_sequence_from_context(context: str, context_start: int) -> str:
    relative_end = context_start + len(context)
    display_start = max(-FLANK, context_start)
    display_end = min(FLANK + 1, relative_end)
    if display_end <= display_start:
        return ""
    return context[
        display_start - context_start : display_end - context_start
    ]


def annotate_site(
    chrom_seq: str,
    pos_1: int,
    strand: str,
    motif_model: str = MAMMALIAN,
    spombe_auxiliary: Optional[Dict[str, List[str]]] = None,
) -> Dict:
    context, context_start = _extract_oriented_context(chrom_seq, pos_1, strand)
    return annotate_context(
        context, context_start, motif_model, spombe_auxiliary
    )


OUTPUT_COLUMNS = [
    "cluster_key",
    "unified_ID",
    "chromosome",
    "strand",
    "mode_site_position",
    "cluster_start",
    "cluster_end",
    "sequence",
    "pas_motif",
    "pas_position",
    "pas_type",
    "search_level",
    "motif_model",
    "cleavage_context",
    "u_rich_fraction",
    "n_samples",
    "motif_support",
    "pa_level",
    "pa_level_basis",
]


def assign_pas_level(n_samples: int, motif_support: str) -> tuple[str, str]:
    recurrent = n_samples >= 2
    primary_motif = motif_support == "primary"

    if recurrent and primary_motif:
        return "PAS-L1", "recurrence_and_motif"
    if recurrent:
        return "PAS-L2", "recurrence_only"
    if primary_motif:
        return "PAS-L2", "motif_only"
    return "PAS-L3", "baseline_call"


def _sample_identity_key(value: Optional[str]) -> str:
    return " ".join((value or "").replace("_", " ").split()).casefold()


def aggregate_clusters(unified_path: str) -> Dict[str, dict]:
    log.info(f"Reading unified file: {unified_path}")
    clusters: Dict[str, dict] = {}
    row_count = 0
    with open(unified_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"cluster_key", "chromosome", "strand", "original_site_position", "mode_site_position"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Unified file lacks required columns {sorted(missing)}: {unified_path}"
            )
        sample_column = next((c for c in ("sample_id", "_sample_uid", "sample") if c in reader.fieldnames), None)
        if sample_column is None:
            raise ValueError("Unified file requires a sample_id or sample column")
        log.info("Counting cluster recurrence using %s", sample_column)
        for row in reader:
            row_count += 1
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed unified-file data row {row_count}")
            cluster_key = row["cluster_key"].strip()
            if not cluster_key:
                raise ValueError(
                    f"Empty cluster_key at unified-file data row {row_count}"
                )
            orig_pos = int(row["original_site_position"])
            mode_pos = int(row["mode_site_position"])
            if orig_pos < 1 or mode_pos < 1 or row["strand"] not in ("+", "-") or not row["chromosome"].strip():
                raise ValueError(f"Invalid coordinate, strand or chromosome at data row {row_count}")
            if cluster_key not in clusters:
                clusters[cluster_key] = {
                    "cluster_key": cluster_key,
                    "unified_ID": row.get("unified_ID", ""),
                    "chromosome": row["chromosome"],
                    "strand": row["strand"],
                    "mode_site_position": mode_pos,
                    "cluster_start": orig_pos,
                    "cluster_end": orig_pos,
                    "samples": set(),
                }
            else:
                c = clusters[cluster_key]
                for column in ("chromosome", "strand", "unified_ID"):
                    if c[column] != row.get(column, ""):
                        raise ValueError(f"Cluster {cluster_key!r} has inconsistent {column}")
                if c["mode_site_position"] != mode_pos:
                    raise ValueError(f"Cluster {cluster_key!r} has inconsistent mode_site_position")
                c["cluster_start"] = min(c["cluster_start"], orig_pos)
                c["cluster_end"] = max(c["cluster_end"], orig_pos)

            sample_key = (row[sample_column].strip() if sample_column != "sample"
                          else _sample_identity_key(row["sample"]))
            if not sample_key:
                raise ValueError(f"Empty {sample_column} at unified-file data row {row_count}")
            clusters[cluster_key]["samples"].add(sample_key)

    for cluster in clusters.values():
        if not cluster["cluster_start"] <= cluster["mode_site_position"] <= cluster["cluster_end"]:
            raise ValueError(f"Representative coordinate outside cluster {cluster['cluster_key']!r}")
        cluster["n_samples"] = len(cluster.pop("samples"))

    log.info(f"  Rows: {row_count:,} | Unique clusters: {len(clusters):,}")
    return clusters


def annotate_clusters(
    clusters: Dict[str, dict],
    fasta_path: str,
    output_path: str,
    motif_model: str,
    spombe_auxiliary: Optional[Dict[str, List[str]]] = None,
    progress_every: int = 1000,
    index_path: Optional[str] = None,
) -> None:
    if progress_every < 1:
        raise ValueError("--progress-every must be positive")
    to_annotate = list(clusters.items())
    total = len(clusters)
    log.info(f"Clusters to annotate: {total:,}")
    output = Path(output_path)
    if output.resolve() == Path(fasta_path).resolve():
        raise ValueError("Annotation output must not overwrite the genome FASTA")
    output.parent.mkdir(parents=True, exist_ok=True)
    if total:
        _build_fasta_index(fasta_path, index_path)
    out_fh = tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
        dir=output.parent, prefix=f".{output.name}.", suffix=".tmp", delete=False)
    temporary = Path(out_fh.name)
    writer = csv.DictWriter(
        out_fh, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
    )
    try:
        writer.writeheader()
        for i, (cluster_key, info) in enumerate(to_annotate, 1):
            chrom = info["chromosome"]
            strand = info["strand"]
            mode_pos = info["mode_site_position"]

            context, context_start = _fetch_oriented_context(
                fasta_path, chrom, mode_pos, strand
            )
            pas = annotate_context(
                context,
                context_start,
                motif_model,
                spombe_auxiliary,
            )
            seq = _display_sequence_from_context(context, context_start)
            motif_support = pas["motif_support"]
            pa_level, pa_level_basis = assign_pas_level(
                info["n_samples"], motif_support
            )

            writer.writerow(
                {
                    "cluster_key": cluster_key,
                    "unified_ID": info["unified_ID"],
                    "chromosome": chrom,
                    "strand": strand,
                    "mode_site_position": mode_pos,
                    "cluster_start": info["cluster_start"],
                    "cluster_end": info["cluster_end"],
                    "sequence": seq,
                    "pas_motif": pas["motif"] or "",
                    "pas_position": pas["position"]
                    if pas["position"] is not None
                    else "",
                    "pas_type": pas["motif_type"],
                    "search_level": pas["search_level"],
                    "motif_model": motif_model.lower(),
                    "cleavage_context": pas["cleavage_context"],
                    "u_rich_fraction": pas["u_rich_fraction"]
                    if pas["u_rich_fraction"] is not None
                    else "",
                    "n_samples": info["n_samples"],
                    "motif_support": motif_support,
                    "pa_level": pa_level,
                    "pa_level_basis": pa_level_basis,
                }
            )

            if i % progress_every == 0 or i == total:
                log.info(f"  {i:,} / {total:,}  ({100 * i / total:.1f}%)")
        out_fh.close()
        os.chmod(temporary, _published_file_mode())
        os.replace(temporary, output)
    finally:
        out_fh.close()
        if temporary.exists():
            temporary.unlink()
        _close_reference(fasta_path)

    log.info(f"Annotation complete. Output: {output_path}")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="annotate_apa",
        description="Annotate APA clusters with PAS signals and ±50 bp sequence context.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument(
        "--unified-file",
        required=True,
        metavar="PATH",
        help="Path to *_unified_apa_sites.txt (TSV with header)",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        metavar="PATH",
        help="Path to reference genome FASTA (.fna / .fa / .fasta)",
    )
    parser.add_argument(
        "--output",
        required=True,
        metavar="PATH",
        help="Output TSV path (e.g. Homo_sapiens_apa_annotation.tsv)",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1000,
        metavar="N",
        help="Log progress every N clusters (default: 1000)",
    )
    parser.add_argument(
        "--species-key",
        metavar="NAME",
        help="Species folder key; inferred from the unified filename when omitted",
    )
    parser.add_argument(
        "--motif-model", type=str.upper,
        choices=sorted(models.IMPLEMENTED_MODELS | {models.NOT_ASSESSED}),
        help="Explicit model override for species not in the bundled table",
    )
    parser.add_argument(
        "--fasta-index", metavar="PATH",
        help="Reusable FASTA index path (default: a hidden .fidx file next to the output)",
    )
    parser.add_argument(
        "--species-model-table",
        default=str(DEFAULT_SPECIES_MODEL_TABLE),
        metavar="PATH",
        help="Species-to-motif-model TSV",
    )
    parser.add_argument(
        "--spombe-motif-table",
        default=str(DEFAULT_SPOMBE_MOTIF_TABLE),
        metavar="PATH",
        help="S. pombe auxiliary motif TSV",
    )
    return parser.parse_args(argv)


def infer_species_key(unified_path: str) -> str:
    name = Path(unified_path).name
    for suffix in ("_dp_segment_unified_apa_sites.txt", "_dp_segment_unified_apa_sites.tsv",
                   "_unified_apa_sites.txt", "_unified_apa_sites.tsv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(unified_path).parent.name


def run_annotation(args):
    for name in ("unified_file", "fasta", "output", "species_model_table", "spombe_motif_table"):
        setattr(args, name, str(Path(getattr(args, name)).expanduser().resolve()))
    for name in ("unified_file", "fasta"):
        if not Path(getattr(args, name)).is_file():
            raise ValueError(f"Input file not found: {getattr(args, name)}")
    if args.progress_every < 1:
        raise ValueError("--progress-every must be positive")
    protected = {args.unified_file, args.fasta, args.species_model_table, args.spombe_motif_table}
    if args.output in protected:
        raise ValueError("Annotation output must not overwrite an input or motif table")
    if args.fasta_index:
        index_path = str(Path(args.fasta_index).expanduser().resolve())
    else:
        source_key = hashlib.sha256(args.fasta.encode()).hexdigest()[:12]
        index_path = str(Path(args.output).parent / f".{Path(args.fasta).name}.{source_key}.fidx")
    if index_path in protected | {args.output}:
        raise ValueError("FASTA index must have a separate path from inputs and output")
    species_key = args.species_key or infer_species_key(args.unified_file)
    species_key = "_".join(species_key.split())
    if args.motif_model:
        motif_model = args.motif_model
    else:
        try:
            motif_model = resolve_species_model(species_key, Path(args.species_model_table))
        except ValueError as exc:
            raise ValueError(f"{exc} Specify --species-key, provide --species-model-table, "
                             "or choose an appropriate --motif-model explicitly.") from exc
    spombe_auxiliary = None
    if motif_model == FISSION_YEAST:
        spombe_auxiliary = load_spombe_auxiliary_motifs(
            Path(args.spombe_motif_table)
        )
    log.info(f"Species: {species_key} | PAS motif model: {motif_model}")

    clusters = aggregate_clusters(args.unified_file)
    annotate_clusters(
        clusters=clusters,
        fasta_path=args.fasta,
        output_path=args.output,
        motif_model=motif_model,
        spombe_auxiliary=spombe_auxiliary,
        progress_every=args.progress_every,
        index_path=index_path,
    )


def main(argv=None):
    args = parse_args(argv)
    try:
        run_annotation(args)
    except (OSError, ValueError, csv.Error) as exc:
        log.error("%s", exc)
        raise SystemExit(1) from None


if __name__ == "__main__":
    main()
