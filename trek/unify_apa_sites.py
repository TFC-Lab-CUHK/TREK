#!/usr/bin/env python3
"""Unify single-sample TREK PAS calls within one species using DP segmentation.

Standalone usage (Python with NumPy and pandas):
    python unify_apa_sites.py /path/to/Species --output-dir /path/to/output

All calculations are contained in this file. Defaults: inclusive cluster span
<=40 nt, adjacent gaps <=24 nt, and an optional strong single-observation rescue
at >=50 reads and >=0.25 pre-filter relative abundance. Coordinates remain
1-based and inclusive.
"""

from __future__ import annotations

import argparse
import bisect
import logging
import math
import os
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

try:
    import numpy as np
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit(
        f"Missing Python dependency: {exc.name}. Install with: "
        "python -m pip install numpy pandas"
    ) from exc


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_BIOTYPES = ("mRNA", "lncRNA", "lnc_RNA", "protein_coding", "")
EXPECTED_COLUMNS = {
    "transcript_id",
    "gene_id",
    "gene_name",
    "chromosome",
    "strand",
    "site_position",
    "site_count",
    "site_abundance",
    "transcript_biotype",
}

# Human and mouse RefSeq mitochondrial accessions. Other species can be kept by
# passing --keep-mitochondrial if their accession is not listed here.
MITOCHONDRIAL_REFSEQ_ACCESSIONS = {
    "NC_012920.1",  # Homo sapiens
    "NC_005089.1",  # Mus musculus
}

GENE_KEY = ["chromosome", "effective_gene_id", "strand"]
JOIN_KEY = GENE_KEY + ["site_position"]
SAMPLE_UID_COLUMN = "sample_id"
AUTO_RESCUE_CUTOFF = "auto"
VERSION = "1.0.0"
CLUSTER_COLUMNS = [
    "chromosome", "effective_gene_id", "strand", "site_position", "cluster_key",
    "unified_ID", "cluster_start", "cluster_end", "cluster_width",
    "n_unique_positions", "mode_site_position", "mode_score", "mode_local_n_samples",
    "mode_local_total_site_count", "mode_exact_n_samples", "mode_exact_total_site_count",
]
DETAILED_COLUMNS = [
    "transcript_id", "gene_id", "gene_name", "effective_gene_id", "chromosome",
    "strand", "cluster_key", "unified_ID", "original_site_position", "cluster_start",
    "cluster_end", "cluster_width", "mode_site_position", "mode_score",
    "mode_local_n_samples", "mode_local_total_site_count", "mode_exact_n_samples",
    "mode_exact_total_site_count", "n_unique_positions", "site_count",
    "transcript_biotype", "sample_attribute", "sample", SAMPLE_UID_COLUMN,
]
SUMMARY_COLUMNS = [
    "cluster_key", "unified_ID", "gene_id", "gene_name", "effective_gene_id",
    "chromosome", "strand", "cluster_start", "cluster_end", "cluster_width",
    "mode_site_position", "mode_score", "mode_local_n_samples",
    "mode_local_total_site_count", "mode_exact_n_samples", "mode_exact_total_site_count",
    "n_unique_positions", "weighted_mean_position", "weighted_position_sd",
    "position_iqr", "n_transcripts", "transcript_ids", "transcript_cluster_counts",
    "total_site_count", "n_samples", "n_tissues", "n_cell_cultures", "sample_list",
    "transcript_biotype", "median_cluster_relative_abundance", "filter_pass_rule",
    "filter_reason",
]


def _parse_rescue_cutoff(value: str) -> int | str:
    """Parse an explicit rescue cutoff or the adaptive cutoff sentinel."""
    if value == AUTO_RESCUE_CUTOFF:
        return AUTO_RESCUE_CUTOFF

    try:
        cutoff = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "strong single-observation cutoff must be auto or a non-negative integer"
        ) from exc

    if cutoff < 0 or str(cutoff) != str(value):
        raise ValueError(
            "strong single-observation cutoff must be auto or a non-negative integer"
        )
    return cutoff


def resolve_strong_single_observation_count(
    value: int | str, effective_sample_files: int
) -> int:
    """Resolve an explicit or sample-count-scaled strong-observation cutoff."""
    if value != AUTO_RESCUE_CUTOFF:
        return int(value)
    if effective_sample_files < 1:
        raise ValueError(
            "adaptive rescue cutoff requires at least one effective sample file"
        )
    return min(effective_sample_files * 10, 50)


# ---------------------------------------------------------------------------
# File discovery and input filtering
# ---------------------------------------------------------------------------


def _sample_name_from_path(fpath: Path, base_dir: Path) -> str:
    sample_dir = fpath.parent
    if sample_dir.name.startswith("trek_rerun") and sample_dir.parent != base_dir:
        sample_dir = sample_dir.parent
    return sample_dir.name


def discover_apa_files(base_dir: Path, pattern: str = "*apa_sites.txt") -> List[dict]:
    """Find input *apa_sites.txt files anywhere under one species directory."""
    files = sorted(p for p in base_dir.rglob(pattern) if p.is_file())
    if not files:
        raise ValueError(f"No input files matching {pattern!r} under {base_dir}")

    discovered: List[dict] = []
    seen = set()
    for fpath in files:
        resolved = fpath.resolve()
        if resolved in seen:
            log.warning("Skipping duplicate path to the same input file: %s", fpath)
            continue
        seen.add(resolved)
        parts = fpath.relative_to(base_dir).parts
        sample_attr = parts[0] if len(parts) >= 2 else "sample"

        discovered.append(
            {
                "path": fpath,
                "sample_attribute": sample_attr,
                "sample": _sample_name_from_path(fpath, base_dir),
                "sample_uid": fpath.relative_to(base_dir).as_posix(),
            }
        )

    log.info("Discovered %d APA site files under %s", len(discovered), base_dir)
    return discovered


def _parse_biotypes(value: str) -> set:
    biotypes = {item.strip() for item in value.split(",")}
    return {"lncRNA" if item == "lnc_RNA" else item for item in biotypes}


def read_and_filter(
    file_infos: List[dict],
    biotypes: set,
    keep_scaffolds: bool,
    keep_mitochondrial: bool,
    exclude_predicted_transcripts: bool,
) -> pd.DataFrame:
    """Read all APA files, concatenate, and apply batch-safe input filters."""
    frames: List[pd.DataFrame] = []
    for info in file_infos:
        fpath = info["path"]
        try:
            header = pd.read_csv(fpath, sep="\t", nrows=0).columns
            if {"unified_ID", "original_site_position"}.issubset(header):
                log.info("Skipping previously unified output: %s", fpath)
                continue
            df = pd.read_csv(
                fpath,
                sep="\t",
                dtype={
                    "transcript_id": str,
                    "gene_id": str,
                    "gene_name": str,
                    "chromosome": str,
                    "strand": str,
                    "transcript_biotype": str,
                },
            )
        except (OSError, ValueError, pd.errors.ParserError) as exc:
            raise ValueError(f"Cannot read PAS input {fpath}: {exc}") from exc

        missing = EXPECTED_COLUMNS - set(df.columns)
        if missing:
            raise ValueError(
                f"PAS input {fpath} is missing columns: {', '.join(sorted(missing))}"
            )

        df["sample_attribute"] = info["sample_attribute"]
        df["sample"] = info["sample"]
        df[SAMPLE_UID_COLUMN] = info["sample_uid"]
        for column in ("site_position", "site_count"):
            values = pd.to_numeric(df[column], errors="coerce")
            invalid = ~np.isfinite(values) | (values <= 0) | (values % 1 != 0)
            if invalid.any():
                row = int(np.flatnonzero(invalid.to_numpy())[0]) + 2
                raise ValueError(
                    f"{fpath}:{row}: {column} must be a finite positive integer"
                )
        invalid = (
            ~df["strand"].isin(["+", "-"])
            | df["chromosome"].fillna("").str.strip().eq("")
            | df["transcript_id"].fillna("").str.strip().eq("")
        )
        if invalid.any():
            row = int(np.flatnonzero(invalid.to_numpy())[0]) + 2
            raise ValueError(f"{fpath}:{row}: invalid strand, chromosome or transcript_id")
        frames.append(df)

    if not frames:
        raise ValueError("No single-sample PAS inputs found; unified outputs are not inputs")

    log.info("Successfully loaded %d APA site files", len(frames))
    combined = pd.concat(frames, ignore_index=True)
    log.info("Total rows loaded: %d", len(combined))

    combined["site_position"] = pd.to_numeric(
        combined["site_position"], errors="coerce"
    )
    combined["site_count"] = pd.to_numeric(combined["site_count"], errors="coerce")
    combined["site_abundance"] = pd.to_numeric(
        combined["site_abundance"], errors="coerce"
    )

    before = len(combined)
    combined = combined.dropna(subset=["site_position", "site_count"]).copy()
    combined["site_position"] = combined["site_position"].astype(int)
    combined["site_count"] = combined["site_count"].astype(float)
    log.info("After numeric cleanup: %d rows (removed %d)", len(combined), before - len(combined))

    combined["transcript_biotype"] = combined["transcript_biotype"].replace(
        {"lnc_RNA": "lncRNA"}
    )
    combined["transcript_biotype"] = combined["transcript_biotype"].fillna("")
    combined["gene_id"] = combined["gene_id"].fillna("").astype(str).str.strip()
    combined["gene_name"] = combined["gene_name"].fillna("").astype(str).str.strip()
    combined["effective_gene_id"] = combined["gene_id"].where(
        combined["gene_id"] != "", combined["gene_name"]
    )
    before = len(combined)
    combined = combined[
        combined["effective_gene_id"].astype(str).str.strip() != ""
    ].copy()
    log.info(
        "After effective gene id cleanup: %d rows (removed %d)",
        len(combined),
        before - len(combined),
    )

    before = len(combined)
    combined = combined[combined["transcript_biotype"].isin(biotypes)].copy()
    log.info(
        "After biotype filter (%s): %d rows (removed %d)",
        ",".join(sorted(biotypes)),
        len(combined),
        before - len(combined),
    )

    if exclude_predicted_transcripts:
        before = len(combined)
        combined = combined[
            ~combined["transcript_id"].str.startswith(("XM_", "XR_"), na=False)
        ].copy()
        log.info(
            "After predicted transcript filter (XM_/XR_ removed): %d rows (removed %d)",
            len(combined),
            before - len(combined),
        )

    if not keep_scaffolds:
        before = len(combined)
        combined = combined[combined["chromosome"].str.startswith("NC_", na=False)].copy()
        log.info(
            "After primary RefSeq accession filter (kept NC_*): %d rows (removed %d)",
            len(combined),
            before - len(combined),
        )

    if not keep_mitochondrial:
        before = len(combined)
        combined = combined[
            ~combined["chromosome"].isin(MITOCHONDRIAL_REFSEQ_ACCESSIONS)
        ].copy()
        log.info(
            "After mitochondrial filter: %d rows (removed %d)",
            len(combined),
            before - len(combined),
        )

    if combined.empty:
        log.warning("No PAS rows remain after input filtering; writing empty output tables")

    contexts = combined[["sample_attribute", "sample", SAMPLE_UID_COLUMN]].drop_duplicates()
    for (attr, sample), group in contexts.groupby(["sample_attribute", "sample"]):
        if len(group) > 1:
            log.warning(
                "%s/%s contains %d input files, counted independently. "
                "Use --input-pattern to exclude old/rerun copies of the same sample.",
                attr, sample, len(group),
            )

    sample_counts = (
        combined.groupby([SAMPLE_UID_COLUMN, "sample"], dropna=False)
        .size()
        .sort_index()
    )
    for (sample_uid, sample_name), count in sample_counts.items():
        log.info(
            "  Sample %-24s (%s): %d records",
            sample_name,
            sample_uid,
            count,
        )

    return combined


# ---------------------------------------------------------------------------
# Mode-seeded clustering
# ---------------------------------------------------------------------------


def _position_evidence(group: pd.DataFrame) -> Tuple[Dict[int, set], Dict[int, float]]:
    sample_sets: Dict[int, set] = {}
    total_counts: Dict[int, float] = {}
    for pos, sub in group.groupby("site_position", sort=True):
        pos_int = int(pos)
        sample_sets[pos_int] = set(
            sub[SAMPLE_UID_COLUMN].dropna().astype(str)
        )
        total_counts[pos_int] = float(sub["site_count"].sum())
    return sample_sets, total_counts


def _interval_stats(
    positions: Sequence[int],
    sample_sets: Dict[int, set],
    total_counts: Dict[int, float],
) -> Tuple[int, float]:
    samples = set()
    total = 0.0
    for pos in positions:
        samples.update(sample_sets[pos])
        total += total_counts[pos]
    return len(samples), total


def _local_stats(
    pos: int,
    positions: Sequence[int],
    sample_sets: Dict[int, set],
    total_counts: Dict[int, float],
    window: int,
) -> Tuple[int, float]:
    left = bisect.bisect_left(positions, pos - window)
    right = bisect.bisect_right(positions, pos + window)
    return _interval_stats(positions[left:right], sample_sets, total_counts)






def _choose_mode(
    cluster_positions: Sequence[int],
    sample_sets: Dict[int, set],
    total_counts: Dict[int, float],
    mode_support_window: int,
    sample_weight: float,
    count_weight: float,
    mode_sample_weight: float,
    mode_count_weight: float,
) -> dict:
    local: Dict[int, Tuple[int, float]] = {
        pos: _local_stats(
            pos, cluster_positions, sample_sets, total_counts, mode_support_window
        )
        for pos in cluster_positions
    }
    max_samples = max((value[0] for value in local.values()), default=1) or 1
    max_count = max((value[1] for value in local.values()), default=1.0) or 1.0
    max_exact_samples = max((len(sample_sets[pos]) for pos in cluster_positions), default=1) or 1
    max_exact_count = max((total_counts[pos] for pos in cluster_positions), default=1.0) or 1.0

    local_scored: Dict[int, Tuple[float, int, float]] = {}
    for pos in cluster_positions:
        local_samples, local_count = local[pos]
        score = (
            sample_weight * local_samples / max_samples
            + count_weight * local_count / max_count
        )
        local_scored[pos] = (score, local_samples, local_count)

    cluster_weights = np.array([total_counts[pos] for pos in cluster_positions], dtype=float)
    if cluster_weights.sum() > 0:
        cluster_center = float(np.average(np.array(cluster_positions), weights=cluster_weights))
    else:
        cluster_center = (cluster_positions[0] + cluster_positions[-1]) / 2.0

    best = None
    best_key = None
    for pos in cluster_positions:
        local_score, local_samples, local_count = local_scored[pos]
        exact_samples = len(sample_sets[pos])
        exact_count = total_counts[pos]
        mode_score = (
            mode_sample_weight * exact_samples / max_exact_samples
            + mode_count_weight * exact_count / max_exact_count
        )
        key = (
            mode_score,
            exact_samples,
            exact_count,
            local_score,
            local_samples,
            local_count,
            -abs(pos - cluster_center),
            -pos,
        )
        if best_key is None or key > best_key:
            best_key = key
            best = {
                "mode_site_position": pos,
                "mode_score": round(mode_score, 4),
                "mode_local_n_samples": local_samples,
                "mode_local_total_site_count": local_count,
                "mode_exact_n_samples": exact_samples,
                "mode_exact_total_site_count": exact_count,
            }

    assert best is not None
    return best






def _position_dp_weights(
    positions: Sequence[int],
    sample_sets: Dict[int, set],
    total_counts: Dict[int, float],
    sample_weight: float,
    count_weight: float,
) -> Dict[int, float]:
    """Build log-scaled evidence weights for coordinate dispersion."""
    max_samples = max((len(sample_sets[pos]) for pos in positions), default=1) or 1
    max_log_count = (
        max((math.log1p(total_counts[pos]) for pos in positions), default=1.0) or 1.0
    )

    weights: Dict[int, float] = {}
    for pos in positions:
        sample_score = len(sample_sets[pos]) / max_samples
        count_score = math.log1p(total_counts[pos]) / max_log_count
        weight = sample_weight * sample_score + count_weight * count_score
        weights[pos] = max(weight, 1e-6)
    return weights


def _segment_dispersion_cost(
    segment_positions: Sequence[int],
    dp_weights: Dict[int, float],
    compactness_weight: float,
) -> float:
    if len(segment_positions) <= 1 or compactness_weight <= 0:
        return 0.0

    values = np.array(segment_positions, dtype=float)
    weights = np.array([dp_weights[pos] for pos in segment_positions], dtype=float)
    center = float(np.average(values, weights=weights))
    variance = float(np.average((values - center) ** 2, weights=weights))
    return compactness_weight * variance


def _observation_counts_for_positions(
    group: pd.DataFrame,
    positions: Sequence[int],
) -> pd.Series:
    sub = group[group["site_position"].isin(positions)]
    if sub.empty:
        return pd.Series(dtype=float)
    return (
        sub.groupby(["transcript_id", SAMPLE_UID_COLUMN], dropna=False)["site_count"]
        .sum()
        .sort_values(ascending=False)
    )


def _cluster_support_stats(
    group: pd.DataFrame,
    positions: Sequence[int],
    transcript_sample_totals: pd.Series,
    min_observation_count: int,
    strong_single_observation_count: int,
    strong_single_min_relative_abundance: float,
) -> dict:
    observation_counts = _observation_counts_for_positions(group, positions)
    if observation_counts.empty:
        return {
            "supported_observations": 0,
            "max_observation_count": 0.0,
            "total_count": 0.0,
            "is_strong": False,
            "is_weak_singleton": len(positions) == 1,
        }

    max_observation = float(observation_counts.iloc[0])
    total_count = float(observation_counts.sum())
    supported = int((observation_counts >= min_observation_count).sum())
    observation_relative = (
        observation_counts
        / transcript_sample_totals.reindex(observation_counts.index).replace(0, np.nan)
    ).fillna(0.0)
    count_eligible_relative = observation_relative[
        observation_counts >= strong_single_observation_count
    ]
    strong_single = (
        strong_single_observation_count > 0
        and max_observation >= strong_single_observation_count
        and not count_eligible_relative.empty
        and float(count_eligible_relative.max())
        >= strong_single_min_relative_abundance
    )

    return {
        "supported_observations": supported,
        "max_observation_count": max_observation,
        "total_count": total_count,
        "is_strong": supported >= 2 or strong_single,
        "is_weak_singleton": (
            len(positions) == 1
            and supported < 1
            and max_observation < min_observation_count
        ),
    }


def _valid_positions(
    positions: Sequence[int],
    max_cluster_width: int,
    max_internal_gap: int,
) -> bool:
    if not positions:
        return False
    if positions[-1] - positions[0] + 1 > max_cluster_width:
        return False
    if len(positions) == 1:
        return True
    return all(
        right - left <= max_internal_gap
        for left, right in zip(positions, positions[1:])
    )


def _valid_segment(
    positions: Sequence[int],
    start_idx: int,
    end_idx: int,
    max_cluster_width: int,
    max_internal_gap: int,
) -> bool:
    if positions[end_idx] - positions[start_idx] + 1 > max_cluster_width:
        return False
    if end_idx == start_idx:
        return True
    return all(
        right - left <= max_internal_gap
        for left, right in zip(
            positions[start_idx:end_idx],
            positions[start_idx + 1 : end_idx + 1],
        )
    )


def _dp_partition_positions(
    positions: Sequence[int],
    dp_weights: Dict[int, float],
    max_cluster_width: int,
    max_internal_gap: int,
    cluster_penalty: float,
    compactness_weight: float,
) -> List[List[int]]:
    """Globally partition one gene's sorted positions into compact clusters."""
    n_positions = len(positions)
    if n_positions == 0:
        return []

    dp_cost = [math.inf] * (n_positions + 1)
    dp_clusters = [10**12] * (n_positions + 1)
    previous = [-1] * (n_positions + 1)
    dp_cost[0] = 0.0
    dp_clusters[0] = 0

    for start_idx in range(n_positions):
        if not math.isfinite(dp_cost[start_idx]):
            continue

        for end_idx in range(start_idx, n_positions):
            if positions[end_idx] - positions[start_idx] + 1 > max_cluster_width:
                break
            if (
                end_idx > start_idx
                and positions[end_idx] - positions[end_idx - 1] > max_internal_gap
            ):
                break
            if not _valid_segment(
                positions,
                start_idx,
                end_idx,
                max_cluster_width,
                max_internal_gap,
            ):
                continue

            segment_positions = positions[start_idx : end_idx + 1]
            segment_cost = _segment_dispersion_cost(
                segment_positions,
                dp_weights,
                compactness_weight,
            )
            candidate_cost = dp_cost[start_idx] + cluster_penalty + segment_cost
            candidate_clusters = dp_clusters[start_idx] + 1
            target = end_idx + 1

            if (
                candidate_cost < dp_cost[target] - 1e-9
                or (
                    abs(candidate_cost - dp_cost[target]) <= 1e-9
                    and candidate_clusters < dp_clusters[target]
                )
            ):
                dp_cost[target] = candidate_cost
                dp_clusters[target] = candidate_clusters
                previous[target] = start_idx

    clusters: List[List[int]] = []
    cursor = n_positions
    while cursor > 0:
        start_idx = previous[cursor]
        if start_idx < 0:
            start_idx = cursor - 1
        clusters.append(list(positions[start_idx:cursor]))
        cursor = start_idx

    clusters.reverse()
    return clusters


def _absorb_weak_singletons(
    clusters: List[List[int]],
    group: pd.DataFrame,
    max_cluster_width: int,
    max_internal_gap: int,
    min_observation_count: int,
    strong_single_observation_count: int,
    strong_single_min_relative_abundance: float,
    enabled: bool,
) -> List[List[int]]:
    if not enabled or len(clusters) < 2:
        return clusters

    transcript_sample_totals = group.groupby(
        ["transcript_id", SAMPLE_UID_COLUMN], dropna=False
    )["site_count"].sum()
    merged = [list(cluster) for cluster in clusters]
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(merged):
            cluster = merged[i]
            stats = _cluster_support_stats(
                group,
                cluster,
                transcript_sample_totals,
                min_observation_count,
                strong_single_observation_count,
                strong_single_min_relative_abundance,
            )
            if not stats["is_weak_singleton"]:
                i += 1
                continue

            candidates = []
            for neighbor_idx in (i - 1, i + 1):
                if neighbor_idx < 0 or neighbor_idx >= len(merged):
                    continue
                neighbor = merged[neighbor_idx]
                neighbor_stats = _cluster_support_stats(
                    group,
                    neighbor,
                    transcript_sample_totals,
                    min_observation_count,
                    strong_single_observation_count,
                    strong_single_min_relative_abundance,
                )
                combined = sorted(cluster + neighbor)
                if not neighbor_stats["is_strong"]:
                    continue
                if not _valid_positions(combined, max_cluster_width, max_internal_gap):
                    continue
                distance = min(
                    abs(cluster[0] - neighbor[0]),
                    abs(cluster[0] - neighbor[-1]),
                )
                candidates.append(
                    (
                        distance,
                        -neighbor_stats["supported_observations"],
                        -neighbor_stats["max_observation_count"],
                        neighbor_idx,
                        combined,
                    )
                )

            if not candidates:
                i += 1
                continue

            _distance, _supported, _max_count, neighbor_idx, combined = min(
                candidates
            )
            keep_idx = min(i, neighbor_idx)
            drop_idx = max(i, neighbor_idx)
            merged[keep_idx] = combined
            del merged[drop_idx]
            changed = True
            i = max(0, keep_idx - 1)

    return merged


def _cluster_one_gene_dp(
    key: Tuple[str, str, str],
    group: pd.DataFrame,
    max_cluster_width: int,
    max_internal_gap: int,
    mode_support_window: int,
    sample_weight: float,
    count_weight: float,
    mode_sample_weight: float,
    mode_count_weight: float,
    dp_cluster_penalty: float,
    dp_compactness_weight: float,
    min_observation_count: int,
    strong_single_observation_count: int,
    strong_single_min_relative_abundance: float,
    absorb_weak_singletons: bool,
) -> List[dict]:
    chrom, effective_gene_id, strand = key
    sample_sets, total_counts = _position_evidence(group)
    positions = sorted(sample_sets)
    dp_weights = _position_dp_weights(
        positions,
        sample_sets,
        total_counts,
        sample_weight,
        count_weight,
    )
    clusters = _dp_partition_positions(
        positions,
        dp_weights,
        max_cluster_width,
        max_internal_gap,
        dp_cluster_penalty,
        dp_compactness_weight,
    )
    clusters = _absorb_weak_singletons(
        clusters,
        group,
        max_cluster_width,
        max_internal_gap,
        min_observation_count,
        strong_single_observation_count,
        strong_single_min_relative_abundance,
        absorb_weak_singletons,
    )

    cluster_rows: List[dict] = []
    for cluster_number, cluster_positions in enumerate(clusters):
        cluster_start = cluster_positions[0]
        cluster_end = cluster_positions[-1]
        cluster_width = cluster_end - cluster_start + 1
        mode = _choose_mode(
            cluster_positions,
            sample_sets,
            total_counts,
            mode_support_window,
            sample_weight,
            count_weight,
            mode_sample_weight,
            mode_count_weight,
        )

        cluster_key = f"{chrom}|{effective_gene_id}|{strand}|{cluster_number}"
        unified_id = f"{chrom}:{cluster_start}-{cluster_end}:{strand}"
        for pos in cluster_positions:
            cluster_rows.append(
                {
                    "chromosome": chrom,
                    "effective_gene_id": effective_gene_id,
                    "strand": strand,
                    "site_position": pos,
                    "cluster_key": cluster_key,
                    "unified_ID": unified_id,
                    "cluster_start": cluster_start,
                    "cluster_end": cluster_end,
                    "cluster_width": cluster_width,
                    "n_unique_positions": len(cluster_positions),
                    **mode,
                }
            )

    return cluster_rows


def build_clusters_dp(
    df: pd.DataFrame,
    max_cluster_width: int,
    max_internal_gap: int,
    mode_support_window: int,
    sample_weight: float,
    count_weight: float,
    mode_sample_weight: float,
    mode_count_weight: float,
    dp_cluster_penalty: float,
    dp_compactness_weight: float,
    min_observation_count: int,
    strong_single_observation_count: int,
    strong_single_min_relative_abundance: float,
    absorb_weak_singletons: bool,
) -> pd.DataFrame:
    """Build a position-to-cluster map with DP segmentation."""
    rows: List[dict] = []
    grouped = df.groupby(GENE_KEY, sort=True, dropna=False)
    total_groups = grouped.ngroups
    log.info("DP clustering %d chromosome/gene/strand groups", total_groups)

    for i, (key, group) in enumerate(grouped, 1):
        rows.extend(
            _cluster_one_gene_dp(
                key,
                group,
                max_cluster_width,
                max_internal_gap,
                mode_support_window,
                sample_weight,
                count_weight,
                mode_sample_weight,
                mode_count_weight,
                dp_cluster_penalty,
                dp_compactness_weight,
                min_observation_count,
                strong_single_observation_count,
                strong_single_min_relative_abundance,
                absorb_weak_singletons,
            )
        )
        if i % 1000 == 0 or i == total_groups:
            log.info("  DP clustered %d / %d gene groups", i, total_groups)

    cluster_map = pd.DataFrame(rows, columns=CLUSTER_COLUMNS)
    log.info(
        "Built %d DP clusters from %d unique positions",
        cluster_map["cluster_key"].nunique(),
        len(cluster_map),
    )
    return cluster_map


def _renumber_retained_cluster_keys(
    detailed: pd.DataFrame,
    filter_info: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if detailed.empty:
        return detailed, filter_info.copy()

    cluster_order = (
        detailed[
            [
                "cluster_key",
                "chromosome",
                "effective_gene_id",
                "strand",
                "cluster_start",
                "cluster_end",
            ]
        ]
        .drop_duplicates()
        .sort_values(
            [
                "chromosome",
                "effective_gene_id",
                "strand",
                "cluster_start",
                "cluster_end",
                "cluster_key",
            ],
            kind="mergesort",
        )
    )

    mapping: Dict[str, str] = {}
    grouped = cluster_order.groupby(
        ["chromosome", "effective_gene_id", "strand"], sort=False, dropna=False
    )
    for (chrom, effective_gene_id, strand), group in grouped:
        for rank, row in enumerate(group.itertuples(index=False), 1):
            mapping[str(row.cluster_key)] = (
                f"{chrom}|{effective_gene_id}|{strand}|{rank}"
            )

    detailed = detailed.copy()
    detailed["cluster_key"] = detailed["cluster_key"].map(
        lambda key: mapping.get(str(key), key)
    )

    output_filter_info = filter_info[
        filter_info["cluster_key"].astype(str).isin(mapping)
    ].copy()
    output_filter_info["cluster_key"] = output_filter_info["cluster_key"].map(
        lambda key: mapping.get(str(key), key)
    )
    return detailed, output_filter_info


# ---------------------------------------------------------------------------
# Filtering and summaries
# ---------------------------------------------------------------------------


def assign_clusters(df: pd.DataFrame, cluster_map: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=DETAILED_COLUMNS)
    detailed = df.merge(cluster_map, on=JOIN_KEY, how="inner")
    detailed = detailed.rename(columns={"site_position": "original_site_position"})

    columns = [
        "transcript_id",
        "gene_id",
        "gene_name",
        "effective_gene_id",
        "chromosome",
        "strand",
        "cluster_key",
        "unified_ID",
        "original_site_position",
        "cluster_start",
        "cluster_end",
        "cluster_width",
        "mode_site_position",
        "mode_score",
        "mode_local_n_samples",
        "mode_local_total_site_count",
        "mode_exact_n_samples",
        "mode_exact_total_site_count",
        "n_unique_positions",
        "site_count",
        "transcript_biotype",
        "sample_attribute",
        "sample",
        SAMPLE_UID_COLUMN,
    ]
    detailed = detailed[columns]
    return detailed


def recalculate_cluster_relative_abundance(detailed: pd.DataFrame) -> pd.DataFrame:
    """Recalculate APA relative abundance after cluster filtering.

    The denominator is the sum of retained cluster counts for the same
    transcript and sample. Exact-position rows inside one cluster share the same
    cluster-level abundance.
    """
    if detailed.empty:
        detailed["cluster_relative_abundance"] = pd.Series(dtype=float)
        return detailed

    keys = ["transcript_id", SAMPLE_UID_COLUMN, "cluster_key"]
    cluster_counts = (
        detailed.groupby(keys, dropna=False)["site_count"]
        .sum()
        .rename("cluster_count_for_abundance")
        .reset_index()
    )
    denominators = (
        cluster_counts.groupby(
            ["transcript_id", SAMPLE_UID_COLUMN], dropna=False
        )[
            "cluster_count_for_abundance"
        ]
        .sum()
        .rename("transcript_sample_total_count")
        .reset_index()
    )
    abundance = cluster_counts.merge(
        denominators,
        on=["transcript_id", SAMPLE_UID_COLUMN],
        how="left",
    )
    abundance["cluster_relative_abundance"] = (
        abundance["cluster_count_for_abundance"]
        / abundance["transcript_sample_total_count"]
    )

    detailed = detailed.merge(
        abundance[
            keys
            + [
                "cluster_count_for_abundance",
                "transcript_sample_total_count",
                "cluster_relative_abundance",
            ]
        ],
        on=keys,
        how="left",
    )
    detailed["cluster_relative_abundance"] = (
        detailed["cluster_relative_abundance"].round(4)
    )
    detailed = detailed.drop(
        columns=["cluster_count_for_abundance", "transcript_sample_total_count"]
    )
    return detailed


def evaluate_cluster_filters(
    detailed: pd.DataFrame,
    min_single_transcript_count: int,
    min_multi_transcript_count: int,
    min_supported_transcripts: int,
    min_observation_count: int,
    min_supported_observations: int,
    strong_single_observation_count: int,
    strong_single_min_relative_abundance: float,
    min_cluster_samples: int,
    filter_strategy: str,
    no_filter: bool,
) -> pd.DataFrame:
    columns = ["cluster_key", "filter_passed", "filter_pass_rule", "filter_reason"]
    if detailed.empty:
        return pd.DataFrame(columns=columns)

    if (
        filter_strategy == "observations"
        and strong_single_observation_count > 0
        and "prefilter_cluster_relative_abundance" not in detailed.columns
    ):
        raise ValueError(
            "strong single-observation rescue requires pre-filter "
            "prefilter_cluster_relative_abundance"
        )

    cluster_keys = pd.Series(
        detailed["cluster_key"].drop_duplicates().to_numpy(), name="cluster_key"
    )
    if no_filter:
        return pd.DataFrame(
            {
                "cluster_key": cluster_keys,
                "filter_passed": True,
                "filter_pass_rule": "no_filter",
                "filter_reason": "",
            }
        )

    if filter_strategy == "observations":
        cluster_stats = (
            detailed.groupby("cluster_key", sort=False, dropna=False)[
                SAMPLE_UID_COLUMN
            ]
            .nunique()
            .rename("n_samples")
            .reset_index()
        )
        observations = (
            detailed.groupby(
                ["cluster_key", "transcript_id", SAMPLE_UID_COLUMN],
                sort=False,
                dropna=False,
            )
            .agg(
                observation_count=("site_count", "sum"),
                observation_relative_abundance=(
                    "prefilter_cluster_relative_abundance",
                    "max",
                ),
            )
            .reset_index()
        )
        supported_observations = (
            observations.loc[
                observations["observation_count"] >= min_observation_count
            ]
            .groupby("cluster_key", sort=False, dropna=False)
            .size()
            .rename("supported_observations")
            .reset_index()
        )
        max_observation = (
            observations.groupby("cluster_key", sort=False, dropna=False)[
                "observation_count"
            ]
            .max()
            .rename("max_observation")
            .reset_index()
        )

        stats = (
            pd.DataFrame({"cluster_key": cluster_keys})
            .merge(cluster_stats, on="cluster_key", how="left")
            .merge(supported_observations, on="cluster_key", how="left")
            .merge(max_observation, on="cluster_key", how="left")
        )
        stats["n_samples"] = stats["n_samples"].fillna(0).astype(int)
        stats["supported_observations"] = (
            stats["supported_observations"].fillna(0).astype(int)
        )
        stats["max_observation"] = stats["max_observation"].fillna(0.0)

        if strong_single_observation_count > 0:
            strong_abundance = (
                observations.loc[
                    observations["observation_count"]
                    >= strong_single_observation_count
                ]
                .groupby("cluster_key", sort=False, dropna=False)[
                    "observation_relative_abundance"
                ]
                .max()
                .rename("max_strong_observation_abundance")
                .reset_index()
            )
            stats = stats.merge(strong_abundance, on="cluster_key", how="left")
            stats["max_strong_observation_abundance"] = stats[
                "max_strong_observation_abundance"
            ].fillna(0.0)
            strong_single = (
                (stats["max_observation"] >= strong_single_observation_count)
                & (
                    stats["max_strong_observation_abundance"]
                    >= strong_single_min_relative_abundance
                )
            )
        else:
            strong_single = pd.Series(False, index=stats.index)

        enough_samples = stats["n_samples"] >= min_cluster_samples
        supported_pass = (
            stats["supported_observations"] >= min_supported_observations
        )
        passed = enough_samples & (supported_pass | strong_single)

        stats["filter_passed"] = passed
        stats["filter_pass_rule"] = ""
        stats.loc[
            enough_samples & supported_pass,
            "filter_pass_rule",
        ] = "supported_transcript_sample_observations"
        stats.loc[
            enough_samples & ~supported_pass & strong_single,
            "filter_pass_rule",
        ] = "strong_single_transcript_sample_observation_with_abundance"

        stats["filter_reason"] = ""
        stats.loc[
            ~enough_samples,
            "filter_reason",
        ] = f"n_samples_lt_{min_cluster_samples}"
        fail_reason = (
            f"fewer_than_{min_supported_observations}_transcript_sample"
            f"_observations_with_count_ge_{min_observation_count}"
        )
        if strong_single_observation_count > 0:
            fail_reason += (
                f"_and_max_observation_lt_{strong_single_observation_count}"
                f"_or_relative_abundance_lt_"
                f"{strong_single_min_relative_abundance:g}"
            )
        stats.loc[
            enough_samples & ~passed,
            "filter_reason",
        ] = fail_reason

        return stats[columns]

    records: List[dict] = []
    for cluster_key, group in detailed.groupby("cluster_key", sort=False):
        n_samples = group[SAMPLE_UID_COLUMN].nunique()
        transcript_counts = (
            group.groupby("transcript_id")["site_count"].sum().sort_values(ascending=False)
        )
        n_transcripts = len(transcript_counts)

        passed = False
        pass_rule = ""
        reason = ""

        if n_samples < min_cluster_samples:
            reason = f"n_samples_lt_{min_cluster_samples}"
        else:
            if n_transcripts == 1:
                count = float(transcript_counts.iloc[0])
                if count >= min_single_transcript_count:
                    passed = True
                    pass_rule = "single_transcript_count"
                else:
                    reason = f"single_transcript_count_lt_{min_single_transcript_count}"
            else:
                supported = transcript_counts[
                    transcript_counts >= min_multi_transcript_count
                ]
                if len(supported) >= min_supported_transcripts:
                    passed = True
                    pass_rule = "multi_transcript_count"
                else:
                    reason = (
                        f"fewer_than_{min_supported_transcripts}_transcripts"
                        f"_with_count_ge_{min_multi_transcript_count}"
                    )

        records.append(
            {
                "cluster_key": cluster_key,
                "filter_passed": passed,
                "filter_pass_rule": pass_rule,
                "filter_reason": reason,
            }
        )

    return pd.DataFrame(records)


def _drop_internal_filter_columns(detailed: pd.DataFrame) -> pd.DataFrame:
    return detailed.drop(
        columns=["prefilter_cluster_relative_abundance"],
        errors="ignore",
    )


def _capitalize_first(value):
    if pd.isna(value):
        return value
    text = str(value)
    if not text:
        return text
    return text[0].upper() + text[1:]


def _capitalize_tissue_sample_names(detailed: pd.DataFrame) -> pd.DataFrame:
    if detailed.empty or "sample" not in detailed or "sample_attribute" not in detailed:
        return detailed

    detailed = detailed.copy()
    tissue_mask = detailed["sample_attribute"].eq("tissue") & detailed[
        "sample"
    ].notna()
    detailed.loc[tissue_mask, "sample"] = detailed.loc[tissue_mask, "sample"].map(
        _capitalize_first
    )
    return detailed


def _format_number(value: float) -> str:
    if pd.isna(value):
        return ""
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:.4f}".rstrip("0").rstrip(".")


def _weighted_mean_and_sd(values: np.ndarray, weights: np.ndarray) -> Tuple[float, float]:
    if len(values) == 0 or weights.sum() <= 0:
        return (math.nan, math.nan)
    mean = float(np.average(values, weights=weights))
    variance = float(np.average((values - mean) ** 2, weights=weights))
    return mean, math.sqrt(variance)


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    if len(values) == 0 or weights.sum() <= 0:
        return math.nan
    order = np.argsort(values)
    sorted_values = values[order]
    sorted_weights = weights[order]
    cumulative = np.cumsum(sorted_weights)
    cutoff = quantile * sorted_weights.sum()
    return float(sorted_values[np.searchsorted(cumulative, cutoff, side="left")])


def _transcript_counts_string(group: pd.DataFrame) -> str:
    counts = group.groupby("transcript_id")["site_count"].sum()
    items = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    return ",".join(f"{tid}:{_format_number(count)}" for tid, count in items)


def summarize_clusters(
    detailed: pd.DataFrame,
    filter_info: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    filter_map = {}
    if filter_info is not None and not filter_info.empty:
        filter_map = filter_info.set_index("cluster_key").to_dict(orient="index")

    records: List[dict] = []
    for cluster_key, group in detailed.groupby("cluster_key", sort=False):
        first = group.iloc[0]
        transcript_ids = sorted(group["transcript_id"].dropna().unique())
        sample_list = sorted(group["sample"].dropna().unique())
        biotypes = sorted(group["transcript_biotype"].dropna().unique())
        positions = group["original_site_position"].to_numpy(dtype=float)
        weights = group["site_count"].to_numpy(dtype=float)
        weighted_mean, weighted_sd = _weighted_mean_and_sd(positions, weights)
        q25 = _weighted_quantile(positions, weights, 0.25)
        q75 = _weighted_quantile(positions, weights, 0.75)
        filter_row = filter_map.get(cluster_key, {})

        records.append(
            {
                "cluster_key": cluster_key,
                "unified_ID": first["unified_ID"],
                "gene_id": first["gene_id"],
                "gene_name": first["gene_name"],
                "effective_gene_id": first["effective_gene_id"],
                "chromosome": first["chromosome"],
                "strand": first["strand"],
                "cluster_start": int(first["cluster_start"]),
                "cluster_end": int(first["cluster_end"]),
                "cluster_width": int(first["cluster_width"]),
                "mode_site_position": int(first["mode_site_position"]),
                "mode_score": first["mode_score"],
                "mode_local_n_samples": int(first["mode_local_n_samples"]),
                "mode_local_total_site_count": _format_number(
                    first["mode_local_total_site_count"]
                ),
                "mode_exact_n_samples": int(first["mode_exact_n_samples"]),
                "mode_exact_total_site_count": _format_number(
                    first["mode_exact_total_site_count"]
                ),
                "n_unique_positions": int(first["n_unique_positions"]),
                "weighted_mean_position": round(weighted_mean, 2)
                if not math.isnan(weighted_mean)
                else "",
                "weighted_position_sd": round(weighted_sd, 2)
                if not math.isnan(weighted_sd)
                else "",
                "position_iqr": _format_number(q75 - q25)
                if not math.isnan(q25) and not math.isnan(q75)
                else "",
                "n_transcripts": len(transcript_ids),
                "transcript_ids": ",".join(transcript_ids),
                "transcript_cluster_counts": _transcript_counts_string(group),
                "total_site_count": _format_number(group["site_count"].sum()),
                "n_samples": group[SAMPLE_UID_COLUMN].nunique(),
                "n_tissues": group.loc[
                    group["sample_attribute"] == "tissue", SAMPLE_UID_COLUMN
                ].nunique(),
                "n_cell_cultures": group.loc[
                    group["sample_attribute"] == "cell_culture", SAMPLE_UID_COLUMN
                ].nunique(),
                "sample_list": ",".join(sample_list),
                "transcript_biotype": ",".join(biotypes),
                "median_cluster_relative_abundance": round(
                    float(group["cluster_relative_abundance"].median()), 4
                )
                if "cluster_relative_abundance" in group
                and group["cluster_relative_abundance"].notna().any()
                else "",
                "filter_pass_rule": filter_row.get("filter_pass_rule", ""),
                "filter_reason": filter_row.get("filter_reason", ""),
            }
        )

    summary = pd.DataFrame(records, columns=SUMMARY_COLUMNS)
    if summary.empty:
        return summary

    return summary.sort_values(
        ["chromosome", "cluster_start", "cluster_end", "strand", "effective_gene_id"],
        ignore_index=True,
    )


def sort_detailed(detailed: pd.DataFrame) -> pd.DataFrame:
    return detailed.sort_values(
        [
            "chromosome",
            "cluster_start",
            "cluster_end",
            "strand",
            "effective_gene_id",
            "transcript_id",
            "sample",
            "original_site_position",
        ],
        ignore_index=True,
    )


def warn_duplicate_coordinate_ids(summary: pd.DataFrame) -> None:
    if summary.empty:
        return
    gene_counts = summary.groupby("unified_ID")["effective_gene_id"].nunique()
    duplicated = gene_counts[gene_counts > 1]
    if not duplicated.empty:
        log.warning(
            "%d coordinate-only unified_ID values are shared by multiple effective gene ids.",
            len(duplicated),
        )


# ---------------------------------------------------------------------------
# CLI and main
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="unify_apa_sites",
        description=(
            "Merge single-sample PAS calls within ONE species using DP segmentation. "
            "All inputs must use the same genome build and transcript annotation."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Example:
    unify_apa_sites /path/to/Species_name --output-dir /path/to/unified
    python unify_apa_sites.py /path/to/Species_name --keep-scaffolds
    unify_apa_sites /path/to/Species_name --input-pattern '*/trek_rerun/*.apa_sites.txt'
        """,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument(
        "base_dir",
        type=Path,
        help="Species-level directory containing *apa_sites.txt files",
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=None,
        help="Directory for output TSV files (default: base_dir)",
    )
    parser.add_argument(
        "-p", "--output-prefix", default=None,
        help="Output prefix (default: {species}_unified_apa)",
    )
    parser.add_argument(
        "--input-pattern", default="*apa_sites.txt",
        help="Recursive input glob relative to the species directory (default: *apa_sites.txt)",
    )
    parser.add_argument(
        "--max-cluster-width",
        type=int,
        default=40,
        help="Maximum inclusive span of one APA cluster in nt (default: 40)",
    )
    parser.add_argument(
        "--max-internal-gap",
        type=int,
        default=24,
        help="Maximum allowed gap between adjacent observed positions inside one cluster (default: 24)",
    )
    parser.add_argument(
        "--mode-support-window",
        type=int,
        default=15,
        help="Window around candidate mode for local support scoring (default: 15)",
    )
    parser.add_argument(
        "--sample-weight",
        type=float,
        default=0.7,
        help="Mode score weight for local unique sample support (default: 0.7)",
    )
    parser.add_argument(
        "--count-weight",
        type=float,
        default=0.3,
        help="Mode score weight for local read count support (default: 0.3)",
    )
    parser.add_argument(
        "--mode-sample-weight",
        type=float,
        default=0.8,
        help="Final mode score weight for exact unique sample support (default: 0.8)",
    )
    parser.add_argument(
        "--mode-count-weight",
        type=float,
        default=0.2,
        help="Final mode score weight for exact read count support (default: 0.2)",
    )
    parser.add_argument(
        "--dp-cluster-penalty",
        type=float,
        default=225.0,
        help="Penalty for each DP cluster (default: 225)",
    )
    parser.add_argument(
        "--dp-compactness-weight",
        type=float,
        default=1.0,
        help="Weight for within-cluster coordinate dispersion (default: 1)",
    )
    parser.add_argument(
        "--no-absorb-weak-singletons",
        action="store_true",
        help="Disable post-DP absorption of weak singleton positions",
    )
    parser.add_argument(
        "--biotypes",
        default=",".join(DEFAULT_BIOTYPES),
        help="Comma-separated transcript biotypes to keep; include a trailing comma to keep blank biotypes (default: mRNA,lncRNA,lnc_RNA,protein_coding,blank)",
    )
    parser.add_argument(
        "--exclude-predicted-transcripts",
        action="store_true",
        help="Exclude RefSeq predicted transcripts with IDs starting XM_ or XR_",
    )
    parser.add_argument(
        "--keep-scaffolds",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Keep all sequence names, including chr1 and scaffolds (default: keep only NC_* accessions)",
    )
    parser.add_argument(
        "--keep-mitochondrial",
        action="store_true",
        help="Keep mitochondrial RefSeq accessions",
    )
    parser.add_argument(
        "--filter-strategy",
        choices=("observations", "transcript-totals"),
        default="observations",
        help=(
            "Confidence filter strategy. 'observations' keeps clusters with enough "
            "transcript-sample observations plus optional strong single-observation "
            "rescue; 'transcript-totals' applies per-transcript read-total "
            "rules instead (default: observations)"
        ),
    )
    parser.add_argument(
        "--min-observation-count",
        type=int,
        default=5,
        help=(
            "For --filter-strategy observations, required reads in one "
            "transcript/sample observation (default: 5)"
        ),
    )
    parser.add_argument(
        "--min-supported-observations",
        type=int,
        default=2,
        help=(
            "For --filter-strategy observations, minimum transcript/sample "
            "observations meeting --min-observation-count (default: 2)"
        ),
    )
    parser.add_argument(
        "--strong-single-observation-count",
        type=_parse_rescue_cutoff,
        default=50,
        metavar="auto|N",
        help=(
            "Strong single-observation rescue count. 'auto' uses min(50, 10 "
            "times the number of effective input files); 0 disables rescue "
            "(default: 50)"
        ),
    )
    parser.add_argument(
        "--strong-single-min-relative-abundance",
        type=float,
        default=0.25,
        help=(
            "For --filter-strategy observations, strong single-observation rescue "
            "also requires this pre-filter cluster relative abundance in the same "
            "transcript/sample. Ignored when --strong-single-observation-count is 0 "
            "(default: 0.25)"
        ),
    )
    parser.add_argument(
        "--min-single-transcript-count",
        type=int,
        default=10,
        help=(
            "For --filter-strategy transcript-totals, keep a single-transcript "
            "cluster only if that transcript has at least this many reads (default: 10)"
        ),
    )
    parser.add_argument(
        "--min-multi-transcript-count",
        type=int,
        default=5,
        help=(
            "For --filter-strategy transcript-totals, required reads per supported "
            "transcript in multi-transcript clusters (default: 5)"
        ),
    )
    parser.add_argument(
        "--min-supported-transcripts",
        type=int,
        default=2,
        help=(
            "For --filter-strategy transcript-totals, minimum transcripts meeting "
            "--min-multi-transcript-count (default: 2)"
        ),
    )
    parser.add_argument(
        "--min-cluster-samples",
        type=int,
        default=1,
        help="Minimum unique samples per cluster (default: 1)",
    )
    parser.add_argument(
        "--no-filter",
        action="store_true",
        help="Disable confidence filtering and keep all clusters",
    )
    parser.add_argument(
        "--write-unfiltered",
        action="store_true",
        help="Also write .unfiltered detailed and summary outputs",
    )
    parser.add_argument(
        "--skip-summary",
        action="store_true",
        help=(
            "Only write the filtered detailed *_unified_apa_sites.txt output; "
            "skip summary/audit files for faster large batch reruns"
        ),
    )
    return parser


def parse_args(argv: Optional[Sequence[str]] = None):
    return build_parser().parse_args(argv)


def _normalize_weights(sample_weight: float, count_weight: float) -> Tuple[float, float]:
    if not all(math.isfinite(w) and w >= 0 for w in (sample_weight, count_weight)):
        raise ValueError("Evidence weights must be finite and non-negative")
    total = sample_weight + count_weight
    if total <= 0:
        log.error("sample-weight + count-weight must be positive.")
        sys.exit(1)
    if abs(total - 1.0) > 1e-9:
        log.warning(
            "Normalizing mode score weights: sample=%s count=%s",
            sample_weight,
            count_weight,
        )
    return sample_weight / total, count_weight / total


def _validate_args(args: argparse.Namespace) -> None:
    """Validate command-line arguments before any file is read or written."""
    for name in ("dp_cluster_penalty", "dp_compactness_weight", "strong_single_min_relative_abundance"):
        if not math.isfinite(getattr(args, name)):
            raise ValueError(f"--{name.replace('_', '-')} must be finite")
    for name in ("min_single_transcript_count", "min_multi_transcript_count", "min_supported_transcripts"):
        if getattr(args, name) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be >= 1")
    if not args.input_pattern or Path(args.input_pattern).is_absolute() or ".." in Path(args.input_pattern).parts:
        raise ValueError("--input-pattern must be a nonempty glob within the species directory")
    if args.output_prefix is not None and (
        not args.output_prefix.strip() or args.output_prefix in (".", "..")
        or "/" in args.output_prefix or "\\" in args.output_prefix
    ):
        raise ValueError("--output-prefix must be a file-name prefix, not a path")
    if args.max_cluster_width < 1:
        log.error("--max-cluster-width must be >= 1")
        sys.exit(1)
    if args.max_internal_gap < 1:
        log.error("--max-internal-gap must be >= 1")
        sys.exit(1)
    if args.max_internal_gap > args.max_cluster_width - 1 and args.max_cluster_width > 1:
        log.warning(
            "--max-internal-gap (%d) is larger than useful for --max-cluster-width (%d)",
            args.max_internal_gap,
            args.max_cluster_width,
        )
    if args.mode_support_window < 0:
        log.error("--mode-support-window must be >= 0")
        sys.exit(1)
    if args.dp_cluster_penalty < 0:
        log.error("--dp-cluster-penalty must be >= 0")
        sys.exit(1)
    if args.dp_compactness_weight < 0:
        log.error("--dp-compactness-weight must be >= 0")
        sys.exit(1)
    if args.min_observation_count < 1:
        log.error("--min-observation-count must be >= 1")
        sys.exit(1)
    if args.min_supported_observations < 1:
        log.error("--min-supported-observations must be >= 1")
        sys.exit(1)
    if not 0 <= args.strong_single_min_relative_abundance <= 1:
        log.error("--strong-single-min-relative-abundance must be between 0 and 1")
        sys.exit(1)
    if args.min_cluster_samples < 1:
        log.error("--min-cluster-samples must be >= 1")
        sys.exit(1)


def _published_file_mode() -> int:
    """Regular-file permissions under the current umask (tempfile defaults to 0600)."""
    current = os.umask(0)
    os.umask(current)
    return 0o666 & ~current


def _write_table(frame: pd.DataFrame, path: Path) -> None:
    """Publish one complete TSV without truncating an existing result on failure."""
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="", dir=path.parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False,
        ) as handle:
            temporary = Path(handle.name)
            frame.to_csv(handle, sep="\t", index=False)
        os.chmod(temporary, _published_file_mode())
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def run_unification(
    args: argparse.Namespace,
    output_prefix: Optional[str] = None,
    capitalize_tissue_samples: bool = False,
) -> dict:
    """Run canonical DP unification for either supported command-line entry point."""
    _validate_args(args)
    base_dir = args.base_dir.expanduser().resolve()
    if not base_dir.is_dir():
        log.error("Base directory does not exist: %s", base_dir)
        sys.exit(1)
    output_dir = (args.output_dir or base_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_weight, count_weight = _normalize_weights(
        args.sample_weight, args.count_weight
    )
    mode_sample_weight, mode_count_weight = _normalize_weights(
        args.mode_sample_weight, args.mode_count_weight
    )
    biotypes = _parse_biotypes(args.biotypes)
    species_name = base_dir.name
    output_prefix = output_prefix or args.output_prefix or f"{species_name}_unified_apa"

    file_infos = discover_apa_files(base_dir, args.input_pattern)
    input_paths = {info["path"].resolve() for info in file_infos}
    output_suffixes = ("_sites.txt", "_summary.txt", "_filtered_out.txt",
                       "_sites.unfiltered.txt", "_summary.unfiltered.txt")
    for suffix in output_suffixes:
        dest = (output_dir / f"{output_prefix}{suffix}").resolve()
        if dest in input_paths:
            header = pd.read_csv(dest, sep="\t", nrows=0).columns
            if not {"unified_ID", "original_site_position"}.issubset(header):
                raise ValueError(f"Output would overwrite a single-sample input: {dest}")
    df = read_and_filter(
        file_infos=file_infos,
        biotypes=biotypes,
        keep_scaffolds=args.keep_scaffolds,
        keep_mitochondrial=args.keep_mitochondrial,
        exclude_predicted_transcripts=args.exclude_predicted_transcripts,
    )
    effective_sample_files = int(df[SAMPLE_UID_COLUMN].nunique())
    resolved_strong_single_count = resolve_strong_single_observation_count(
        args.strong_single_observation_count,
        effective_sample_files,
    ) if effective_sample_files else 0
    cutoff_mode = (
        "auto"
        if args.strong_single_observation_count == AUTO_RESCUE_CUTOFF
        else "override"
    )
    log.info(
        "Effective sample files: %d; strong single-observation cutoff: %d (%s)",
        effective_sample_files,
        resolved_strong_single_count,
        cutoff_mode,
    )

    cluster_map = build_clusters_dp(
        df=df,
        max_cluster_width=args.max_cluster_width,
        max_internal_gap=args.max_internal_gap,
        mode_support_window=args.mode_support_window,
        sample_weight=sample_weight,
        count_weight=count_weight,
        mode_sample_weight=mode_sample_weight,
        mode_count_weight=mode_count_weight,
        dp_cluster_penalty=args.dp_cluster_penalty,
        dp_compactness_weight=args.dp_compactness_weight,
        min_observation_count=args.min_observation_count,
        strong_single_observation_count=resolved_strong_single_count,
        strong_single_min_relative_abundance=(
            args.strong_single_min_relative_abundance
        ),
        absorb_weak_singletons=not args.no_absorb_weak_singletons,
    )
    detailed_all = assign_clusters(df, cluster_map)
    detailed_all = recalculate_cluster_relative_abundance(detailed_all)
    detailed_all = detailed_all.rename(
        columns={"cluster_relative_abundance": "prefilter_cluster_relative_abundance"}
    )

    filter_info = evaluate_cluster_filters(
        detailed=detailed_all,
        min_single_transcript_count=args.min_single_transcript_count,
        min_multi_transcript_count=args.min_multi_transcript_count,
        min_supported_transcripts=args.min_supported_transcripts,
        min_observation_count=args.min_observation_count,
        min_supported_observations=args.min_supported_observations,
        strong_single_observation_count=resolved_strong_single_count,
        strong_single_min_relative_abundance=args.strong_single_min_relative_abundance,
        min_cluster_samples=args.min_cluster_samples,
        filter_strategy=args.filter_strategy,
        no_filter=args.no_filter,
    )
    pass_keys = set(filter_info.loc[filter_info["filter_passed"], "cluster_key"])
    fail_keys = set(filter_info.loc[~filter_info["filter_passed"], "cluster_key"])

    detailed = detailed_all[detailed_all["cluster_key"].isin(pass_keys)].copy()
    detailed, output_filter_info = _renumber_retained_cluster_keys(
        detailed,
        filter_info,
    )
    detailed = recalculate_cluster_relative_abundance(detailed)
    if capitalize_tissue_samples:
        detailed = _capitalize_tissue_sample_names(detailed)
    dropped = detailed_all[detailed_all["cluster_key"].isin(fail_keys)].copy()

    out_detailed = output_dir / f"{output_prefix}_sites.txt"
    output_detailed = _drop_internal_filter_columns(detailed.copy())
    output_detailed = sort_detailed(output_detailed)
    _write_table(output_detailed, out_detailed)
    log.info(
        "Wrote filtered detailed output: %s (%d rows)",
        out_detailed,
        len(output_detailed),
    )

    retained_cluster_count = detailed["cluster_key"].nunique()
    if args.skip_summary:
        log.info("Skipped summary/audit outputs (--skip-summary)")
    else:
        summary = summarize_clusters(detailed, output_filter_info)
        filtered_out = summarize_clusters(dropped, filter_info)
        warn_duplicate_coordinate_ids(summary)

        out_summary = output_dir / f"{output_prefix}_summary.txt"
        out_filtered = output_dir / f"{output_prefix}_filtered_out.txt"

        _write_table(summary, out_summary)
        _write_table(filtered_out, out_filtered)

        retained_cluster_count = len(summary)
        log.info(
            "Wrote filtered summary output: %s (%d clusters)",
            out_summary,
            len(summary),
        )
        log.info(
            "Wrote filtered-out cluster audit: %s (%d clusters)",
            out_filtered,
            len(filtered_out),
        )

    if args.write_unfiltered:
        unfiltered_detailed = recalculate_cluster_relative_abundance(detailed_all.copy())
        if capitalize_tissue_samples:
            unfiltered_detailed = _capitalize_tissue_sample_names(unfiltered_detailed)
        unfiltered_summary = summarize_clusters(unfiltered_detailed, filter_info)
        unfiltered_detailed = _drop_internal_filter_columns(unfiltered_detailed)
        unfiltered_detailed = sort_detailed(unfiltered_detailed)
        out_unfiltered_detailed = output_dir / f"{output_prefix}_sites.unfiltered.txt"
        out_unfiltered_summary = output_dir / f"{output_prefix}_summary.unfiltered.txt"
        _write_table(unfiltered_detailed, out_unfiltered_detailed)
        _write_table(unfiltered_summary, out_unfiltered_summary)
        log.info(
            "Wrote unfiltered detailed output: %s (%d rows)",
            out_unfiltered_detailed,
            len(unfiltered_detailed),
        )
        log.info(
            "Wrote unfiltered summary output: %s (%d clusters)",
            out_unfiltered_summary,
            len(unfiltered_summary),
        )

    log.info(
        "Done. Retained %d / %d clusters.",
        retained_cluster_count,
        filter_info["cluster_key"].nunique(),
    )
    return {
        "detailed_path": out_detailed,
        "retained_clusters": retained_cluster_count,
        "total_clusters": int(filter_info["cluster_key"].nunique()),
        "effective_sample_files": effective_sample_files,
        "strong_single_observation_count": resolved_strong_single_count,
    }


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    try:
        run_unification(args, capitalize_tissue_samples=True)
    except (ValueError, OSError) as exc:
        log.error("%s", exc)
        raise SystemExit(1) from None


if __name__ == "__main__":
    main()
