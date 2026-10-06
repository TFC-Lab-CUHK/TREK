#!/usr/bin/env python3
"""Lineage-aware sequence models for PA-signal annotation."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Iterable, Optional


MAMMALIAN = "MAMMALIAN"
LAND_PLANT = "LAND_PLANT"
BUDDING_YEAST = "BUDDING_YEAST"
FISSION_YEAST = "FISSION_YEAST"
NOT_ASSESSED = "NOT_ASSESSED"

IMPLEMENTED_MODELS = {MAMMALIAN, LAND_PLANT, BUDDING_YEAST, FISSION_YEAST}

CANONICAL_HEXAMERS = ["AATAAA", "ATTAAA"]
VARIANT_HEXAMERS = [
    "AGTAAA",
    "TATAAA",
    "CATAAA",
    "GATAAA",
    "AATATA",
    "AATACA",
    "AATAGA",
    "ACTAAA",
    "AAGAAA",
    "AATGAA",
    "TTTAAA",
]
MAMMALIAN_UPSTREAM_AUX = ["TGTA", "TATA", "ATAT", "TTTT"]
MAMMALIAN_DOWNSTREAM_AUX = ["TGTG", "GTGT", "TTTT", "GGGG"]

PLANT_NUE = ["AATAAA", "ATTAAA", "AGTAAA", "TATAAA", "AATGAA"]
PLANT_FUE = ["TTTGTA", "TTGTAT", "TTGTAA", "TTGTA", "TTGTT", "TGTGTA"]

BUDDING_EE = ["TATATA", "TATGTA", "TACATA", "TAAATA", "TATTTA", "ATATAT"]
BUDDING_PE = ["AATAAA", "AAAAAAAA", "TTAAGAAC", "AAGAA", "AATAATGA"]

FISSION_A_RICH = ["AATAAA", "TTAATA", "TAATAT", "ATTAAT"]


def load_species_models(table_path: Path) -> Dict[str, str]:
    models: Dict[str, str] = {}
    with Path(table_path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"species_key", "motif_model_group", "model_status"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"Species motif table lacks required columns {sorted(required)}: "
                f"{table_path}"
            )
        for row in reader:
            species_key = row["species_key"].strip()
            if not species_key:
                raise ValueError(f"Species motif table contains an empty species key: {table_path}")
            if species_key in models:
                raise ValueError(
                    f"Species {species_key!r} occurs more than once in motif model table: "
                    f"{table_path}"
                )
            group = row["motif_model_group"].strip()
            status = row["model_status"].strip()
            if status == "reference_supported" and group in IMPLEMENTED_MODELS:
                models[species_key] = group
            else:
                models[species_key] = NOT_ASSESSED
    return models


def resolve_species_model(species_key: str, table_path: Path) -> str:
    models = load_species_models(table_path)
    if species_key in models:
        return models[species_key]
    raise ValueError(f"Species {species_key!r} is absent from motif model table: {table_path}")


def load_spombe_auxiliary_motifs(table_path: Path) -> Dict[str, list[str]]:
    groups = {
        "UAG-containing": [],
        "GUA-containing": [],
        "GUA-UAG-containing": [],
    }
    with Path(table_path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"motif_dna", "motif_family"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"S. pombe motif table lacks required columns {sorted(required)}: "
                f"{table_path}"
            )
        for row in reader:
            family = row["motif_family"].strip()
            motif = row["motif_dna"].strip().upper()
            if family not in groups:
                raise ValueError(f"Unknown S. pombe motif family {family!r} in {table_path}")
            if motif:
                groups[family].append(motif)
    if any(not motifs for motifs in groups.values()):
        missing = sorted(family for family, motifs in groups.items() if not motifs)
        raise ValueError(f"S. pombe motif table has empty families: {missing}")
    return groups


def _slice_relative(
    sequence: str,
    context_start: int,
    rel_start: int,
    rel_end: int,
) -> tuple[str, int]:
    """Return a half-open relative window and its actual start after clipping."""
    context_end = context_start + len(sequence)
    actual_start = max(context_start, rel_start)
    actual_end = min(context_end, rel_end)
    if actual_end <= actual_start:
        return "", actual_start
    left = actual_start - context_start
    right = actual_end - context_start
    return sequence[left:right], actual_start


def _first_motif(
    sequence: str,
    context_start: int,
    rel_start: int,
    rel_end: int,
    motifs: Iterable[str],
) -> Optional[tuple[str, int]]:
    window, actual_start = _slice_relative(
        sequence, context_start, rel_start, rel_end
    )
    for motif in motifs:
        offset = window.find(motif)
        if offset >= 0:
            return motif, actual_start + offset
    return None


def _u_rich_hexamer(
    sequence: str,
    context_start: int,
    rel_start: int,
    rel_end: int,
) -> Optional[tuple[str, int]]:
    window, actual_start = _slice_relative(
        sequence, context_start, rel_start, rel_end
    )
    candidates = []
    for offset in range(max(0, len(window) - 5)):
        motif = window[offset : offset + 6]
        t_count = motif.count("T")
        if t_count >= 5:
            position = actual_start + offset
            candidates.append((-t_count, abs(position), position, motif))
    if not candidates:
        return None
    _negative_t_count, _distance, position, motif = min(candidates)
    return motif, position


def _context_metrics(sequence: str, context_start: int) -> tuple[str, Optional[float]]:
    cleavage, _ = _slice_relative(sequence, context_start, -1, 1)
    u_window, _ = _slice_relative(sequence, context_start, -10, 16)
    u_fraction = round(u_window.count("T") / len(u_window), 4) if u_window else None
    return cleavage, u_fraction


def _result(
    motif: Optional[str],
    position: Optional[int],
    motif_type: str,
    search_level: str,
    motif_support: str,
    cleavage_context: str,
    u_rich_fraction: Optional[float],
) -> dict:
    return {
        "motif": motif,
        "position": position,
        "motif_type": motif_type,
        "search_level": search_level,
        "motif_support": motif_support,
        "cleavage_context": cleavage_context,
        "u_rich_fraction": u_rich_fraction,
    }


def annotate_context(
    sequence: str,
    context_start: int,
    model: str,
    spombe_auxiliary: Optional[Dict[str, list[str]]] = None,
) -> dict:
    sequence = sequence.upper()
    cleavage_context, u_fraction = _context_metrics(sequence, context_start)

    if model == NOT_ASSESSED:
        return _result(
            None, None, "not_assessed", "not_assessed", "not_assessed", "", None
        )

    if model == MAMMALIAN:
        hit = _first_motif(sequence, context_start, -40, 0, CANONICAL_HEXAMERS)
        if hit:
            return _result(*hit, "canonical", "mammalian_hexamer", "primary", cleavage_context, None)
        hit = _first_motif(sequence, context_start, -40, 0, VARIANT_HEXAMERS)
        if hit:
            return _result(*hit, "variant", "mammalian_hexamer", "primary", cleavage_context, None)
        hit = _first_motif(sequence, context_start, -40, 0, MAMMALIAN_UPSTREAM_AUX)
        if hit:
            return _result(*hit, "upstream", "mammalian_upstream", "auxiliary", cleavage_context, None)
        hit = _first_motif(sequence, context_start, 1, 41, MAMMALIAN_DOWNSTREAM_AUX)
        if hit:
            return _result(*hit, "downstream", "mammalian_downstream", "auxiliary", cleavage_context, None)
        return _result(None, None, "none", "none", "none", cleavage_context, None)

    if model == LAND_PLANT:
        hit = _first_motif(sequence, context_start, -30, -9, PLANT_NUE)
        if hit:
            return _result(*hit, "plant_nue", "plant_nue", "primary", cleavage_context, u_fraction)
        hit = _first_motif(sequence, context_start, -150, -29, PLANT_FUE)
        if hit:
            return _result(*hit, "plant_fue", "plant_fue", "auxiliary", cleavage_context, u_fraction)
        return _result(None, None, "none", "none", "none", cleavage_context, u_fraction)

    if model == BUDDING_YEAST:
        efficiency = _first_motif(sequence, context_start, -60, -34, BUDDING_EE)
        positioning = _first_motif(sequence, context_start, -30, -9, BUDDING_PE)
        if efficiency and positioning:
            return _result(*positioning, "yeast_ee_pe", "yeast_ee_pe", "primary", cleavage_context, None)
        if positioning:
            return _result(*positioning, "yeast_positioning", "yeast_pe", "auxiliary", cleavage_context, None)
        if efficiency:
            return _result(*efficiency, "yeast_efficiency", "yeast_ee", "auxiliary", cleavage_context, None)
        hit = _u_rich_hexamer(sequence, context_start, -10, 11)
        if hit:
            return _result(*hit, "yeast_u_rich", "yeast_u_rich", "auxiliary", cleavage_context, None)
        return _result(None, None, "none", "none", "none", cleavage_context, None)

    if model == FISSION_YEAST:
        hit = _first_motif(sequence, context_start, -30, -9, FISSION_A_RICH)
        if hit:
            return _result(*hit, "pombe_a_rich", "pombe_a_rich", "primary", cleavage_context, None)
        if spombe_auxiliary is None:
            raise ValueError("S. pombe auxiliary motif families are required")
        combined = spombe_auxiliary["GUA-UAG-containing"]
        hit = _first_motif(sequence, context_start, -80, -30, combined)
        if not hit:
            hit = _first_motif(sequence, context_start, 15, 61, combined)
        if hit:
            return _result(*hit, "pombe_gua_uag", "pombe_gua_uag", "auxiliary", cleavage_context, None)
        hit = _first_motif(
            sequence,
            context_start,
            -80,
            -30,
            spombe_auxiliary["UAG-containing"],
        )
        if hit:
            return _result(*hit, "pombe_uag", "pombe_uag", "auxiliary", cleavage_context, None)
        hit = _first_motif(
            sequence,
            context_start,
            15,
            61,
            spombe_auxiliary["GUA-containing"],
        )
        if hit:
            return _result(*hit, "pombe_gua", "pombe_gua", "auxiliary", cleavage_context, None)
        hit = _u_rich_hexamer(sequence, context_start, -20, 21)
        if hit:
            return _result(*hit, "pombe_u_rich", "pombe_u_rich", "auxiliary", cleavage_context, None)
        return _result(None, None, "none", "none", "none", cleavage_context, None)

    raise ValueError(f"Unsupported motif model: {model!r}")
