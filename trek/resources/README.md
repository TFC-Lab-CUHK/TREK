# PAS motif reference tables

These tables are bundled so `annotate_apa` works without any external files.

- `species_motif_groups.tsv`: species-to-model assignments for 99 species.
  Columns: `species_key` (underscored binomial, e.g. `Homo_sapiens`),
  `species_name`, `major_lineage`, `motif_model_group`, `model_status`,
  `pa_confidence_policy`. A model is applied only when `model_status` is
  `reference_supported` and `motif_model_group` is one of `MAMMALIAN`,
  `LAND_PLANT`, `BUDDING_YEAST` or `FISSION_YEAST`; every other row resolves to
  `NOT_ASSESSED`. The table is a curated snapshot, not automatic taxonomic
  inference.
- `spombe_auxiliary_motifs.tsv`: the 53 auxiliary motifs used by the
  fission-yeast model, grouped into `UAG-containing`, `GUA-containing` and
  `GUA-UAG-containing` families with their search windows (Stroup and Ji,
  *Genome Res* 2024).

To annotate a species that is not listed, either pass a custom
`--species-model-table` with the same columns or choose `--motif-model`
explicitly. Use `NOT_ASSESSED` when no lineage model is appropriate; the
annotation then records sequence context and recurrence only.

Matching rules, motif lists and search windows are implemented in
`pas_motif_models.py`.
