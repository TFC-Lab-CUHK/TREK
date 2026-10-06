# trek (package)

This directory is the Python package behind the `trek`, `unify_apa_sites` and
`annotate_apa` commands. User documentation is the
[repository README](../README.md).

Module map:

| Module | Role |
| --- | --- |
| `run_trek.py` | `trek` command: annotation → alignment → read assignment → PAS calling → internal-priming filter |
| `gtf_processor.py`, `gff_processor.py` | Parse GTF / GFF3 into transcript models |
| `alignment_processor.py` | Run minimap2 and assign reads to transcripts |
| `apa_finder.py` | Gaussian-mixture clustering of read 3′ ends and peak filtering |
| `internal_priming_filter.py` | Remove A-rich (internal-priming) candidate sites |
| `unify_apa_sites.py` | `unify_apa_sites` command: merge per-sample PAS tables within a species |
| `annotate_apa.py`, `pas_motif_models.py`, `resources/` | `annotate_apa` command: sequence context, PAS motifs, evidence levels |
