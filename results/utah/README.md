# New external evaluation: complete Utah humeri

Source: Henninger HB, University of Utah. *3D models of the human scapula and humerus with defined anatomic landmarks*, v1.11, 17 March 2026. https://doi.org/10.5281/zenodo.19077748 . Source licence: CC BY 4.0, https://creativecommons.org/licenses/by/4.0/ . Copyright remains with the source authors. The accompanying case inventory, selected metadata and numerical analysis are derived from that release; the cavity workflow was not authored or endorsed by the source authors. Original source meshes are not redistributed in the journal review archive.

## Selection and results

All five unpaired and all 120 paired-side complete humeri were inventoried before construction. There were no metadata exclusions or duplicate raw vertex sets. Five unpaired humeri plus 60 right sides form the primary one-bone-per-donor subset (65 source-defined donors). The 60 left sides are secondary observations on those same paired donors.

| Subset | Attempted | Source-evaluable | Sampled passes | Final accepted |
| --- | ---: | ---: | ---: | ---: |
| Primary, one bone per donor | 65 | 64 | 64 | 63 |
| Secondary left sides | 60 | 59 | 59 | 59 |
| Both subsets, descriptive only | 125 | 123 | 123 | 122 |

Primary acceptance per attempt was 63/65 (96.92%; exact two-sided 95% Clopper-Pearson CI 89.32-99.63%). Conditional on source evaluability it was 63/64 (98.44%; 91.60-99.96%). Do not compute a binomial interval treating the 125 bilateral bones as independent donors, or pool historical and revision stages.

The open sources were utah_paired_047_humerus_r and utah_paired_056_humerus_l. The cavity in utah_paired_008_humerus_r met sampled checks at 0.149906842 mm after contraction, but the continuous edge-edge minimum was 0.149529358 mm, below 0.150 mm, and the wrapper withheld it. None of these adverse cases was removed. Accepted primary outputs had minimum continuous clearance 0.360743 mm and required no contraction.

The local protocol and source selection were fixed before method outcomes; this was not a public preregistration. The existing v3 generator, preparation method, normalized external profile and audit were unchanged. Logs from a first preparation setup attempt without SciPy are retained separately. No anatomical or physical validation, public software release, author approval or journal acceptance is implied.

## Reproduction

Obtain database_v1.11.zip from the cited Zenodo version. Its published MD5 is fd87208ea7484b27b6ce0ebfba6711fb. Use an isolated Python environment with the versions in requirements.txt and Blender 5.0.1. The portable helper verifies the archive and source hashes, performs the preserved preparation and profile materialization, and invokes the checked launcher.

From the review repository root:

```sh
python code/reproduce_utah.py --archive /data/database_v1.11.zip --subset primary --dry-run
python code/reproduce_utah.py --archive /data/database_v1.11.zip --blender /path/to/blender --out /data/new_utah_run
```

The default run includes both subsets. Use --subset primary or --subset secondary_left to retain the fixed side selections, or repeat --case-id for a documented reproduction subset. New output must be empty and is kept separate from the archived evaluation. Exit status 2 means that at least one selected case was withheld or stopped; inspect each result rather than assuming a technical error or success. The original protocol, freeze records, source inventory, preparation statuses and per-case decisions are supplied alongside the numerical results.

Source hashes in freeze records refer to the original local files. Packaging removes local home paths from documentary JSON/logs and records original and packaged hashes in ARCHIVE_MANIFEST.json. The preserved numerical core, preparation helpers and normalized profile retain their source bytes. The portable reproduction helper does not reinterpret modified documentary paths as inputs.
