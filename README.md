# Synthetic bone cavity construction and checking

Review package for BPEX-106985, Generating and checking synthetic cavities in bone surface meshes with Blender and Python, by Simon Arvati, François Luc, Laurent Delbreilh, Agnès Liard-Zmuda and Alexis Lubet.

Version 1.1.0-review. This research-software repository contains code and numerical results for the associated manuscript, which is under peer review. The independently authored workflow code is available under the MIT License. Source-data and derived-geometry terms are documented separately. This software release does not imply journal acceptance.

## What the workflow does

The workflow constructs a configurable synthetic cavity from an aligned outer bone envelope. The checked entry point requires topology, intersection, sampled containment and sampled clearance checks followed by a whole-triangle minimum-distance audit. Candidate exports are distinct from accepted outputs. A result failing either stage is withheld.

The cavity is a prescribed engineering geometry. The accompanying paper reports computational tests and a digital hollow-segment demonstration. Printing, anatomical endosteal recovery, physical instrumentation, mechanical properties, radiographic realism and training outcomes require separate studies.

## Requirements

- Blender 5.0.1 build a3db93c5b259 was used for the recorded analyses, with embedded Python 3.11.13 observed in the revision. Blender includes NumPy and the geometry APIs used by the core.
- External Python 3.10 or later for the checked launcher; NumPy for independent winding-number calculations. The Utah preparation environment is pinned in results/utah/requirements.txt (Python 3.12.14, NumPy, trimesh, SciPy and NetworkX).

## Input preparation

Use triangulated outer geometry in millimetres. Follow the archived preparation protocol: largest connected component, deterministic tie-breaking, merged coincident vertices, consistent orientation, and principal-component alignment with the longitudinal axis on z. No hole repair, smoothing or anatomy completion was used in the external study.

The included portable profile records normalized settings. The historical prepare_mesh_v1.py helper uses NumPy, trimesh, SciPy and NetworkX to prepare geometry and a report. make_case_config_v1.py materializes the profile using the prepared length. A case configuration contains a single key matching the supplied case identifier. It must not be applied to a differently scaled or differently aligned input without repeating preparation and materialization.

```sh
python code/historical/prepare_mesh_v1.py --case-id example --source /data/example.stl --output /data/example_prepared.stl --report /data/example_prepared.json
python code/historical/make_case_config_v1.py --case-id example --prepared-report /data/example_prepared.json --profile configs/EXTERNAL_PARAMETER_PROFILE_v1.json --output /data/example_config.json
python code/run_checked.py --blender /path/to/blender --source /data/example_prepared.stl --config /data/example_config.json --case-id example --out results/example
```

The launcher refuses a nonempty output directory. It writes command and decision records and retains execution logs. Only a successful final decision creates an accepted directory containing outer.stl and cavity.stl. A withheld result exits with status 2. Inspect its reason; do not reinterpret candidate files as accepted results.

To rerun the archived sensitivity configurations after acquiring and preparing the six public inputs, use `python code/reproduce_sensitivity.py --prepared-inputs /data/prepared --blender /path/to/blender --out /data/new_run`. The prepared folder must contain `CASE_ID/prepared.stl` or `CASE_ID.stl`. Use `--dry-run` to inspect commands, or `--case-id` and `--variant` to select a subset. New output must be kept distinct from the recorded study results.

## Reusing the study data

Obtain source geometries from their original repositories and observe their licences:

- Keast, Bonacci and Fox, Statistical Shape Model of the Tibia, Figshare version 4, doi:10.6084/m9.figshare.20454462.v4. The archived source record specifies CC BY 4.0.
- Bertsatos and colleagues, Human Long Bone Dataset for csg-toolkit, doi:10.5281/zenodo.6614647. Bilateral long bones are from one donor; retain this shared donor identifier and check the source licence.
- MCM-Fischer, VSDFullBodyBoneModels version 3.0, doi:10.5281/zenodo.8316967. The source README specifies CC BY-NC-SA 4.0 for mesh MAT/XLSX data and EUPL v1.2 for its own MATLAB/Python source. Those licences do not automatically apply to this independently authored code.

- Henninger HB, University of Utah, 3D models of the human scapula and humerus with defined anatomic landmarks, v1.11, doi:10.5281/zenodo.19077748. CC BY 4.0. Complete humerus groups contain five unpaired specimens and 60 bilateral pairs; preserve shared donor IDs. Source metadata and attribution accompany results/utah/README.md.

No original public mesh or local patient-derived development mesh is included. The example hollow segment is an output derived from the CC BY Figshare tibia 102480, with attribution supplied in examples/ATTRIBUTION.md.

## New fixed-method Utah cohort

The source inventory contains 125 complete humeri from 65 source-defined donors. The fixed primary subset uses five unpaired humeri plus the right side of each pair: 63/65 attempts accepted, 63/64 source-evaluable inputs. The 60 secondary left sides yielded 59/60 attempts and 59/59 evaluable inputs. Two open sources and one continuous-audit rejection remain in the attempted denominators. Do not count the bilateral bones as 125 independent donors. Selection and method were locally fixed before construction, without retuning or repair.

Acquire database_v1.11.zip from the cited Zenodo release, install the preparation requirements in an isolated environment, and run:

```sh
python code/reproduce_utah.py --archive /data/database_v1.11.zip --subset primary --dry-run
python code/reproduce_utah.py --archive /data/database_v1.11.zip --blender /path/to/blender --out /data/new_utah_run
```

The portable route checks archive and source hashes, prepares each source, materializes the preserved profile and invokes the two-stage checked runner. A final exit status 2 includes source/method rejections; inspect per-case results. The original local orchestration is archived separately. See results/utah/README.md for protocol, paired-side analysis, attribution and documentary hash interpretation.

## Analyses and interpretation

results/sensitivity_results.csv contains 90 configuration attempts on six previously observed meshes from five donors. It is a paired exploratory sensitivity analysis, not a new independent cohort. One extended-coverage construction was rejected and two further outputs failed continuous clearance; all adverse cases remain in the derived results.

The historical initial Figshare result, amended Figshare reanalysis, one-donor Zenodo development and separate VSD evaluation have different roles and must not be pooled. Repeated executions do not add biological subjects. The 0.150-mm threshold is an engineering separation rule, not a validated manufacturing minimum.

The historical sampled guard permits 0.0001 mm below the requested threshold, while the preserved continuous auditor uses a comparison allowance of 0.000000001 mm. The checked wrapper requires both decisions. Neither allowance is a certified numerical error bound.

The independent winding-number comparison evaluates selected inside–outside decisions using a separate solid-angle implementation in external NumPy. Its agreement does not establish superiority over pseudonormals, fast winding numbers or other robust geometry systems. Historical within-environment byte identity must not be generalized across exporters: the revision baseline replay had identical serialized triangle coordinates, but different headers and normal records.

The revision scripts retain their analysis-specific case and configuration mappings. Local preparation scripts requiring the original study layout are identified as archival reproduction tools; the checked entry point above is portable and accepts explicit paths. The archive includes individual configurations so that the six public sensitivity cases can be prepared and rerun without access to local development anatomy.

## Repository layout

- code/core: unchanged numerical generator and continuous-distance auditor
- code: checked launchers and portable reproduction helpers
- code/historical: preserved source preparation and configuration materialization
- code/revision: numerical analysis and qualification sources from the original study layout
- configs: normalized profile and sensitivity/Utah configurations
- results: numerical tables, source decisions and audit records
- sources: published dataset metadata and source assessment
- examples: attributed digital hollow-segment STL and design scenario

Article drafts, correspondence, original submission archives, local source anatomies and patient-derived development geometry are excluded. Numerical records retain their original cohort roles; bilateral sides are not independent donors. The code licence and third-party attribution are described in LICENSE and LICENSE_STATUS.md.

The checked launchers and reproduction helpers accept explicit paths. Sources under code/revision document analyses using the original study layout and require the obtained geometries and locally resolved documentary paths. Source meshes and geometry-dependent query containers must be regenerated from the cited datasets under their respective licences.
