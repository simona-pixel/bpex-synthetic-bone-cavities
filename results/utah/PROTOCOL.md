# Additional external evaluation: Utah full humeri

This protocol is fixed locally before preparation or cavity construction on the new data. It is not a public preregistration. Dataset discovery and the repository README were inspected before this document was written; no cavity outcomes were available.

## Source and selection

- Henninger HB, University of Utah. *3D models of the human scapula and humerus with defined anatomic landmarks*, v1.11, 17 March 2026. DOI: https://doi.org/10.5281/zenodo.19077748. Zenodo API identifies the licence as CC BY 4.0.
- Use the dedicated `single_humeri_STLs` and `paired_humeri_STLs` groups, described as 5 unpaired complete humeri and 60 pairs of complete humeri. Do not select from the shoulder, scapular or kinematic groups; these contain partial bones and may overlap.
- Include every STL in these two full-humerus groups. Metadata inspection may identify an explicitly cropped/incomplete model or an exact duplicate; record such a model in the acquisition inventory and explain its exclusion before running the method. Do not exclude on age, sex, pathology, geometry quality or eventual outcome.
- Treat the two sides in a paired group as one donor. The primary analysis uses the right side of each pair plus every unpaired humerus, chosen without looking at method outcomes. The left sides are a secondary paired-side evaluation, not additional independent donors. Expected counts are 65 primary humeri and 60 secondary humeri, subject to metadata and duplicate verification.
- Inspect the source metadata and byte/geometry duplicates before construction, and preserve all case IDs and provenance. These groups were not used in the historical Figshare, VSD or ABH074 evaluations.

## Frozen method

- Preserve the historical `pipeline_v3.py`, continuous triangle-surface audit and normalized `EXTERNAL_PARAMETER_PROFILE_v1.json` without changes.
- Use the existing `prepare_mesh_v1.py` and `make_case_config_v1.py`. Only the historically allowed component selection, coincident-vertex merge, normal orientation and rigid PCA alignment are permitted. No hole filling, smoothing, remeshing, manual editing or case-specific parameter adjustment.
- Use the portable checked runner. Both the historical sampled pipeline and the continuous triangle-surface clearance audit must pass. The minimum clearance threshold is 0.150 mm, with the existing audit comparison tolerance of 1e-9 mm. Retain the historical guard's separate 1e-4 mm tolerance, unchanged.
- Run one evaluation per bone; reproducibility repetitions were tested separately and are not new observations.
- Resource setting: at most two concurrent Blender processes with two threads per process. A resource/technical failure may be retried with the same source/configuration; preserve the first failure and record the reason. No outcome-based tuning is permitted.

## Outcomes and denominators

- Report attempted inputs, source-evaluable inputs, pipeline rejections, continuous-audit rejections, technical failures, final accepted outputs, and the continuous minimum clearance for accepted outputs.
- The primary independent-donor estimate uses the one-bone-per-donor subset. Report exact two-sided Clopper-Pearson 95% binomial intervals for acceptance per attempted and per evaluable input. Do not use a 125-bone interval as if the two sides were independent subjects.
- Report the additional left-side results separately and, when applicable, donor-level acceptance of both sides.
- Preserve every failure and reason. A poor result is evidence about scope, not a reason to change the method or select a more successful dataset.
- Geometry acceptance supports digital portability on these surface meshes. It does not establish anatomical reconstruction of the medullary canal, mechanical/haptic performance, physical printability, clinical utility or production-grade reliability.

The JSON freeze record stores UTC time and SHA-256 hashes of this protocol, code, profile and source metadata before any preparation/run.
