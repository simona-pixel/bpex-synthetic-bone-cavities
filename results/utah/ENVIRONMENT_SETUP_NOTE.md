# Preparation setup retry

The first preparation attempt stopped on all 125 inputs before source evaluation because trimesh diagnostics required SciPy, which was absent from the isolated Python environment. Logs and status records are retained in setup_attempt_missing_scipy/. No cavity construction occurred. SciPy is installed before an identical preparation retry; source selection, preparation code, profile, generator and checking criteria remain fixed. These are environment setup failures, not mesh or method outcomes.
