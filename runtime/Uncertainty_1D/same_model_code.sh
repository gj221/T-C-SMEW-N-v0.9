#!/usr/bin/env bash
# The ensemble's lower and upper sets run in this checkout; its original is
# the Base or SSP campaign result in another checkout. The three compare
# only if they ran the same model, so every file that decides the numbers
# must be byte for byte the same in both: the model code and its inputs,
# the driver, and the calibrated set. Tests, the forcing-preprocessing
# pipeline (Inputs/Preprocess, whose outputs the driver reads only as the
# finished forcing files), plotting, campaign wrappers and documentation
# may differ.
#
# Usage, from this checkout's root:
#   bash runtime/Uncertainty_1D/same_model_code.sh ORIGINAL_ROOT

set -euo pipefail

original="${1:?Usage: same_model_code.sh ORIGINAL_ROOT}"
here="$(pwd -P)"
if [[ ! -f "$here/v2/run/Run_One_Subfield.m" ]]; then
    echo "Run this from the T&C-SMEW-N v0.9 checkout root." >&2
    exit 2
fi
if [[ ! -f "$original/v2/run/Run_One_Subfield.m" ]]; then
    echo "No T&C-SMEW-N v0.9 checkout at $original." >&2
    exit 2
fi

status=0
diff -rq --exclude=Tests --exclude=Preprocess "$here/v2/T&C_Code" "$original/v2/T&C_Code" || status=1
diff -rq "$here/v2/run" "$original/v2/run" || status=1
for file in runtime/Base_1D/Campaign_Parameters.m runtime/Base_1D/Calcite_Case.m; do
    if ! cmp -s "$here/$file" "$original/$file"; then
        echo "Files $here/$file and $original/$file differ"
        status=1
    fi
done
if (( status != 0 )); then
    echo "The model differs between $here and $original." >&2
    exit 1
fi
echo "Model code and inputs: identical to $original"
