#!/bin/bash

# Establecemos que el código de retorno de un pipeline sea el del último programa con código de retorno
# distinto de cero, o cero si todos devuelven cero.
set -o pipefail

# Put here the program (maybe with path)
GETF0="get_pitch"

# Array of threshold values
thresholds=(45 43 50 52 55)

# Iterate over each threshold value
for threshold in "${thresholds[@]}"; do
    echo "Running get_pitch with threshold: $threshold"
    for fwav in pitch_db/train/*.wav; do
        ff0=${fwav/.wav/.f0}
        # echo "$GETF0 -t $threshold $fwav $ff0 ----"
        $GETF0 -t $threshold $fwav $ff0 > /dev/null || ( echo -e "\nError in $GETF0 -t $threshold $fwav $ff0" && exit 1 )
    done

    pitch_evaluate pitch_db/train/*.f0ref
done

exit 0
