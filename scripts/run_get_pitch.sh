#!/bin/bash

# Establecemos que el código de retorno de un pipeline sea el del último programa con código de retorno
# distinto de cero, o cero si todos devuelven cero.
set -o pipefail

# Put here the program (maybe with path)
GETF0="get_pitch"

# Array of threshold values
thresholds1=(0)
thresholds2=(0)

# Iterate over each combination of thresholds
for threshold1 in "${thresholds1[@]}"; do
    for threshold2 in "${thresholds2[@]}"; do
        echo "Running get_pitch with thresholds: $threshold1 and $threshold2"
        
        for fwav in pitch_db/train/*.wav; do
            ff0=${fwav/.wav/.f0}
            # echo "$GETF0 -t $threshold1 $fwav $ff0 ----"
            $GETF0 -t $threshold1 -r $threshold2 $fwav $ff0 > /dev/null || ( echo -e "\nError in $GETF0 -t $threshold1 -t $threshold2 $fwav $ff0" && exit 1 )
        done
        pitch_evaluate pitch_db/train/*.f0ref
    done
done

exit 0
0
