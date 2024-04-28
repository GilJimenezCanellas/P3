#!/bin/bash

# Establecemos que el código de retorno de un pipeline sea el del último programa con código de retorno
# distinto de cero, o cero si todos devuelven cero.
set -o pipefail

# Put here the program (maybe with path)
GETF0="get_pitch"

# Array of threshold values
th_rlags=(0.39 0.4 0.41)
th_r1s=(0.94 0.96 0.98)
th_zs=(32 33 34 35)
th_ms=(1 3)
th_cs=(0.0043 0.0045 0.0047 0.0049)

# Iterate over each combination of thresholds
for th_r1s in "${th_r1s[@]}"; do
    for th_rlags in "${th_rlags[@]}"; do
        for th_zs in "${th_zs[@]}"; do
            for th_ms in "${th_ms[@]}"; do
                for th_cs in "${th_cs[@]}"; do 
                    echo "Running get_pitch with thresholds: $th_r1s, $th_rlags, $th_zs, $th_ms and $th_cs"
                    
                    for fwav in pitch_db/train/*.wav; do
                        ff0=${fwav/.wav/.f0}
                        # echo "$GETF0 -t $threshold1 $fwav $ff0 ----"
                        $GETF0 -t $th_rlags -r $th_r1s -z $th_zs -m $th_ms -c $th_cs $fwav $ff0 > /dev/null || ( echo -e "\nError in $GETF0 -t $th_r1s -t $th_rlags -z $th_zs -m $th_ms -c $th_cs $fwav $ff0" && exit 1 )
                    done
                    pitch_evaluate pitch_db/train/*.f0ref
                done
            done
        done
    done
done

exit 0
0
