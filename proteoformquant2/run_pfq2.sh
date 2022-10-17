#!/bin/bash



# for filename in "Data/"*.mzid; do
#     echo $filename
#     python3 proteoformquant2.py -i "$filename" -s "${filename%.mzid}.mgf"
# done


MAX_THREADS=4
for filename in "Data/"*.mzid; do
    python3 proteoformquant2.py -i "$filename" -s "${filename%.mzid}.mgf" &
    while [ $( jobs | wc -l ) -ge "$MAX_THREADS" ]; do
        sleep 0.1
    done
done
 