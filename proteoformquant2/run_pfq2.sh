#!/bin/bash



for filename in "Data/"*.mzid; do
    echo $filename
    python3 proteoformquant2.py -i "$filename" -s "${filename%.mzid}.mgf" -d "/Output"
done


#MAX_THREADS=1
#for filename in "Data/"*.mzid; do

#    python3 proteoformquant2.py -i "$filename" -s "${filename%.mzid}.mgf" -d "/Output" &
    
#    while [ $( jobs | wc -l ) -ge "$MAX_THREADS" ]; do
#        sleep 0.1
#    done
#done
 
