#!/bin/bash



# for filename in "Data/"*.mzid; do
#     echo $filename
#     python3 proteoformquant2.py -i "$filename" -s "${filename%.mzid}.mgf"
# done


MAX_THREADS=8

myArray=(0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.90 0.95)

for filename in "Data/"*.mzid; do
    for p in ${myArray[@]}; do
		
		file_out=${filename%.mzid}_param_${p//./_}
		file_out=${file_out##*/}
		echo $file_out
	    python3 proteoformquant2.py -i "$filename" -s "${filename%.mzid}.mgf" -o "$file_out" -min_ep_fit ${p} &
	    while [ $( jobs | wc -l ) -ge "$MAX_THREADS" ]; do
		sleep 0.1
	    done
	done
done
 
