#!/bin/bash



# for filename in "Data/"*.mzid; do
#     echo $filename
#     python3 proteoformquant.py -i "$filename" -s "${filename%.mzid}.mgf"
# done


MAX_THREADS=5

myArray=(1 2 5 10)

for filename in "Data/"*.mzid; do
    for p in ${myArray[@]}; do
		
		file_out=${filename%.mzid}_param_${p//./_}
		file_out=${file_out##*/}
		echo $file_out
	    python3 proteoformquant.py -i "$filename" -s "${filename%.mzid}.mgf" -o "$file_out" -d Output -max_rank ${p} &
	    while [ $( jobs | wc -l ) -ge "$MAX_THREADS" ]; do
		sleep 0.1
	    done
	done
done
 
