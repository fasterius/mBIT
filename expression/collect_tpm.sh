#!/bin/bash

# Combines TPM from an arbitrary number of files (ENSTID in first column).
# Both Salmon and Kallisto outputs can be collected, but the input directories
# should be on this format: <sample>.[kallisto|salmon]/[abundance.tsv/quant.sf]

# Check for missing input
if [ -z ${1+x} ]; then
    echo "ERROR: missing input file; aborting."
    echo "Please provide input directory to use as the first argument!"
    exit 1
fi

# Find files
FILES=$(find $1 -name 'abundance.t*' -or -name 'quant.sf' | xargs)

# Get replicate names for column header
printf "%s" 'ENSTID'
for FILE in $FILES; do
    NAME=${FILE/.\//}
    NAME2=$(echo $NAME | cut -d "." -f 1)
	printf "\t%s" $NAME2
done
printf "\n"

FILE=$(echo $FILES | xargs -n 1 | head -1)
COL=$(($(head -1 $FILE | xargs -n 1 | nl | grep -i "TPM" | cut -f 1)))

gawk -v COL=$COL \
    'BEGIN { OFS="\t" } 
	{ vals[$1,ARGIND]=$COL; keys[$1] } 
	END {
    		for (key in keys) {
        		printf "%s%s", key, OFS
        		for (colNr=1; colNr<=ARGIND; colNr++) {
        			printf "%s%s", vals[key,colNr], (colNr<ARGIND?OFS:ORS)
			    }
    		} 
	}' $FILES
