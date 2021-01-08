!#/bin/bash

#To be used after running Cutadapt. Moves files that are too short 'cutadapt2Short_fw' directory

for f in $(ls *.fastq.gz); do  var=$(wc -l ${f}); var2=$(echo ${var} | cut -f 1 -d " "); if [[((${var2} -lt 17))]]; then mv ${f} ./cutadapt2Short_fw/; fi; done
