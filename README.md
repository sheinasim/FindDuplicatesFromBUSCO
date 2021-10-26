# FindDuplicatesFromBUSCO 
Print list of contigs identified as "duplicates" using the output from BUSCO run in -geno mode. 

This script is designed to be run on the full\_table.tsv output of BUSCO. 

Dependencies:

* Python3 
* Biopython
* Python pandas
* Python numpy  
* Python re
* Python os
* Python Scipy
* Python json

> export PATH=$PATH:[PATH TO FindDuplicatesFromBUSCO]  

### Usage
  
> python findDuplicates.py \<path to Blob dir\> \<filename or path to full\_table.tsv\> \<assembly.fasta\>

If no arguments are provided, the script will return help message.

## Outputs

* \{fasta prefix\}\_duplicates.txt 

### Citation

If this script is useful to you, please cite the following in your publication:

```
@software{FindDuplicatesFromBUSCO,
  author = {Sim, Sheina B.},
  title = {FindDuplicatesFromBUSCO},
  url = {https://github.com/sheinasim/FindDuplicatesFromBUSCO}
}
```

Sheina B. Sim  
USDA-ARS  
US Pacific Basin Agricultural Research Service  
Hilo, Hawaii, 96720 USA  
sheina.sim@usda.gov  

This script is in the public domain in the United States per 17 U.S.C. § 105
