
#!/bin/bash

for i in {1..5}
do
   #Rscript 00_tidy_domains.R "$i" &
   #Rscript 03_proteomics_models.R "$i" &
   #Rscript 03_biochemistry_models.R "$i" &
   #Rscript 03_clinical_proteomics_models.R "$i" &
   #Rscript 03_clinical_biochemistry_models.R "$i" &
   Rscript 03_clinical_proteomics_biochemistry_models.R "$i" &
done
