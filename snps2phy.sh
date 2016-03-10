#!/bin/bash 

#Tyler K. Chafin 
#December 4 2015 
#Converts .snps file from pyRAD to phylip format 
#Email: tkchafin@uark.edu with issues

if [ $1 ]; then 
  file="$1"; 
else 
  printf "\nUsage: $0 <.snps>\n\n";  
  exit 1; 
fi; 

#Format to phylip
sed -r 's/([0-9]+[A-Z]+[0-9]+) *([A-Z_-]+)/\1\t\2/g' $file | sed 's/ //g' | sed 's/_//g' >> $file.phy;
#Replace header line
sed -i -r 's/##([0-9]+).+,.*,([0-9]+).*/\1\t\2/g' $file.phy;


