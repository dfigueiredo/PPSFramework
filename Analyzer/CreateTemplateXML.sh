#!/bin/bash

search_dir="/eos/cms/store/group/phys_exotica/PPS-Exo"
today="$( date +"%Y%m%d" )" 

printf "\n[XML Condor]\n\tCreating a file to be used with the Condor Tool!\n"
printf "<?xml version=\"1.0\"?>\n<samples>\n" > condor_template_$today.xml

for entry in "$search_dir"/*_2021-02-*/*/*/*/*
do
  printf "\t<dataset>\n\t\t<name>"$entry"/</name>\n\t\t<era></era>\n\t\t<mode></mode>\n\t\t<xangle></xangle>\n\t\t<datatype></datatype>\n\t\t<output></output>\n\t\t<parameters>--none</parameters>\n\t\t<enable></enable>\n\t</dataset>\n\n" >> condor_template_$today.xml
done

printf "</samples>" >> condor_template_$today.xml

printf "\tTemplate file is condor_template_"$today".xml\n\n"

