#!/bin/bash

tar -xzvf eschColiREL606-tRNAs.tar.gz eschColiREL606-tRNAs_name_map.txt
let lines=$(wc -l eschColiREL606-tRNAs_name_map.txt|cut -f 1 -d " ")-1
tail -$lines eschColiREL606-tRNAs_name_map.txt|cut -f 2|cut -f 2,3 -d '-'|sed 's/-/,/'|sed 's/fMet/Met/'|sed 's/Ile2/Ile/' > REL606tRNAPre.txt
