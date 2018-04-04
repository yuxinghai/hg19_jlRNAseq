#!/usr/bin/bash
awk '{if($3 != "gene") print $0}' gencode.v19.annotation.gtf| grep -v "^#" | gtfToGenePred /dev/stdin /dev/stdout | genePredToBed stdin gencode.v19.annotation.bed

