#!/usr/bin/bash
nohup snakemake -j 30 --snakefile hg19_Snakefile.py >snake.log &
