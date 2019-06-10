#!/bin/bash
set -e

# Test for reproducibility.
se_result="02_rsem/se/test_single.genes.results"
pe_result="02_rsem/pe/test_paired.genes.results"

snakemake -s startup.smk -p
snakemake -s ../../Snakefile --configfile config.yaml --config result_dir=result0 -p
snakemake -s ../../Snakefile --configfile config.yaml --config result_dir=result1 -p

se_a=$(md5sum result0/${se_result} | cut -d' ' -f1)
se_b=$(md5sum result1/${se_result} | cut -d' ' -f1)
pe_a=$(md5sum result0/${pe_result} | cut -d' ' -f1)
pe_b=$(md5sum result1/${pe_result} | cut -d' ' -f1)

if [ "$se_a" == "$se_b" ] && [ "$pe_a" == "$pe_b" ]; then
    curl https://gist.githubusercontent.com/dohlee/3ea154d52932b27d042566605a2cb9e2/raw/update_reproducibility.sh -H 'Cache-control: no-cache' | bash /dev/stdin -y
else
    curl https://gist.githubusercontent.com/dohlee/3ea154d52932b27d042566605a2cb9e2/raw/update_reproducibility.sh -H 'Cache-control: no-cache' | bash /dev/stdin -n
fi

echo "Test exited with $?."
