set -e

# Test for reproducibility.
se_result="02_rsem/se/test_single.genes.results"
pe_result="02_rsem/pe/test_paired.genes.results"

(
cd tests/test_reproducibility
snakemake -s startup.smk -p
snakemake -s ../../Snakefile --configfile config.yaml --config result_dir=result0 -p
snakemake -s ../../Snakefile --configfile config.yaml --config result_dir=result1 -p

a=$(md5sum result0/${se_result} | cut -d' ' -f1)
b=$(md5sum result1/${se_result} | cut -d' ' -f1)
[[ $a = $b ]]

a=$(md5sum result0/${se_result} | cut -d' ' -f1)
b=$(md5sum result1/${se_result} | cut -d' ' -f1)
[[ $a = $b ]]
)

echo "Test exited with $?."
