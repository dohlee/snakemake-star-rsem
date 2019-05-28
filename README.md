# snakemake-star-rsem

STAR-RSEM pipeline in snakemake.

## Creating conda environment
You can create a new conda environment from `environment.yaml` file, 
```shell
$ conda env create -f environment.yaml && conda activate star-rsem
```

or you can use `--use-conda` option when executing `snakemake`.
```shell
$ snakemake --use-conda -p -j 32
```

## Removing conda environment
If you created a new conda environment, remove the environment when you are done with the pipeline.
```shell
$ conda env remove star-rsem
```
