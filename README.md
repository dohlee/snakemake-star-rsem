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

## Configuring config.yaml
Every parameter for tools should be configured in `config.yaml`. You should not modify snakemake rules in `rules` directory, unless you have good reason to tweak it.

## Preparing manifest file
All the samples that you are going to process should be specified in `manifest.csv` file. The first row should be a header, having 'name' and 'library\_layout' as mandatory fields.

You can change the name of the manifest file, but you have to change the value of `manifest` field in `config.yaml`.


