# snakemake-star-rsem

[![Build Status](https://travis-ci.org/dohlee/snakemake-star-rsem.svg?branch=master)](https://travis-ci.org/dohlee/snakemake-star-rsem)
![Reproduction Status](https://img.shields.io/endpoint.svg?url=https://gist.githubusercontent.com/dohlee/75791354350a7c1b52358581a45af9ab/raw/star-rsem.json)

STAR-RSEM pipeline in snakemake.

## Reproducible pipeline
This pipeline guarantees reproducible results, as long as one uses the same random seed for the pipeline.
To set another random seed to get slightly different (but still reasonable) result, please change relevant random seed parameters in `config.yaml`.

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
$ conda env remove -n star-rsem
```

## Configuring config.yaml
Every parameter for tools should be configured in `config.yaml`. You should not modify snakemake rules in `rules` directory, unless you have good reason to tweak it.

## Preparing manifest file
All the samples that you are going to process should be specified in `manifest.csv` file. The first row should be a header, having 'name' and 'library\_layout' as mandatory fields.

You can change the name of the manifest file, but you have to change the value of `manifest` field in `config.yaml`.

## Running the pipeline
After setting all the configurations and manifest, you can dry-run the pipeline with `-n` option and see if the pipeline works appropriately.
Please note that even the pipeline seems OK for now, some unexpected errors may get you in trouble at runtime.
In this case, you should inspect log files in `logs` directory, and troubleshoot the error.
```shell
$ snakemake -n
```
If it seems OK, provide the pipeline with appropriate number of cores with `-j` option, and it will find optimal execution scenario that maximizes the usage of cores. 
`-p` option makes the actual commands to be printed out when they are to be executed.
```shell
$ snakemake -j 32 -p 
```

## Estimated runtime
It may be useful to have some rough estimates about the running time of the pipeline in mind.

| Rule | Input FASTQ size | Wallclock time (hms) |
| --- |:---:| ---:|
| star\_2\_pass\_single | 820M | 0:11:53 |
| star\_2\_pass\_paired | 1.8G,1.7G | 0:17:46 |
| rsem\_calculate\_expression\_se | 820M | 0.25:19 |
| rsem\_calculate\_expression\_pe | 1.8G,1.7G | 0:45:06 |

