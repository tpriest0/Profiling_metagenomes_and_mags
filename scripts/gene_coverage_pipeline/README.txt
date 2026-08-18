A snakemake workflow that profiles a gene catalogue across metagenomes/metatranscriptomes.

If you have both metagenomes and metatranscriptomes available (which will allow for expression profiles to be generated later), then it is recommended that you run the workflow twice, once for each of the type of data that you have.

The best way to do this is to first run the workflow on metagenome samples (and adjust the config.yaml file so that the OUTPUT_DIR ends with something like 'MG_analysis') and then run on metatranscriptome samples (adjusting the OUTPUT_DIR to 'MT_analysis').
