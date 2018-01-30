# Drop-seq aligment pipeline

This is a Snakemake implementation of the Drop-seq alignment pipeline based on Drop-seq tools provided by the [McCarroll lab](http://mccarrolllab.com/dropseq/).

Paired-end raw data is expected to be stored in fastq format in directories named after samples. Currently the pipeline expects that files in these directories end with ".\*sample_1_sequence.txt.gz", respectively ".\*sample_2_sequence.txt.gz" for the first and second read. E.g.: "exp1/dropseq_exp1_1_sequence.txt.gz" & "exp1/dropseq_exp1_2_sequence.txt.gz". If this does not suit your setup, you can easily change this by adapting the fastq_to_bam input function in the Snakefile.

For further information on Drop-seq tools requirements, please refer to the Drop-seq Alignment Cookbook found under: http://mccarrolllab.com/dropseq/.

A good way to use this pipeline is in a conda environment, as this ensures that the correct dependencies are available. To set up the conda environment (assuming that [conda](https://conda.io/docs/user-guide/install/index.html#regular-installation) is installed) and to run the pipeline, follow these steps:

```
# clone pipeline into working directory
git clone https://github.com/argschwind/dropseq_alignment_pipeline.git path/to/workdir
cd path/to/workdir

# download Drop-seq tools
curl -O http://mccarrolllab.com/download/1276/Drop-seq_tools-1.13.zip
unzip Drop-seq_tools-1.13.zip && rm Drop-seq_tools-1.13.zip

# install other dependencies into an isolated environment
conda env create -n my-workflow-name -f environment.yml

# activate environment
source activate my-workflow-name

# adapt pipeline to your needs
vim Snakefile

# edit config file to add samples, alignment reference, and change parameters
vim config.yml

# execute entire pipeline for all samples in config file (-n = dryrun)
snakemake -n

# or execute individual processes by specifying target files
snakemake -n exp1/desired.output.file

# deactivate conda environment after work is done
source deactivate
```
