# Pandora analysis pipeline

Snakemake pipeline to run pandora with and without denovo enabled.

The version used in the pandora paper has tag `pandora_paper_tag1`.

# Running

## Requirements

### Dependencies
* python 3.6+;
* singularity 2.4.1+;

### Setting up virtualenv
`./setup.sh`

## Running on the sample example:

1. Download sample data (TODO: add link);
2. `unzip sample_data.zip`
3. `source venv/bin/activate`
4. `bash scripts/run_pipeline_local.sh -j8`

## Running on the paper data:

## WARNING

`make_prg` memory limit on LSF clusters is recognised by killing the job with `KeyboardInterrupt`,
so try to avoid killing jobs when running with `--snakefile Snakefile_get_denovo_updated_prg` as much as possible.

1. `git checkout pandora_paper_tag1`
2. `source venv/bin/activate`

### If you want to run local:

3. `bash scripts/run_pipeline_local.sh -j <NB_OF_THREADS> --configfile config.pandora_paper_tag1.yaml`

### If you want to run on an LSF cluster:

3. `bash scripts/submit_lsf.sh --configfile config.pandora_paper_tag1.yaml`

# Troubleshooting

If you get this error:
```
Building DAG of jobs...
Pulling singularity image docker://leandroishilima/subsampler:pandora_paper_tag1.
WorkflowError:
Failed to pull singularity image from docker://leandroishilima/subsampler:pandora_paper_tag1:
WARNING: pull for Docker Hub is not guaranteed to produce the
WARNING: same image on repeated pull. Use Singularity Registry
WARNING: (shub://) to pull exactly equivalent images.
ERROR: Image file exists, not overwriting.

  File "/hps/nobackup2/iqbal/leandro/pandora_paper_tag1/subsampler/venv/lib/python3.7/site-packages/snakemake/deployment/singularity.py", line 88, in pull
```

pass to the running script the default location where singularity images are store.
For example, in the EBI cluster, it would be `--singularity-prefix /nfs/leia/singularity/leandro/`.

