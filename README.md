# Pandora analysis pipeline

Snakemake pipeline to run pandora with and without denovo enabled.

# Running

## Requirements

### Dependencies
You need `singularity` and `python 3.6+`.

### Setting up virtualenv
```
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Running on the paper data:

Be sure you are on tag `pandora_paper_tag1`

If you are on a *LSF* cluster, run:
`bash scripts/submit_lsf.sh`

Otherwise, you can run locally as:
`bash scripts/run_pipeline_local.sh`
