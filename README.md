# SomaticMutationPipeline
A snakemake pipeline for calling somatic mutations from single-cell transcriptomes

## Installation and requirements
### Creating a new environment
You need to create the main environment where we will be running the snakemake pipeline
```
# Clone this repository
git clone https://github.com/MorganResearchLab/SomaticMutationPipeline.git
cd SomaticMutationPipeline

# Create the environment
conda env create -f environment.yml

# Activate the environment
conda activate somatic-mutation-pipeline
```

On top of this environment, you need to also clone the SComatic repository, as you will need to use the scripts and the data for your analysis. If you have access to the Morgan Lab shareds folder (`/uoa/scratch/shared/Morgan_Lab`) you can skip this process.

```
git clone https://github.com/cortes-ciriano-lab/SComatic.git
```

### Setting up the snakemake configuration
In order to run the workflow, you have to update the `config/config.yaml` according to your dataset. Most of the time you will just need to update these three parameters

#### BAM directory
```
paths:
    ...
    bam_dir: 
```

#### Cell type annotation file

#### Sample batch information


## Usage

## FAQs - Frequently asked questions

## Contacts
