### Create conda environments

mkdir -p envs/conda

conda create --name mluci -channel conda-forge --channel bioconda --yes git=2.21.0 htop=2.2.0 pandas=0.24.2 r=3.5.1 snakemake=5.4.2
source activate mluci
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/mluci.yaml

source deactivate

# bamtools2_5_1
conda create --name bamtools2_5_1 --channel bioconda --yes bamtools=2.5.1
source activate bamtools2_5_1
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/bamtools2_5_1.yaml

source deactivate

# seqtk1_3
conda create --name seqtk1_3 --channel bioconda --yes seqtk=1.3
source activate seqtk1_3
conda env export --no-builds | sed '$d' | sed '$d' > envs/conda/seqtk1_3.yaml

source deactivate
