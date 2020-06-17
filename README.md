# spec2vec_gnps_data_analysis
Analysis and benchmarking of mass spectra similarity measures using gnps data set.

## Create environment
```
conda create --name spec2vec_analysis python=3.7
conda activate spec2vec_analysis
conda install --channel nlesc --channel bioconda --channel conda-forge matchms
conda install --channel nlesc --channel bioconda --channel conda-forge spec2vec
pip install jupyter
```

## Clone this repository and run notebooks
```
git clone https://github.com/iomega/spec2vec_gnps_data_analysis
cd spec2vec_gnps_data_analysis
jupyter notebook
```
