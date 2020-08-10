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

## Download data
- Original data was obtained from GNPS: https://gnps-external.ucsd.edu/gnpslibrary/ALL_GNPS.json
- Cleaned and processed GNPS dataset for positive mode spectra (raw data accessed on 2020-05-11), can be found on zenodo: https://zenodo.org/record/3978072

## Download pre-trained models
Pretrained Word2Vec models to be used with Spec2Vec can be found on zenodo.
- Model trained on __UniqueInchikey__ subset (12,797 spectra): https://zenodo.org/record/3978054
- Model trained on __AllPositive__ set of all positive ionization mode spectra (after filtering): https://zenodo.org/record/3978070
