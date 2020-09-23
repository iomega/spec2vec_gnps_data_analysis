# spec2vec_gnps_data_analysis
Analysis and benchmarking of mass spectra similarity measures using gnps data set.

If you use **spec2vec** for your research, please cite the following references:

F Huber, L Ridder, S Rogers, JJJ van der Hooft, "Spec2Vec: Improved mass spectral similarity scoring through learning of structural relationships", bioRxiv, https://doi.org/10.1101/2020.08.11.245928 

(and if you use **matchms** as well:
F. Huber, S. Verhoeven, C. Meijer, H. Spreeuw, E. M. Villanueva Castilla, C. Geng, J.J.J. van der Hooft, S. Rogers, A. Belloum, F. Diblen, J.H. Spaaks, (2020). matchms - processing and similarity evaluation of mass spectrometry data. Journal of Open Source Software, 5(52), 2411, https://doi.org/10.21105/joss.02411 )

Thanks!

## Create environment
Current spec2vec works with Python 3.7 or 3.8, it might also work with earlier versions but we haven't tested.
```
conda create --name spec2vec_analysis python=3.7  # or 3.8 if you prefer
conda activate spec2vec_analysis
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
