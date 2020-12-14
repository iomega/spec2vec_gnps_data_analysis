import os
import sys
import numpy as np
import pytest
from matchms import Spectrum
from spec2vec import Spec2Vec
from spec2vec import SpectrumDocument
#path_root = os.path.dirname(os.path.__file__)
path_root = os.path.dirname(os.getcwd())
sys.path.insert(0, os.path.join(path_root, "custom_functions"))
from custom_functions.library_search import library_matching


def test_library_matching():
    spectrum_1 = Spectrum(mz=np.array([100, 150, 200.]),
                          intensities=np.array([0.7, 0.2, 0.1]),
                          metadata={'precursor_mz': 500.5})
    spectrum_2 = Spectrum(mz=np.array([100, 140, 190.]),
                          intensities=np.array([0.4, 0.2, 0.1]),
                          metadata={'precursor_mz': 500.11})
    spectrum_3 = Spectrum(mz=np.array([100, 140, 190.]),
                          intensities=np.array([0.3, 0.5, 0.2]),
                          metadata={'precursor_mz': 501.1})
    spectrum_4 = Spectrum(mz=np.array([97.5, 137.5, 200.]),
                          intensities=np.array([0.8, 0.5, 0.4]),
                          metadata={'precursor_mz': 500.1})
    documents_library = [SpectrumDocument(s) for s in [spectrum_1, spectrum_2, spectrum_3]]
    documents_query = [SpectrumDocument(spectrum_4)]
    found_matches = library_matching(documents_query, documents_library,
                                     model=None,
                                     presearch_based_on=["precursor_mz"],
                                     include_scores=["cosine", "modcosine"],
                                     ignore_non_annotated=False,
                                     intensity_weighting_power=0.5,
                                     allowed_missing_percentage=5.0,
                                     cosine_tol=2.0,
                                     mass_tolerance=2.0,
                                     mass_tolerance_type="Dalton")

    scores_cosine = found_matches[0].values[:,0]
    expected_scores_cosine = np.array([0.05312127152597306, 0.0, 0.0])
    scores_modcos = found_matches[0].values[:,2]
    expected_scores_modcos = np.array([0.05312127152597306, 0.0, 0.7757282939050968])
    assert list(scores_cosine) == [pytest.approx(x, 1e-6) for x in expected_scores_cosine], \
        "Expected different scores."
    assert list(scores_modcos) == [pytest.approx(x, 1e-6) for x in expected_scores_modcos], \
        "Expected different mod. cosine scores."
    assert np.all(found_matches[0].values[:,3] == np.array([1, 0, 2])), \
        "Expected different number of matches"
    assert np.all(found_matches[0].values[:,4]), "Expected all mass matches to be True"
