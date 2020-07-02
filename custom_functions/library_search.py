from typing import List
import numpy as np
import pandas as pd
from matchms.similarity import CosineGreedyNumba, ModifiedCosine
from spec2vec import SpectrumDocument
from spec2vec import Spec2Vec
from spec2vec import Spec2VecParallel
from future_matchms import ParentmassMatchParallel


def library_matching(documents_query: List[SpectrumDocument],
                     documents_library: List[SpectrumDocument],
                     model,
                     presearch_based_on=["parentmass", "spec2vec-top10"],
                     ignore_non_annotated: bool = True,
                     include_scores=["spec2vec", "cosine", "modcosine"],
                     intensity_weighting_power: float = 0.5,
                     allowed_missing_percentage: float = 0,
                     cosine_tol: float = 0.005,
                     mass_tolerance: float = 1.0):
    """Selecting potential spectra matches with spectra library.

    Suitable candidates will be selected by 1) top_n Spec2Vec similarity, and 2)
    same precursor mass (within given mz_ppm tolerance(s)).
    For later matching routines, additional scores (cosine, modified cosine)
    are added as well.

    Args:
    --------
    documents_query:
        List containing all spectrum documents that should be queried against the library.
    documents_library:
        List containing all library spectrum documents.
    model:
        Pretrained word2Vec model.
    top_n: int, optional
        Number of entries witht the top_n highest Spec2Vec scores to keep as
        found matches. Default = 10.
    ignore_non_annotated: bool, optional
        If True, only annotated spectra will be considered for matching.
        Default = True.
    cosine_tol: float, optional
        Set tolerance for the cosine and modified cosine score. Default = 0.005
    mass_tolerance
        Specify tolerance for a parentmass match.
    """

    # Initializations
    found_matches = []
    m_mass_matches = None
    m_spec2vec_similarities = None

    def get_metadata(documents):
        metadata = []
        for doc in documents:
            metadata.append(doc._obj.get("smiles"))
        return metadata

    library_spectra_metadata = get_metadata(documents_library)
    if ignore_non_annotated:
        # Get array of all ids for spectra with smiles
        library_ids = np.asarray([i for i, x in enumerate(library_spectra_metadata) if x])
    else:
        library_ids = np.arange(len(documents_library))

    msg = "Presearch must be done either by 'parentmass' and/or 'spec2vec-topX'"
    assert "parentmass" in presearch_based_on or np.any(["spec2vec" in x for x in presearch_based_on]), msg

    # 1. Search for top-n Spec2Vec matches ------------------------------------
    if np.any(["spec2vec" in x for x in presearch_based_on]):
        top_n = int([x.split("top")[1] for x in presearch_based_on if "spec2vec" in x][0])
        print("Pre-selection includes spec2vec top {}.".format(top_n))
        spec2vec = Spec2VecParallel(model=model, intensity_weighting_power=intensity_weighting_power,
                                    allowed_missing_percentage=allowed_missing_percentage)
        m_spec2vec_similarities = spec2vec([documents_library[i] for i in library_ids], documents_query)
    
        # Select top_n similarity values:
        selection_spec2vec = np.argpartition(m_spec2vec_similarities, -top_n, axis=0)[-top_n:, :]
    else:
        selection_spec2vec = np.empty((0, len(documents_query)), dtype="int")

    # 2. Search for parent mass based matches ---------------------------------
    if "parentmass" in presearch_based_on:
        mass_matching = ParentmassMatchParallel(mass_tolerance)
        m_mass_matches = mass_matching([documents_library[i]._obj for i in library_ids],
                                       [x._obj for x in documents_query])
        selection_massmatch = np.where(m_mass_matches[:, i] == 1)[0]
    else:
        selection_massmatch = np.empty((0, len(documents_query)), dtype="int")

    # 3. Combine found matches ------------------------------------------------
    for i in range(len(documents_query)):
        s2v_top_ids = selection_spec2vec[:, i]
        mass_match_ids = selection_massmatch[:, i]

        all_match_ids = np.unique(np.concatenate((s2v_top_ids, mass_match_ids)))
        
        if len(all_match_ids) > 0:
            if "modcosine"in include_scores:
                # Get cosine score for found matches
                cosine_similarity = CosineGreedyNumba(tolerance=cosine_tol)
                cosine_scores = []
                for match_id in library_ids[all_match_ids]:
                    cosine_scores.append(cosine_similarity(documents_library[match_id]._obj,
                                                           documents_query[i]._obj))
            else:
                cosine_scores = len(all_match_ids) * ["not calculated"]

            if "cosine"in include_scores:
                # Get modified cosine score for found matches
                mod_cosine_similarity = ModifiedCosine(tolerance=cosine_tol)
                mod_cosine_scores = []
                for match_id in library_ids[all_match_ids]:
                    mod_cosine_scores.append(mod_cosine_similarity(documents_library[match_id]._obj,
                                                                   documents_query[i]._obj))
            else:
                mod_cosine_scores = len(all_match_ids) * ["not calculated"]

            matches_df = pd.DataFrame({"cosine_score": [x[0] for x in cosine_scores],
                                      "cosine_matches": [x[1] for x in cosine_scores],
                                      "mod_cosine_score": [x[0] for x in mod_cosine_scores],
                                      "mod_cosine_matches": [x[1] for x in mod_cosine_scores]},
                                      index=library_ids[all_match_ids])

            if m_mass_matches is not None:
                matches_df["mass_match"] = m_mass_matches[all_match_ids, i]

            if m_spec2vec_similarities is not None:
                matches_df["s2v_score"] = m_spec2vec_similarities[all_match_ids, i]
            elif "spec2vec"in include_scores:
                spec2vec_similarity = Spec2Vec(model=model, intensity_weighting_power=intensity_weighting_power,
                                               allowed_missing_percentage=allowed_missing_percentage)
                spec2vec_scores = []
                for match_id in library_ids[all_match_ids]:
                    spec2vec_scores.append(spec2vec_similarity(documents_library[match_id]._obj,
                                                               documents_query[i]._obj))
                matches_df["s2v_score"] = spec2vec_scores
            found_matches.append(matches_df.fillna(0))
        else:
            found_matches.append([])

    return found_matches
