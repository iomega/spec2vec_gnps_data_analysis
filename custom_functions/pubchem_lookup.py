import logging
import re
import pubchempy as pcp
import numpy as np
from matchms.utils import is_valid_inchikey


def pubchem_metadata_lookup(spectrum_in, name_search_depth=10, formula_search=False,
                            min_formula_length=6, formula_search_depth=25, verbose=1):
    """

    Parameters
    ----------
    spectrum_in
        Matchms type spectrum as input.
    name_search_depth: int
        How many of the most relevant name matches to explore deeper. Default = 10.

    """
    if spectrum_in is None:
        return None

    spectrum = spectrum_in.clone()
    if is_valid_inchikey(spectrum.get("inchikey")):
        return spectrum

    def _plausible_name(compound_name):
        return (isinstance(compound_name, str) and len(compound_name) > 4)

    compound_name = spectrum.get("compound_name")
    if not _plausible_name(compound_name):
        return spectrum

    # Start pubchem search
    inchi = spectrum.get("inchi")
    parent_mass = spectrum.get("parent_mass")
    if isinstance(parent_mass, np.ndarray):
        parent_mass = parent_mass[0]
    formula = spectrum.get("formula")

    # 1) Search for matching compound name
    results_pubchem = pubchem_name_search(compound_name, name_search_depth=name_search_depth,
                                          verbose=verbose)

    if len(results_pubchem) > 0:

        # 1a) Search for matching inchi
        if likely_has_inchi(inchi):
            inchi_pubchem, inchikey_pubchem, smiles_pubchem = find_pubchem_inchi_match(results_pubchem, inchi,
                                                                                       verbose=verbose)
        # 1b) Search for matching mass
        if not likely_has_inchi(inchi) or inchikey_pubchem is None:
            inchi_pubchem, inchikey_pubchem, smiles_pubchem = find_pubchem_mass_match(results_pubchem, parent_mass,
                                                                                      verbose=verbose)

        if inchikey_pubchem is not None and inchi_pubchem is not None:
            logging.info("Matching compound name: %s", compound_name)
            if verbose >= 1:
                print(f"Matching compound name: {compound_name}")
            spectrum.set("inchikey", inchikey_pubchem)
            spectrum.set("inchi", inchi_pubchem)
            spectrum.set("smiles", smiles_pubchem)
            return spectrum

        elif verbose >= 2:
            print(f"No matches found for compound name: {compound_name}")

    # 2) Search for matching formula
    if formula_search and formula and len(formula) >= min_formula_length:
        results_pubchem = pubchem_formula_search(formula, formula_search_depth=formula_search_depth,
                                                 verbose=verbose)

        if len(results_pubchem) > 0:

            # 2a) Search for matching inchi
            if likely_has_inchi(inchi):
                inchi_pubchem, inchikey_pubchem, smiles_pubchem = find_pubchem_inchi_match(results_pubchem, inchi)
            # 2b) Search for matching mass
            if inchikey_pubchem is None:
                inchi_pubchem, inchikey_pubchem, smiles_pubchem = find_pubchem_mass_match(results_pubchem, parent_mass)

            if inchikey_pubchem is not None and inchi_pubchem is not None:
                logging.info("Matching formula: %s", formula)
                if verbose >= 1:
                    print(f"Matching formula: {formula}")
                spectrum.set("inchikey", inchikey_pubchem)
                spectrum.set("inchi", inchi_pubchem)
                spectrum.set("smiles", smiles_pubchem)
                return spectrum

            elif verbose >= 2:
                print(f"No matches found for formula: {formula}")

    return spectrum


def likely_has_inchi(inchi):
    """Quick test to avoid excess in-depth testing"""
    if inchi is None:
        return False
    inchi = inchi.strip('"')
    regexp = r"(InChI=1|1)(S\/|\/)[0-9, A-Z, a-z,\.]{2,}\/(c|h)[0-9]"
    if not re.search(regexp, inchi):
        return False
    return True


def likely_inchi_match(inchi_1, inchi_2, min_agreement=3):
    """Try to match defective inchi to non-defective ones.

    Compares inchi parts seperately. Match is found if at least the first
    'min_agreement' parts are a good enough match.
    The main 'defects' this method accounts for are missing '-' in the inchi.
    In addition, differences between '-', '+', and '?'will be ignored.

    Parameters
    ----------
    inchi_1: str
        inchi of molecule.
    inchi_2: str
        inchi of molecule.
    min_agreement: int
        Minimum number of first parts that MUST be a match between both input
        inchi to finally consider it a match. Default is min_agreement=3.
    """
    if min_agreement < 2:
        print("Warning! 'min_agreement' < 2 has no discriminative power. Should be => 2.")
    if min_agreement == 2:
        print("Warning! 'min_agreement' == 2 has little discriminative power",
              "(only looking at structure formula. Better use > 2.")
    agreement = 0

    # Remove spaces and '"' to account for different notations.
    # Remove everything with little discriminative power.
    ignore_lst = ['"', ' ', '-', '+', '?']
    for ignore in ignore_lst:
        inchi_1 = inchi_1.replace(ignore, '')
        inchi_2 = inchi_2.replace(ignore, '')

    # Split inchi in parts.
    inchi_1_parts = inchi_1.split('/')
    inchi_2_parts = inchi_2.split('/')

    # Check if both inchi have sufficient parts (seperated by '/')
    if len(inchi_1_parts) >= min_agreement and len(
            inchi_2_parts) >= min_agreement:
        # Count how many parts agree well
        for i in range(min_agreement):
            agreement += (inchi_1_parts[i] == inchi_2_parts[i])

    if agreement == min_agreement:
        return True
    else:
        return False


def likely_inchikey_match(inchikey_1, inchikey_2, min_agreement=1):
    """Try to match inchikeys.

    Compares inchikey parts seperately. Match is found if at least the first
    'min_agreement' parts are a good enough match.

    Parameters
    ----------
    inchikey_1: str
        inchikey of molecule.
    inchikey_2: str
        inchikey of molecule.
    min_agreement: int
        Minimum number of first parts that MUST be a match between both input
        inchikey to finally consider it a match. Default is min_agreement=1.
    """
    if min_agreement not in [1, 2, 3]:
        print("Warning! 'min_agreement' should be 1, 2, or 3.")
    agreement = 0

    # Harmonize strings
    inchikey_1 = inchikey_1.upper().replace('"', '').replace(' ', '')
    inchikey_2 = inchikey_2.upper().replace('"', '').replace(' ', '')

    # Split inchikey in parts.
    inchikey_1_parts = inchikey_1.split('-')
    inchikey_2_parts = inchikey_2.split('-')

    # Check if both inchikey have sufficient parts (seperated by '/')
    if len(inchikey_1_parts) >= min_agreement and len(
            inchikey_2_parts) >= min_agreement:
        # Count how many parts mostly agree
        for i in range(min_agreement):
            agreement += (inchikey_1_parts[i] == inchikey_2_parts[i])

    return agreement == min_agreement


def pubchem_name_search(compound_name: str, name_search_depth=10, verbose=1):
    """Search pubmed for compound name"""
    results_pubchem = pcp.get_compounds(compound_name,
                                        'name',
                                        listkey_count=name_search_depth)
    if verbose >=2:
        print("Found at least", len(results_pubchem),
              "compounds of that name on pubchem.")
    return results_pubchem


def pubchem_formula_search(compound_formula: str, formula_search_depth=25, verbose=1):
    """Search pubmed for compound formula"""
    sids_pubchem = pcp.get_sids(compound_formula,
                                'formula',
                                listkey_count=formula_search_depth)

    results_pubchem = []
    for sid in sids_pubchem:
        result = pcp.Compound.from_cid(sid['CID'])
        results_pubchem.append(result)

    if verbose >=2:
        print(f"Found at least {len(results_pubchem)} compounds of with formula: {compound_formula}.")
    return results_pubchem


def find_pubchem_inchi_match(results_pubchem,
                             inchi,
                             min_inchi_match=3,
                             verbose=1):
    """Searches pubmed matches for inchi match.
    Then check if inchi can be matched to (defective) input inchi.


    Outputs found inchi and found inchikey (will be None if none is found).

    Parameters
    ----------
    results_pubchem: List[dict]
        List of name search results from Pubchem.
    inchi: str
        Inchi (correct, or defective...). Set to None to ignore.
    min_inchi_match: int
        Minimum number of first parts that MUST be a match between both input
        inchi to finally consider it a match. Default is min_inchi_match=3.
    """

    inchi_pubchem = None
    inchikey_pubchem = None
    smiles_pubchem = None

    # Loop through first 'name_search_depth' results found on pubchem. Stop once first match is found.
    for result in results_pubchem:
        inchi_pubchem = '"' + result.inchi + '"'
        inchikey_pubchem = result.inchikey
        smiles_pubchem = result.isomeric_smiles
        if smiles_pubchem is None:
            smiles_pubchem = result.canonical_smiles

        match_inchi = likely_inchi_match(inchi, inchi_pubchem,
                                         min_agreement=min_inchi_match)

        if match_inchi:
            logging.info("Matching inchi: %s", inchi)
            if verbose >= 1:
                print(f"Found matching compound for inchi: {inchi} (Pubchem: {inchi_pubchem}")
            break

    if not match_inchi:
        inchi_pubchem = None
        inchikey_pubchem = None
        smiles_pubchem = None

        if verbose >= 2:
            print("No matches found for inchi", inchi, "\n")

    return inchi_pubchem, inchikey_pubchem, smiles_pubchem


def find_pubchem_mass_match(results_pubchem,
                            parent_mass,
                            mass_tolerance=2.0,
                            verbose=1):
    """Searches pubmed matches for inchi match.
    Then check if inchi can be matched to (defective) input inchi.


    Outputs found inchi and found inchikey (will be None if none is found).

    Parameters
    ----------
    results_pubchem: List[dict]
        List of name search results from Pubchem.
    parent_mass: float
        Spectrum"s guessed parent mass.
    mass_tolerance: float
        Acceptable mass difference between query compound and pubchem result.
    """
    inchi_pubchem = None
    inchikey_pubchem = None
    smiles_pubchem = None

    for result in results_pubchem:
        inchi_pubchem = '"' + result.inchi + '"'
        inchikey_pubchem = result.inchikey
        smiles_pubchem = result.isomeric_smiles
        if smiles_pubchem is None:
            smiles_pubchem = result.canonical_smiles

        pubchem_mass = results_pubchem[0].exact_mass
        match_mass = (np.abs(pubchem_mass - parent_mass) <= mass_tolerance)

        if match_mass:
            logging.info("Matching molecular weight %s vs parent mass of %s",
                         str(np.round(pubchem_mass,1)),
                         str(np.round(parent_mass,1)))
            if verbose >= 1:
                print(f"Matching molecular weight ({pubchem_mass:.1f} vs parent mass of {parent_mass:.1f})")
            break

    if not match_mass:
        inchi_pubchem = None
        inchikey_pubchem = None
        smiles_pubchem = None

        if verbose >= 2:
            print(f"No matches found for mass {parent_mass} Da")

    return inchi_pubchem, inchikey_pubchem, smiles_pubchem
