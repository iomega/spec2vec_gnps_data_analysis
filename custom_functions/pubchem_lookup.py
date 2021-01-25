import pubchempy as pcp


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


def pubchem_lookup(spectrum_in):
    
    spectrum = spectrum_in.clone()

    



def find_pubchem_match(compound_name,
                       inchi,
                       inchikey=None,
                       mode='and',
                       min_inchi_match=3,
                       min_inchikey_match=1,
                       name_search_depth=10,
                       formula_search=False,
                       formula_search_depth=25):
    """Searches pubmed for compounds based on name.
    Then check if inchi and/or inchikey can be matched to (defective) input inchi and/or inchikey.

    In case no matches are found: For formula_search = True, the search will continue based on the
    formula extracted from the inchi.

    Outputs found inchi and found inchikey (will be None if none is found).

    Parameters
    ----------
    compound_name: str
        Name of compound to search for on Pubchem.
    inchi: str
        Inchi (correct, or defective...). Set to None to ignore.
    inchikey: str
        Inchikey (correct, or defective...). Set to None to ignore.
    mode: str
        Determines the final matching criteria (can be se to 'and' or 'or').
        For 'and' and given inchi AND inchikey, a match has to be a match with inchi AND inchikey.
        For 'or' it will be sufficient to find a good enough match with either inchi OR inchikey.
    min_inchi_match: int
        Minimum number of first parts that MUST be a match between both input
        inchi to finally consider it a match. Default is min_inchi_match=3.
    min_inchikey_match: int
        Minimum parts of inchikey that must be equal to be considered a match.
        Can be 1, 2, or 3.
    name_search_depth: int
        How many of the most relevant name matches to explore deeper. Default = 10.
    formula_search: bool
        If True an additional search using the chemical formula is done if the name
        did not already give a good match. Makes the search considerable slower.
    formula_search_depth: int
        How many of the most relevant formula matches to explore deeper. Default = 25.
    """
    if inchi is None:
        match_inchi = True
        mode = 'and'  # Do not allow 'or' in that case.
    else:
        match_inchi = False

    if inchikey is None:
        match_inchikey = True
        mode = 'and'  # Do not allow 'or' in that case.
    else:
        match_inchikey = False

    if mode == 'and':
        operate = operator.and_
    elif mode == 'or':
        operate = operator.or_
    else:
        print("Wrong mode was given!")

    inchi_pubchem = None
    inchikey_pubchem = None

    # Search pubmed for compound name:
    results_pubchem = pcp.get_compounds(compound_name,
                                        'name',
                                        listkey_count=name_search_depth)
    print("Found at least", len(results_pubchem),
          "compounds of that name on pubchem.")


    # Loop through first 'name_search_depth' results found on pubchem. Stop once first match is found.
    for result in results_pubchem:
        inchi_pubchem = '"' + result.inchi + '"'
        inchikey_pubchem = result.inchikey

        if inchi is not None:
            match_inchi = likely_inchi_match(
                inchi,
                inchi_pubchem,
                min_agreement=min_inchi_match)
        if inchikey is not None:
            match_inchikey = likely_inchikey_match(
                inchikey,
                inchikey_pubchem,
                min_agreement=min_inchikey_match)

        if operate(
                match_inchi, match_inchikey
                ):  # Found match for inchi and/or inchikey (depends on mode = 'and'/'or')
            print("--> FOUND MATCHING COMPOUND ON PUBCHEM.")
            if inchi is not None:
                print("Inchi ( input ): " + inchi)
                print("Inchi (pubchem): " + inchi_pubchem + "\n")
            if inchikey is not None:
                print("Inchikey ( input ): " + inchikey)
                print("Inchikey (pubchem): " + inchikey_pubchem + "\n")
            break

    if not operate(match_inchi, match_inchikey):
        if inchi is not None and formula_search:
            # Do additional search on Pubchem with the formula

            # Get formula from inchi
            inchi_parts = inchi.split('InChI=')[1].split('/')
            if len(inchi_parts) >= min_inchi_match:
                compound_formula = inchi_parts[1]

                # Search formula on Pubchem
                sids_pubchem = pcp.get_sids(compound_formula,
                                            'formula',
                                            listkey_count=formula_search_depth)
                print("Found at least", len(sids_pubchem),
                      "compounds with formula", compound_formula,
                      "on pubchem.")

                results_pubchem = []
                for sid in sids_pubchem:
                    result = pcp.Compound.from_cid(sid['CID'])
                    results_pubchem.append(result)

                for result in results_pubchem:
                    inchi_pubchem = '"' + result.inchi + '"'
                    inchikey_pubchem = result.inchikey

                    if inchi is not None:
                        match_inchi = likely_inchi_match(
                            inchi,
                            inchi_pubchem,
                            min_agreement=min_inchi_match)
                    if inchikey is not None:
                        match_inchikey = likely_inchikey_match(
                            inchikey,
                            inchikey_pubchem,
                            min_agreement=min_inchikey_match)
                    # Found match for inchi and/or inchikey (depends on mode = 'and'/'or')
                    if operate(match_inchi, match_inchikey):
                        print("--> FOUND MATCHING COMPOUND ON PUBCHEM.")
                        if inchi is not None:
                            print("Inchi ( input ): " + inchi)
                            print("Inchi (pubchem): " + inchi_pubchem + "\n")
                        if inchikey is not None:
                            print("Inchikey ( input ): " + inchikey)
                            print("Inchikey (pubchem): " + inchikey_pubchem +
                                  "\n")
                        break

    if not operate(match_inchi, match_inchikey):
        inchi_pubchem = None
        inchikey_pubchem = None

        if inchi is not None and inchikey is not None:
            print("No matches found for inchi", inchi, mode, " inchikey",
                  inchikey, "\n")
        elif inchikey is None:
            print("No matches found for inchi", inchi, "\n")
        else:
            print("No matches found for inchikey", inchikey, "\n")

    return inchi_pubchem, inchikey_pubchem
