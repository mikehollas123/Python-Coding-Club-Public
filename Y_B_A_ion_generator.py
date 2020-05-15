import numpy as np
import re
dict_atomic_mass = {
    # mono isotopic mass for each atom - can expand for isotopes
    "C": 12,
    "H": 1.007825,
    "N": 14.003074,
    "O": 15.994915,
    "S": 31.972072}

def mono_mass_calc(seq, charge=1):
    """
    Returns the monoisotopic mass for a peptide sequence in the form GAP(15.994915)GP(15.994915)DGNNG

    mass modifications can be placed after the amino acid code that is modified inside brackets (XXX.XXXXXX) up to 6 d.p.
    and to a max mass of 999.999999.For positive modifications a + symbol is not required
    however for negative modifications a - symbol must be before the mass

    The charge can be specified or will default as +1

    if the input sequence contains the previous and next residues (such as R.GAPQGVQGGK.G) they will be ignored

"""


    peptide_comp = np.array([0, 2, 0, 1, 0])
    for i in dict_aminoacid_composition:
        for match in re.finditer(i, seq):
            peptide_comp += dict_aminoacid_composition[i]
    v = peptide_comp[0]  # C
    w = peptide_comp[1]  # H
    x = peptide_comp[2]  # N
    y = peptide_comp[3]  # O
    z = peptide_comp[4]  # S

    main_mass_peak = v * dict_atomic_mass["C"] + w * dict_atomic_mass["H"] + x * \
                     dict_atomic_mass["N"] + y * dict_atomic_mass["O"] + z * \
                     dict_atomic_mass["S"]

    if charge != 0:
        main_mass_peak = float(((main_mass_peak + charge * 1.007276466) / charge))

    return float(main_mass_peak)

from proteomic_tools import dict_aminoacid_composition


def Y_B_A_ion_generator(seq,max_charge=3,no_A_ions = 0):

    """
    Returns a dictionary of B anf Y ions for a given peptide sequence (seq) for charges +1 up to +max-charge
    if max charge is not specified, the default is +3.

    the dictonary is in the form:
    dict['Y' or 'B'][charge][ion number] i.e. ion b4 +1 would be dict['B'][1][4].


    """

    if max_charge < 2:

        max_charge = 2
    if seq[1] == ".":
        seq = seq[2:-2]
    new_string = list(seq)

    for x in range(len(new_string)):
        if new_string[x] == "[":
            new_string.pop(x)
            new_string.insert(x, "(")


        elif new_string[x] == "]":
            new_string.pop(x)
            new_string.insert(x, ")")
    seq = "".join(new_string)

    b_peptides = [seq[:i+1]for i in range(0,len(seq)-1) if seq[i+1] in dict_aminoacid_composition]
    b_peptides.append(seq)

    dict_ions = {}
    dict_ions["B"] = {}

    for i in range(1,max_charge):

        dict_ions["B"][i] = {}
        for j in range(0,len(b_peptides)):

            dict_ions["B"][i][j+1] = mono_mass_calc(b_peptides[j],i) - (15.994915 + 2*1.007825)

    dict_ions["Y"] = {}
    y_peptides = [seq[i-1:]for i in range(1,len(seq)+1) if seq[i-1] in dict_aminoacid_composition ]
    for i in range(1,max_charge):
        dict_ions["Y"][i] = {}
        for j in range(0,len(y_peptides)):
            dict_ions["Y"][i][len(y_peptides)-j] = mono_mass_calc(y_peptides[j],i)

    if no_A_ions == 0:
        dict_ions["A"] = {}
        for i in range(1, max_charge):
            dict_ions["A"][i] = {}
            for j in range(0, len(b_peptides)):
                dict_ions["A"][i][j + 1] = mono_mass_calc(b_peptides[j], i) - (15.994915 + 2 * 1.007825) - 27.994915


    return dict_ions


if __name__ == '__main__':
    pep_seq = "ILDLGITGPEGHVLSRPEEVEAEAVNR"


    dict_ions = Y_B_A_ion_generator(pep_seq,4)


    Y_ions = [dict_ions['Y'][i][j]  for i in dict_ions['Y'] for j in dict_ions['Y'][i]]
    B_ions = [dict_ions['B'][i][j]  for i in dict_ions['Y'] for j in dict_ions['Y'][i]]
    A_ions = [dict_ions['A'][i][j]  for i in dict_ions['A'] for j in dict_ions['A'][i]]

    print (dict_ions['Y'])


