import math
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from numpy import linspace


dict_aminoacid_composition = { # format = (C,H,N,O,S,C13,D)
    "A": (3,5,1,1,0,0,0),
    "R": (6,12,4,1,0,0,0),
    "N": (4,6,2,2,0,0,0),
    "D": (4,5,1,3,0,0,0),
    "C": (3,5,1,1,1,0,0),
    "E": (5,7,1,3,0,0,0),
    "Q": (5,8,2,2,0,0,0),
    "G": (2,3,1,1,0,0,0),
    "H": (6,7,3,1,0,0,0),
    "I": (6,11,1,1,0,0,0),
    "L": (6,11,1,1,0,0,0),
    "K": (6,12,2,1,0,0,0),
    "M": (5,9,1,1,1,0,0),
    "F": (9,9,1,1,0,0,0),
    "P": (5,7,1,1,0,0,0),
    "S": (3,5,1,2,0,0,0),
    "T": (4,7,1,2,0,0,0),
    "U": (5,5,1,2,0,0,0),
    "W": (11,10,2,1,0,0,0),
    "Y": (9,9,1,2,0,0,0),
    "V": (5,9,1,1,0,0,0)}

def Isotopic_Mass_pattern_generator(seq,charge=1,mods="",plot=0,line_width=0.01,print_output=0,max_only=0,mono_only=0):

    """
    This function takes a peptide sequence, the charge of the peptide, and any modifications (optional)
    and returns two lists of the same length: isotopic mass and relative intensity. yielding the isotopic
    distribution of a peptide. Using the various arguments described below, this data can be plotted for a visual
    representation of the distribution (either modelled using a gaussian lineshape or as stick representation of centroids).


    Arguments:
    seq - sequence of peptide using the codes in the dict_aminoacid_composition above
    charge - The charge of the peptide, defaults as 1
    mods - should be in the format of a string using the letters C,H,D,N,O,S and C[13] followed by the number of atoms (order does not matter)
    plot - determines if a plot of the isotopic distribution will be displayed. 0 = no plot, 1 = gaussian, 2 = centroid
    line_width - determines the line width used for the gaussian plot, only used if plot = 1
    print_output - if selected (print_output = 1) will print the determined chemical composition and the monoisotopic mass
    max_only - if selected (max_only = 1) the function will return a tuple containing the mass of the highest peak
    and number of peaks after the monoisotopic
    mono_only - if selected (mono_only = 1)the function will return the monoisotopic mass only
    """


    """Regex code to take the mods string and search for C12,C13,H,D,N,O,S atom count and add them to the array
"""
    m = re.search("C([\d]{1,3})", mods)
    if m != None:
        c12 = int(m.group(1))
    else:
        c12 = 0
    m = re.search("H([\d]{1,3})", mods)
    if m != None:
        H = int(m.group(1))
    else:
        H = 0
    m = re.search("N([\d]{1,3})", mods)
    if m != None:
        N = int(m.group(1))
    else:
        N = 0
    m = re.search("O([\d]{1,3})", mods)
    if m != None:
        O = int(m.group(1))
    else:
        O = 0
    m = re.search("S([\d]{1,3})", mods)
    if m != None:
        S = int(m.group(1))
    else:
        S = 0
    m = re.search("D([\d]{1,3})", mods)
    if m != None:
        D = int(m.group(1))
    else:
        D = 0
    m = re.search("C\[13\]([\d]{1,3})", mods)
    if m != None:
        c13 = int(m.group(1))
    else:
        c13 = 0

    """Create array of atoms in the format (C,H,N,O,S,C13,D)"""
    peptide_comp = np.array([c12, H + 2, N, O + 1, S,c13,D])  # start composition,  H2O + modification

    for i in seq:
        if i in dict_aminoacid_composition:
            peptide_comp = peptide_comp + np.array(
                dict_aminoacid_composition[i])  # add numpy arrays together to give total composition
        else:
            sys.exit("Unrecognised amino acid residue") # if unrecognized amino acid function exits with error



    v = peptide_comp[0] #C
    w = peptide_comp[1] #H
    x = peptide_comp[2] #N
    y = peptide_comp[3] #O
    z = peptide_comp[4] #S
    m = peptide_comp[5] #C13 from mods
    n = peptide_comp[6] #D from mods





    # probablities of isotopes
    PrC13 = 0.0107
    PrC12 = 1 - PrC13
    Ph1 = 0.99988
    Ph2 = 1 - Ph1
    Pn14 = 0.99636
    Pn15 = 1 - Pn14
    Po17 = 0.00038
    Po18 = 0.00205
    Po16 = 1 - Po17 - Po18  # 0.9976
    Ps32 = 0.9499
    Ps33 = 0.0075
    Ps34 = 0.0425
    Ps36 = 0.0001


    dict_atomic_abundance = { # (mono isotopic mass, natural abundance) in order of most common to least - some of this is redundant
        "C" : ((12 , 0.9893), (13.003355,0.0107)),
        "H" : ((1.007825 , 0.99988), (2.014102,0.00012)),
        "N" : ((14.003074 , 0.99636), (15.001090,0.00364)),
        "O" : ((15.994915 , 0.99757), (16.999131,0.00038),(17.999160,0.00205)),
        "S" : ((31.972071 , 0.9499), (32.971459,0.0075),(33.967867,0.0425),(35.967081,0.0001))}

    # mass differences from most abundant isotope for each isotope - may need to check


    add_H2 = dict_atomic_abundance["H"][1][0] - dict_atomic_abundance["H"][0][0]

    add_C13 = dict_atomic_abundance["C"][1][0] - dict_atomic_abundance["C"][0][0]

    add_N15 = dict_atomic_abundance["N"][1][0] - dict_atomic_abundance["N"][0][0]

    add_O17 = dict_atomic_abundance["O"][1][0] - dict_atomic_abundance["O"][0][0]
    add_O18 = dict_atomic_abundance["O"][2][0] - dict_atomic_abundance["O"][0][0]
    add_S33 = dict_atomic_abundance["S"][1][0] - dict_atomic_abundance["S"][0][0]
    add_S34 = dict_atomic_abundance["S"][2][0] - dict_atomic_abundance["S"][0][0]
    add_S36 = dict_atomic_abundance["S"][3][0] - dict_atomic_abundance["S"][0][0]
    #determine monoisotopic mass
    main_mass_peak = v * dict_atomic_abundance["C"][0][0] + w * dict_atomic_abundance["H"][0][0] + x * \
                         dict_atomic_abundance["N"][0][0] + y * dict_atomic_abundance["O"][0][0] + z * \
                         dict_atomic_abundance["S"][0][0] + m * dict_atomic_abundance["C"][1][0] + n * dict_atomic_abundance["H"][1][0]


    """
    Create isotopic patterns for each isotope type
    k = count of heavy atoms (0-15) 
    
    
    
    """

    number_of_iterations = 10


    C_peaks = [(k * add_C13,
                (math.factorial(v) / (math.factorial(k) * math.factorial(v - k))) * (PrC12 ** (v - k)) * (PrC13 ** k))
               for k in range(number_of_iterations) if (v - k) >= 0]


    H_peaks = [
        (k * add_H2, (math.factorial(w) / (math.factorial(k) * math.factorial(w - k))) * (Ph1 ** (w - k)) * (Ph2 ** k))
        for k in range(number_of_iterations) if (w - k) >= 0]

    N_peaks = [(k * add_N15,
                (math.factorial(x) / (math.factorial(k) * math.factorial(x - k))) * (Pn14 ** (x - k)) * (Pn15 ** k)) for
               k in range(number_of_iterations) if (x - k) >= 0]
    O_peaks = [
        (k * add_O17 + j * add_O18,
         (math.factorial(y) / (math.factorial(k) * math.factorial(j) * math.factorial(y - k - j))) * (
                     Po16 ** (y - k - j)) * (Po17 ** k) * (Po18 ** j)) for k in
        range(number_of_iterations) for j in range(number_of_iterations) if (y - k - j) >= 0]
    S_peaks = [
        (l * add_S36 + k * add_S33 + j * add_S34, (math.factorial(z) / (
                    math.factorial(l) * math.factorial(k) * math.factorial(j) * math.factorial(z - k - j - l))) * (
                     Ps32 ** (z - k - j)) * (Ps33 ** k) * (Ps34 ** j) * (Ps36 ** l)) for k in
        range(number_of_iterations) for j in range(number_of_iterations) for l in range(number_of_iterations) if z - k - j - l >= 0]


    """ 
    starting with the Carbon peaks, iterate through each peak and split each peak into another isotopic pattern.
    i.e. each carbon peak is split into the isotopic pattern of the H atoms, creating a new list of peaks, this list 
    is then split by the N atoms isotope pattern etc. this is ALOT of perumutations, so we ignore any values 
    below a threshold (0.0001 here)
    """

    peaks = [(C_peaks[i][0] + H_peaks[j][0], C_peaks[i][1] * H_peaks[j][1]) for j in range(len(H_peaks)) for i in
             range(len(C_peaks)) if C_peaks[i][1] * H_peaks[j][1] > 0.000001]

    peaks_2 = [(peaks[i][0] + N_peaks[j][0], peaks[i][1] * N_peaks[j][1]) for j in range(len(N_peaks)) for i in
               range(len(peaks)) if peaks[i][1] * N_peaks[j][1] > 0.000001]

    peaks_3 = [(peaks_2[i][0] + O_peaks[j][0], peaks_2[i][1] * O_peaks[j][1]) for j in range(len(O_peaks)) for i in
               range(len(peaks_2)) if peaks_2[i][1] * O_peaks[j][1] > 0.000001]
    peaks_4 = [(peaks_3[i][0] + S_peaks[j][0], peaks_3[i][1] * S_peaks[j][1]) for j in range(len(S_peaks)) for i in
               range(len(peaks_3)) if peaks_3[i][1] * S_peaks[j][1] > 0.000001]



    mean_mass_changes = [[(x[0],x[1])for x in peaks_4 if b-0.4 <= x[0] <= b+0.4] for b in range(9)]


    weighted_mass_changes = [sum([mean_mass_changes[l][o][0]*(mean_mass_changes[l][o][1]/sum([mean_mass_changes[l][x][1] for x in range(len(mean_mass_changes[l]))]))for o in range(len(mean_mass_changes[l]))])for l in range(len(mean_mass_changes))]


    """
    this leads to a lot of peaks, many are overlapping (or at least close) with each other i.e. +1 can come from 
    D,C13,N15,O17, S33 but with varying amounts. THese can be added together to create a centroid or all modelled
    with a gaussian and the curves added together. this creates the isotopic pattern. 
    """

    X = [(peaks_4[x][0] + main_mass_peak + charge * 1.007276466) / charge for x in range(len(peaks_4))] # add the monoisotopic mass and adjust for the charge

    Y = [peaks_4[x][1] for x in range(len(peaks_4))] # deconvelute peaks 4 into X and Y values for plotting
    Y = [(Y[x]/max(Y))*100 for x in range(len(Y))] # normalize Y values

    centroid_peaks = []
    for j in range(9):
        Y_values = 0.0
        X_values = []
        center = float(((main_mass_peak + charge * 1.007276466 + weighted_mass_changes[j]) / charge))

        for i in range(len(X)):

            if (center - (0.4/int(charge))) <= float(X[i]) <= (center + (0.4/int(charge))):
                Y_values += float(Y[i])
                X_values.append(float(X[i]))


        centroid_peaks.append((center, Y_values)) #


    cent_X = [centroid_peaks[x][0] for x in range(len(centroid_peaks))]
    cent_Y = [centroid_peaks[x][1] for x in range(len(centroid_peaks))]
    cent_Y = [ (cent_Y[x]/max(cent_Y))*100 for x in range(len(cent_Y)) ]

    if plot == 2: #centroid

        plt.title(seq + " +" + str(charge))
        plt.vlines(cent_X, 0,cent_Y)
        plt.xlabel("m/z")
        plt.ylabel("Relative Abundance")
        plt.show()


    if plot == 1: #gaussian
        def _1gaussian(x, amp1, cen1, sigma1):
            return amp1 * (1 / (sigma1 * (np.sqrt(2 * np.pi)))) * (np.exp(-((x - cen1) ** 2) / ((2 * sigma1) ** 2)))
        X_gaus = linspace(X[0]-5/charge, X[-1]+5/charge, 10000)
        comp_list = [_1gaussian(X_gaus, Y[i], X[i], line_width) for i in range(len(X))]

        master_list = [sum(x) for x in zip(*comp_list)]
        max_value = max([master_list[x] for x in range(len(master_list))])
        master_list = (master_list / max_value) * 100
        plt.title(seq + " +" + str(charge))
        plt.plot(X_gaus,master_list,"r")
        plt.xlabel("m/z")
        plt.ylabel("Relative Abundance")
        plt.show()


    if print_output == 1:
        print("peptide composition: C: " + str(peptide_comp[0]) + " H: " + str(peptide_comp[1]) + " N: " + str(
            peptide_comp[2]) + " O: " + str(peptide_comp[3]) + " S: " + str(peptide_comp[4]) + " D: " + str(
            peptide_comp[5]) + " C13: " + str(peptide_comp[6]))
        # print the peptide composition
        print( "monoisotopic mass: {0:.4f}".format((main_mass_peak + (charge) * 1.007276466) / charge))

    if max_only ==1:
        for i in range(len(cent_Y)):
            if cent_Y[i] == 100:
                return cent_X[i], i
    if mono_only ==1:
        return (main_mass_peak + (charge)*1.007276466)/charge

    return cent_X,cent_Y

if __name__ == '__main__':


    """
        open_file = open("peptides.js","r")
        lines = open_file.read().split("\n")
            
        start_time = time.clock()
        dict_peptides = {}
        for i in range(len(lines)):
            dict_peptides[lines[i]] = Isotopic_Mass_pattern_generator(lines[i],charge=1,plot=0,print_output=0)
    
        total_time = time.clock() - start_time
        average_time = total_time/len(lines)
    
        print total_time
        print average_time
        print 1/average_time
    
        print dict_peptides["RPGGWVEKETYY"]
    """

    seq = "LSDQCTGLQGFLVFHSFGGGTGSGFTSLLMER"
    print(Isotopic_Mass_pattern_generator(seq,charge=3,mods="C6H7O2N1", plot=1,print_output=1))







