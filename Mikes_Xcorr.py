



from Y_B_A_ion_generator import Y_B_A_ion_generator
from matplotlib import pyplot as plt
import numpy as np
import math


def divide_chunks(l, n):
    """
    :param l: input list to be chunked
    :param n: the length of each chunk
    :return: the starting list normalized to 50 within each chunk
    """
    out_list = []
    # looping till length l
    for i in range(0, len(l), n):
        bin_list = l[i:i + n]

        bin_max = max(bin_list)
        if bin_max > 0:
            out_list = out_list + [(bin_list[x] / bin_max * 50) for x in range(len(bin_list))]
        else:
            out_list = out_list + [(0) for x in range(len(bin_list))]
    return out_list


def gaussian(x, amp, cen, wid):
    return amp * np.exp(-(x - cen) ** 2 / wid)

def Xcorr_peptide(scan_data,peptide,pre_mass,pre_charge=2,plot=0):

    """



    :param scan_data: input scan data in the form [(mass,int)]
    :param peptide: peptide sequence to compare to in the form GFP[+15.994915]GTSLPGPSGR
    :param pre_mass: mass of the precursor
    :param pre_charge: charge of the precursor
    :param plot: set to 1 to see the Xcorr plot
    :return: returns the Xcorr score for the peptide sequence with the scan specified
    """
    scan_data = sorted(scan_data)

    scan_no_pre = [scan_data[i] for i in range(len(scan_data)) if float(scan_data[i][0]) + 3 < pre_mass or float(scan_data[i][0]) - 3 > pre_mass if float(scan_data[i][0]) >0] #remove precursor

    X_gaus = np.linspace(200, (scan_no_pre[-1][0]+100), 4096)  # Create a continuous function with 4096 points

    comp_list_data = [gaussian(X_gaus, scan_no_pre[i][1], scan_no_pre[i][0], 0.05) for i in
                      range(len(scan_no_pre))]  # same here for the data
    master_list_data = [sum(x) for x in zip(*comp_list_data)]

    bin_size = int(math.ceil(float(len(master_list_data)) / 10))  # define the size of bins required for 10 (mostly) equal bins
    binned_data = list(divide_chunks(master_list_data, bin_size))  # place data into bins, and normalize to 50


    expected_ions = Y_B_A_ion_generator(peptide,pre_charge) # generate the expeceted ions a, b,and y

    Y_ions = [(expected_ions['Y'][1][i], 50 )for i in range(1,len(expected_ions['Y'][1])) if expected_ions['Y'][1][i] > 200] +  [(expected_ions['Y'][1][i]+1.006276746, 25 )for i in range(1,len(expected_ions['Y'][1])) if expected_ions['Y'][1][i] > 200]
    B_ions = [(expected_ions['B'][1][i], 50 )for i in range(1,len(expected_ions['B'][1])) if expected_ions['B'][1][i] > 200] +  [(expected_ions['B'][1][i]+1.006276746, 25 )for i in range(1,len(expected_ions['B'][1]))if expected_ions['B'][1][i] > 200]
    A_ions = [(expected_ions['A'][1][i], 10 )for i in range(1,len(expected_ions['A'][1])) if expected_ions['A'][1][i] > 200] + [(expected_ions['A'][1][i]+1.006276746, 10 )for i in range(1,len(expected_ions['A'][1])) if expected_ions['A'][1][i] > 200]

    """
    B and Y ions are given an intensity of 50, with the +1 isotopes given a intensity of 25 
    A ions are given an intensity of 10
    """
    all_ions = Y_ions + B_ions + A_ions
    all_ions.sort()



    comp_list_ions = [gaussian(X_gaus,all_ions[i][1],all_ions[i][0],0.05) for i in range(len(all_ions))] # create a list of all the gaussian modelled mass spec peaks
    master_list_ions =[sum(x) for x in zip(*comp_list_ions)] # and combine them into one spectrum


    xcorr_0 = np.correlate(master_list_ions,binned_data)/10000 # corrlation with no shift (/10000 is arbitrary lol)

    # correlation while shifting the ion spectrum over the data spectrum
    xcorr_1to75 = [(np.correlate(binned_data,[0]*i + master_list_ions[:-i]))/10000 for i in range(1,76)]
    xcorr_minus_1to75 = [(np.correlate(binned_data, master_list_ions[-i:]+[0]*-i))/10000 for i in range(-76,-1,1)]

    if plot ==1: # if enabled will show the Xcorr plot
        all = xcorr_minus_1to75 + [xcorr_0] + xcorr_1to75
        plt.vlines(range(-76,76,1),0,all)
        plt.show()

        plt.plot(X_gaus,binned_data)
        plt.plot(X_gaus,master_list_ions, "r")
        plt.show()

    return float(xcorr_0 - ((sum(xcorr_1to75)+sum(xcorr_minus_1to75))/150)) # returns the Xcorr score - the correlation at lag 0 - the mean correlation for lags -75 to +75 (not at 0)


if __name__ == '__main__':
    scan_data = [(394.2099, 1556.2), (326.1456, 898.1), (451.7497, 871.1), (884.4491, 4372.7), (207.4778, 893.6),
                 (226.1196, 4089.8), (210.1243, 786.5), (255.1462, 849.0), (283.1417, 4039.8), (616.3312, 1719.8),
                 (333.1782, 1892.4), (495.6889, 912.6), (956.48, 1280.0), (428.2158, 1823.6), (334.0955, 869.6),
                 (418.1049, 881.0), (920.8442, 867.1), (318.132, 1059.1), (208.109, 5429.6), (555.3085, 28018.8),
                 (543.6057, 971.6), (556.3107, 1561.4), (457.7762, 841.8), (337.1895, 1110.8), (526.7762, 1419.6),
                 (477.213, 998.6), (554.3516, 3528.8), (447.2299, 4493.9), (212.1038, 4464.6), (240.1356, 1379.9),
                 (446.2268, 42805.9), (300.1202, 6099.8), (730.3748, 1184.7), (955.4868, 1279.2), (422.5382, 810.5),
                 (559.3134, 1186.5), (411.238, 1327.8), (554.3073, 65517.0), (885.4547, 1437.9), (422.2055, 1047.5),
                 (828.4251, 1002.7), (555.1774, 878.0), (269.1264, 1533.6), (497.7496, 821.8), (456.9819, 767.5),
                 (956.9329, 1021.3), (258.1463, 1001.2), (827.4281, 2969.6), (206.4518, 782.0), (333.1782, 1892.4),
                 (446.2268, 42805.9), (884.4491, 4372.7), (885.4547, 1437.9), (956.48, 1280.0), (730.3748, 1184.7),
                 (828.4251, 1002.7), (226.1196, 4089.8), (447.2299, 4493.9), (616.3312, 1719.8), (283.1417, 4039.8),
                 (955.4868, 1279.2), (827.4281, 2969.6), (559.3134, 1186.5)]
    scan_data.sort()
    peptide = "GPAGPNGIP[+15.994915]GEK"

    pre_mass = 555.29
    pre_charge = 2
    print(Xcorr_peptide(scan_data,peptide,pre_mass,pre_charge,1))