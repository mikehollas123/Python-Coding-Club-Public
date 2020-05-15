import glob
from Y_B_A_ion_generator import mono_mass_calc
import pickle
import numpy as np
from matplotlib import pyplot as plt

def generate():
    files = glob.glob("C:\Data\covid\\*.fasta")
    dict_seq = {}

    for file in files:
        with open(file,"r") as openfile:
            splitfiles = openfile.read().split(">")
            for protein in splitfiles[1:]:
                seq = protein.split("\n")[1]
                dict_seq[seq] = mono_mass_calc(seq,0)
    with open("C:\Data\covid\\x_outdict.pkl","wb") as p:
        pickle.dump(dict_seq,p)
    print(dict_seq)


def output(datafile="C:\Data\covid\\x_outdict.pkl"):

    with open(datafile, "rb") as pkl:
        dict_data = pickle.load(pkl)


    masses = [dict_data[x] for x in dict_data]

    binwidth = 1.0
    frq, edges = np.histogram(masses, np.arange(min(masses), max(masses) + binwidth, binwidth))

    unique = 0
    two_or_less = 0
    five_or_less = 0
    ten_or_less = 0
    greaterthan_ten = 0

    for i in frq:
        if i == 1:
            unique += 1
        if i <= 2:
            two_or_less += 1

        if i <= 5:
            five_or_less += 1

        if i <= 10:
            ten_or_less += 1
        if i > 10:
            greaterthan_ten += 1

    print("-----")
    print(unique)
    print(two_or_less)
    print(five_or_less)
    print(ten_or_less)
    print(greaterthan_ten)


    plt.hist(masses, np.arange(min(masses), max(masses) + binwidth, binwidth))
    plt.show()

if __name__ == '__main__':

    #generate()
    output()