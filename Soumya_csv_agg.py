import glob
from Examples_CSV_file_reader import csv_to_dict



file_dir = "C:\Data\Soumya"


files = glob.glob(file_dir + "\*.csv")


for file in files:

    dict_temp = csv_to_dict(file)
    try:
        split_filename = file.split("_")
        gender = split_filename[8]
        age = split_filename[9]
        sol = split_filename[10]
    except:
        print(split_filename)


    mono_masses = dict_temp["Monoisotopic_Mass"]

    dict_temp["gender"] = gender
    dict_temp["age"] = age
    dict_temp["sol"] = sol

    with open(file_dir + "\\combined.csv","a") as write_file:

        for i in range(len(dict_temp["Monoisotopic_Mass"])):
            mass = dict_temp['Monoisotopic_Mass'][i]
            write_file.writelines("{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9} \n ".format(dict_temp['Monoisotopic_Mass'][i],dict_temp['Sum_Intensity'][i],dict_temp['Monoisotopic_Mass'][i],dict_temp['Monoisotopic_Mass'][i],dict_temp['Monoisotopic_Mass'][i],dict_temp['Monoisotopic_Mass'][i],dict_temp['Monoisotopic_Mass'][i],dict_temp['Monoisotopic_Mass'][i],dict_temp['Monoisotopic_Mass'][i],dict_temp['Monoisotopic_Mass'][i]))












