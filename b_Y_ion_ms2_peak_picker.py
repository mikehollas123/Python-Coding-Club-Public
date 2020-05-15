
from Y_B_A_ion_generator import Y_B_A_ion_generator
from matplotlib import pyplot




def ms2_matched_ions_plot(scan_data,peptide_seq,charge,show_plot=0,ppm=20.0,no_A_ions=0):

    expected_ions = Y_B_A_ion_generator(peptide_seq,charge,no_A_ions)
    ion_data = scan_data
    envelope_mass_diffs = [0.0]
    for i in range(1,4):
        envelope_mass_diffs.append(1.007276466/i)

    dict_found_ions = {}
    dict_found_ions = {}
    dict_found_ions["noise"] = {}
    dict_found_ions["ions"] = {}
    for i in range(len(ion_data)):
        for j in expected_ions:
            for k in expected_ions[j]:
                for l in expected_ions[j][k]:
                    if ion_data[i][0] > 50:
                        for m in range(len(envelope_mass_diffs)):
                            if ((expected_ions[j][k][l] - envelope_mass_diffs[m]) - ((expected_ions[j][k][l] - envelope_mass_diffs[m]) * (ppm / 1000000.0))) <= ion_data[i][0] <= ((expected_ions[j][k][l] - envelope_mass_diffs[m]) + ((expected_ions[j][k][l] - envelope_mass_diffs[m]) * (ppm / 1000000.0))):
                                dict_found_ions["ions"][ion_data[i][0]] = {}
                                dict_found_ions["ions"][ion_data[i][0]]["ion_type"] = "{0}_{1}+{3}".format(j, l, ion_data[i][0],k, (expected_ions[j][k][l]-ion_data[i][0]))
                                dict_found_ions["ions"][ion_data[i][0]]["error"] = ((ion_data[i][0] - (expected_ions[j][k][l] - envelope_mass_diffs[m])) * 1000000) / ion_data[i][0]
                                dict_found_ions["ions"][ion_data[i][0]]["int"] = ion_data[i][1]
                                dict_found_ions["ions"][ion_data[i][0]]["isotope"] = "-" + str(envelope_mass_diffs[m])

                            elif ((expected_ions[j][k][l] + envelope_mass_diffs[m]) - ((expected_ions[j][k][l] + envelope_mass_diffs[m]) * (ppm / 1000000.0))) <= ion_data[i][0] <= ((expected_ions[j][k][l] + envelope_mass_diffs[m]) + ((expected_ions[j][k][l] + envelope_mass_diffs[m]) * (ppm / 1000000.0))):
                                dict_found_ions["ions"][ion_data[i][0]] = {}
                                dict_found_ions["ions"][ion_data[i][0]]["error"] = ((ion_data[i][0] - (expected_ions[j][k][l] + envelope_mass_diffs[m])) * 1000000) / ion_data[i][0]
                                dict_found_ions["ions"][ion_data[i][0]]["int"] = ion_data[i][1]
                                dict_found_ions["ions"][ion_data[i][0]]["isotope"] = "+" + str(envelope_mass_diffs[m])
                                dict_found_ions["ions"][ion_data[i][0]]["ion_type"] = "{0}_{1}+{3}+{2}".format(j, l,  m,k, (expected_ions[j][k][l] -ion_data[i][0]))

        if ion_data[i][0] not in dict_found_ions["ions"]:
            dict_found_ions["noise"][ion_data[i][0]] = ion_data[i][1]

    dict_found_ions["found_ion_count"] = len([dict_found_ions["ions"][x] for x in dict_found_ions["ions"] if dict_found_ions["ions"][x]['isotope'] == '-0.0' ])
    dict_found_ions["expected_ion_count"] = len([expected_ions[i][j][k]  for i in expected_ions for j in expected_ions[i] for k in expected_ions[i][j] ])-2*(charge-1)

    if show_plot ==1:

        noise_X = [i for i in dict_found_ions["noise"]]
        noise_Y = [dict_found_ions["noise"][i] for i in dict_found_ions["noise"]]

        signal_X =[ i for i in dict_found_ions["ions"]]
        signal_Y = [dict_found_ions["ions"][i]['int'] for i in dict_found_ions["ions"]]
        norm_signal_Y = [(signal_Y[i]/max(signal_Y))*100 for i in range(len(signal_Y))]
        norm_noise_Y = [(noise_Y[i]/max(noise_Y))*100 for i in range(len(noise_Y))]
        pyplot.vlines(noise_X,0,norm_noise_Y,label="noise")
        pyplot.vlines(signal_X,0,norm_signal_Y,"r",label="signal")
        pyplot.title(peptide_seq)


        text = []
        for i in range(len(norm_signal_Y)):
            if len(dict_found_ions["ions"][signal_X[i]]["ion_type"].split("_")[1].split("+")) == 3:
                charge = dict_found_ions["ions"][signal_X[i]]["ion_type"].split("_")[1].split("+")[1]
                isotope = "[M+{}]".format(dict_found_ions["ions"][signal_X[i]]["ion_type"].split("_")[1].split("+")[2])
            else:
                isotope = ""
                charge = dict_found_ions["ions"][signal_X[i]]["ion_type"].split("_")[1].split("+")[1]
            ion = dict_found_ions["ions"][signal_X[i]]["ion_type"].split("_")[0]
            ion_no = dict_found_ions["ions"][signal_X[i]]["ion_type"].split("_")[1].split("+")[0]
            pyplot.text(signal_X[i], norm_signal_Y[i] ,"$" +ion + "_{" +ion_no+ "}^{+" + charge + "}$" + isotope)

            text.append(ion + "<sub>" +ion_no+ "</sub><sup>+" + charge + "</sup> " + isotope)
        pyplot.show()

    return dict_found_ions

if __name__ == '__main__':

    ms2_matched_ions_plot() # (scan_data,peptide_seq,charge,show_plot=0):