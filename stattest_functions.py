import numpy as np
import scipy.stats as stats
from scipy.optimize import curve_fit
import os

def calculate_chi2r_p(reduced_chi_square,dof):
    """
    calculate p-value from reduced chi-square value
    """
    p_value = stats.chi2.sf(reduced_chi_square*dof, dof)
    return p_value

def calculate_chi2_p(chi_square,dof):
    """
    calculate p-value from chi-square value
    """
    p_value = stats.chi2.sf(chi_square, dof)
    return p_value

def get_header_footer(file):
    """
    get number of headerlines and footerlines of data file
    """

    header,footer = 0,0
    f = open(file)
    try:
        lines = f.readlines()
    except:
        print('Error: cannot read lines of file. Do you have some special characters in the file? Try removing them and rerun')
        print('file: %s' % file)

    CONTINUE_H,CONTINUE_F = True,True
    j = 0
    while CONTINUE_H or CONTINUE_F:
        line_h = lines[j]
        #print(line_h)
        line_f = lines[-1-j]
        tmp_h = line_h.split()
        tmp_f = line_f.split()
        try:
            NAN = 0
            for i in range(len(tmp_h)):
                1/float(tmp_h[i]) # divide to ensure non-zero values
                if np.isnan(float(tmp_h[i])):
                    NAN = 1
            if not tmp_h:
                NAN = 1 #empty line

            if NAN:
                header+=1
            else:
                CONTINUE_H = False
        except:
            header+=1
        try:
            NAN = 0
            for i in range(len(tmp_f)):
                1/float(tmp_f[i]) # divide to ensure non-zero values
                if np.isnan(float(tmp_f[i])):
                    NAN = 1
            if not tmp_h:
                NAN = 1 #empty line
                
            if NAN:
                footer+=1
            else:   
                CONTINUE_F = False
        except:
            footer+=1
        j+=1

    return header,footer