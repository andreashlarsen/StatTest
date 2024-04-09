import numpy as np
import scipy.stats as stats

import sys

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

def convert_input(x_in,covert_to_float):
    """
    convert input string to list and remove empty entries
    """

    try:
        x_tmp = x_in.split(' ')      
        x = []
        for i in range(len(x_tmp)):
            if not x_tmp[i] in ['',' ','  ','   ','    ','     ','      ','       ','        ']:
                if covert_to_float:
                    x.append(float(x_tmp[i]))
                else:
                    x.append(x_tmp[i])
    except:
        print("ERROR: could not find input fit (fitfile) or k (number of free parameters). Try to:")
        print("    * add option -f fit.dat OR -f \"fit1.dat fit2.dat\" for multiple fits")
        print("    * add option -k k1 OR -k \"k1 k2\" for multiple fits")
        sys.exit(1) 

    return x 

def calculate_chi2_p(chi_square,DOF):
    """
    calculate p-value from chi-square value (one-tailed)
    """
    if chi_square > DOF:
        p_value = stats.chi2.sf(chi_square,DOF)
    else:
        p_value = stats.chi2.cdf(chi_square,DOF)
    
    return p_value

def get_chi2(R,DOF):
    """
    calculate chi2 statistics
    chi2: chi-square
    chi2r: reduced chi-square
    chi2_sigma: uncertainty on chi-square
    chi2_p: p value (two-tailed)
    """
    chi2 = np.sum(R**2)
    chi2r = chi2/DOF
    chi2_var = 2*DOF
    chi2_sigma = np.sqrt(chi2_var)
    chi2_p = 2*calculate_chi2_p(chi2,DOF) # two-tailed: factor 2

    return chi2,chi2r,chi2_sigma,chi2_p
    
def get_runs_histogram(R):
    """
    get histogram of run lengths (h)
    """
    R_prev = R[0]    
    run,h = 1,np.zeros(len(R))
    for d in R[1:]:
        if d > 0:
            if R_prev > 0:
                run += 1
            else:
                h[run] += 1
                run = 1
        else:
            if R_prev < 0:
                run += 1
            else:
                h[run] += 1
                run = 1
        R_prev = d
    h[run] += 1
    return h

def get_RN(h,R,DOF):
    """
    get number of runs
    Schilling runs test

    RN: number of runs
    RNr: reduced RN
    RN_exp: expected RN
    RN_sigma: uncertainty for RN
    RN_p: p-value for RN (two-tailed)
    """
    RN = np.sum(h)
    Np = len(np.where(R>0)[0])
    Nm = len(R)-Np
    RN_exp = 1+2*Np*Nm/DOF
    RNr = RN_exp/RN
    RN_var = (RN-1)*(RN-2)/(DOF-1)
    RN_sigma = np.sqrt(RN_var)
    xx = np.linspace(RN_exp-8*RN_sigma,RN_exp+8*RN_sigma,500)
    d = np.exp(-(RN_exp-xx)**2.0/RN_var/2.0)
    idx = np.where(xx<=RN)
    RN_p = np.sum(d[idx])/np.sum(d)
    if RN_p < 0.5:
        RN_p = 2.*RN_p
    else:
        RN_p = 2.*(1.-RN_p)     

    return  RN,RNr,RN_exp,RN_sigma,RN_p

def get_RL(h,DOF):
    """
    get longest run 
    Wald–Wolfowitz runs test

    RL: longest run
    RLr: reduced RL
    RL_exp: expected RL
    RL_sigma: uncertainty for RL
    RL_p: p-value for RL (two-tailed)
    """
    RL = np.max(np.where(h>0))
    RL_exp = np.log2(DOF)-1
    RLr = RL/RL_exp
    RL_var = np.pi**2/6*np.log(2)**2+1/12
    RL_sigma = np.sqrt(RL_var)
    RL_p = 1-np.exp(-0.5**(RL - RL_exp+1))
    if(RL_p < 0.5):
            RL_p = 2.*RL_p
    else:
            RL_p = 2.*(1-RL_p)

    return RL,RLr,RL_exp,RL_sigma,RL_p              

def get_F(chi2r1,chi2r2,dof1,dof2):
    """
    F-test for model comparison
    F0: F statistics
    F_p: associated p-value (two-tailed)
    """
    F0 = chi2r1/chi2r2
    F_p = 1-stats.f.cdf(F0,dof1,dof2) # find p-value of F statistic (one-tailed)
    F_p *= 2 # two tailed p-value 
    return F0,F_p