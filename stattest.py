#!/usr/bin/python3

"""
StatTest - statitiscal tests for SAXS or SANS data    
===================================================    

version beta0.3     

see description on the GitHub page: https://github.com/andreashlarsen/StatTest/tree/main       

"""

## importing python packages
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
from scipy.stats import f

plot=1

## import helper funcitons
try:
    from stattest_functions import *
except:
    print("ERROR: stat-test tried to import functions from files stattest_functions.py")
    print("this file should be in the same directory as stattest.py\n")
    from stattest_functions import *

if __name__ == "__main__":

    ## input values
    #parser = argparse.ArgumentParser(description="""stat-test - statitiscal tests for SAXS or SANS data""",usage="python stattest.py -d data1.dat -f fit.dat <OPTIONAL ARGUMENTS>" )
    parser = argparse.ArgumentParser(description=__doc__,usage="python stattest.py -d data1.dat -f fit.dat <OPTIONAL ARGUMENTS>", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-d", "--data", help="Datafile. Include path and file extension. columns: x,y,sigma_y")
    parser.add_argument("-f", "--fit", help="fitfiles (format: fit1.dat OR \"fit1.dat fit.dat\" for multiple fits). Include path and file extension. columns: x,y_fit")
    parser.add_argument("-k", "--k", help="Number of free parameters (format: k1 OR \"k1 k2\" for multiple fits) (can be floats or integeres)")
    parser.add_argument("-logx", "--logx", action="store_true", help="Plot data with log-x axis", default=False)
    parser.add_argument("-nology", "--nology", action="store_true", help="Plot data with lin-y axis", default=False)
    parser.add_argument("-col", "--fit_column", help="column of fit in fitfile (integer, first column is 1, not 0)", default=2)
    parser.add_argument("-xlabel", "--xlabel", help="label on x-axis", default='q')
    parser.add_argument("-ylabel", "--ylabel", help="label on y-axis", default='I(q)')
    parser.add_argument("-rlabel", "--reslabel", help="label on residual y-axis", default='$\Delta I/\sigma$')
    parser.add_argument("-p", "--path", help="Common path for data and fit file(s)",default='')
    parser.add_argument("-htest", "--htest", action="store_true", help="Apply h-tests", default=False)
    args = parser.parse_args()
    
    ## convert fitfile string to list and remove empty entries
    fit_in = args.fit
    try:
        fit_tmp = fit_in.split(' ')      
        fit = []
        for i in range(len(fit_tmp)):
            if not fit_tmp[i] in ['',' ','  ','   ','    ','     ','      ','       ','        ']:
                fit.append(fit_tmp[i])
    except:
        print("ERROR: could not find fit. Try with option -f \"fit1.dat fit2.dat\"")
        sys.exit(1)

    k_in = args.k 
    try:
        k_tmp = k_in.split(' ')
        k = []
        for i in range(len(k_tmp)):
            if not k_tmp[i] in ['',' ','  ','   ','    ','     ','      ','       ','        ']:
                k.append(float(k_tmp[i]))
    except:
        print("ERROR: could not find k (number of free parameters fitted). Try with option -k k1 OR -k \"k1 k2\" for multiple fits")
        sys.exit(1)
    Nk = len(k)

    datafile = "%s/%s" % (args.path,args.data)
    if not os.path.exists(datafile):
         datafile = "%s%s" % (args.path,args.data)
    header,footer = get_header_footer(datafile)
    x,y,dy = np.genfromtxt(datafile,skip_header=header,skip_footer=footer,usecols=[0,1,2],unpack=True)
    N = len(x)

    print('----------------------')
        
    print('data file: %s' % datafile)
    print('Number of points: %d' % N)

    print('----------------------')

    figsize_x,figsize_y = 12,10
    if Nk == 1:
        fig,ax = plt.subplots(2,1,gridspec_kw={'height_ratios': [4,1]},figsize=(figsize_x,figsize_y))
    elif Nk == 2:
         fig,ax = plt.subplots(3,1,gridspec_kw={'height_ratios': [4,1,1]},figsize=(figsize_x,figsize_y))
    elif Nk == 3:
         fig,ax = plt.subplots(4,1,gridspec_kw={'height_ratios': [5,1,1,1]},figsize=(figsize_x,figsize_y))
    elif Nk == 4:
         fig,ax = plt.subplots(5,1,gridspec_kw={'height_ratios': [5,1,1,1,1]},figsize=(figsize_x,figsize_y))
    elif Nk > 4:
         print("ERROR: currently, only up to 4 fits are supported at a time, contact developers if you need more.")
         sys.exit(-1)
    ax[0].errorbar(x,y,yerr=dy,marker='.',color='red',linestyle='none',label='data: %s' % args.data,zorder=0)

    colors = ('black','green','blue','grey')

    fit_column = int(args.fit_column)-1
    chi2r_array,dof_array = np.zeros(Nk),np.zeros(Nk)

    for kk in range(Nk):

        fitfile =  "%s/%s" % (args.path,fit[kk])
        if not os.path.exists(fitfile):
            fitfile = "%s%s" % (args.path,fit[kk])

        K,color = k[kk],colors[kk]

        # import fit
        header,footer = get_header_footer(fitfile)
        yfit = np.genfromtxt(fitfile,skip_header=header,skip_footer=footer,usecols=[fit_column],unpack=True)


        DOF = N-K
        dof_array[kk] = DOF

        DY = y-yfit
        d_prev = DY[0]
        if d_prev > 0:
            Np,Nm = 1,0
        else:
            Np,Nm = 0,1
        
        run,h = 1,np.zeros(N)
        for d in DY[1:]:

            if d > 0:
                Np += 1
                if d_prev > 0:
                    run += 1
                else:
                    h[run] += 1
                    run = 1
            else:
                Nm += 1
                if d_prev < 0:
                    run += 1
                else:
                    h[run] += 1
                    run = 1
            d_prev = d
        h[run] += 1

        RN = np.sum(h)
        RL = np.max(np.where(h>0))
        
        RN_exp = 1+2*Np*Nm/DOF
        RNr = RN_exp/RN
        RL_exp = np.log2(DOF)-1
        RLr = RL/RL_exp

        RL_var = np.pi**2/6*np.log(2)**2+1/12
        RL_sigma = np.sqrt(RL_var)
        RL_res = (RL - RL_exp)/RL_sigma
        RL_p = 1-np.exp(-0.5**(RL - RL_exp+1))
        if(RL_p < 0.5):
                RL_p = 2.*RL_p
        else:
                RL_p = 2.*(1-RL_p)
        
        RN_var = (RN-1)*(RN-2)/(DOF-1)
        RN_sigma = np.sqrt(RN_var)
        RN_res = (RN - RN_exp)/RN_sigma

        xx = np.linspace(RN_exp-8*RN_sigma,RN_exp+8*RN_sigma,500)
        d = np.exp(-(RN_exp-xx)**2.0/RN_var/2.0)
        idx = np.where(xx<=RN)
        RN_p = np.sum(d[idx])/np.sum(d)
        if RN_p < 0.5:
            RN_p = 2.*RN_p
        else:
            RN_p = 2.*(1.-RN_p)

        R = DY/dy
        chi2 = np.sum(R**2)
        chi2r = chi2/DOF
        chi2r_array[kk] = chi2r
        chi2_var = 2*DOF
        chi2_sigma = np.sqrt(chi2_var)

        chi2_p = 2*calculate_chi2_p(chi2,DOF) # two-tailed: factor 2
        
        print('fit %d file: %s' % (kk+1,fitfile))
        print('----------------------')
        print('Number of free parameters: %d' % K)
        print('Degrees of freedom: %d' % DOF)

        print('---')

        print('chi2: %1.0f' % chi2)
        print('Expected chi2: %1.0f +- %1.0f' % (DOF,chi2_sigma))
        print('Reduced chi2: %1.1f' % chi2r)
        print('p-value chi2: %1.4f' % chi2_p)

        print('---')

        print('N+: %d' % Np)
        print('N-: %d' % Nm)

        print('---')

        print('Longest run: %d' % RL)
        print('Expected Longest run: %1.1f +- %1.1f' % (RL_exp,RL_sigma))
        print('Reduced Longest run: %1.1f' % RLr)
        print('p-value Longest run: %1.4f' % RL_p)

        print('---')

        print('Number of runs: %d' % RN)
        print('Expected Number of runs: %1.0f +- %1.0f' % (RN_exp,RN_sigma))
        print('Reduced Number of runs: %1.1f' % RNr)
        print('p-value Number of run: %1.4f' % RN_p)

        print('----------------------')

        if args.htest:
            print('------------------------------ h test and hplusminus test ----------------------------------')
            from htest.evaluate import *
            from htest.io_hplus import *

            results = all_statistical_tests(R)
            print_pvalues_to_screen(results)
    
            print('--------------------------------------------------------------------------------------------')

        Rmax = np.ceil(np.amax(abs(R)))
        Rmin = -Rmax

        ax[0].plot(x,yfit,color=color,label=r'fit %d: %s, $\chi^2_r$: %1.1f, $R^L_r$: %1.1f, $R^N_r$: %1.1f' % (kk+1,fit[kk],chi2r,RLr,RNr),zorder=kk+1)

        ax[kk+1].plot(x,R,marker='.',color='red',linestyle='none')
        ax[kk+1].plot(x,x*0,color=color)
        ax[kk+1].set_ylim(Rmin,Rmax)
        if Rmax > 4.5 and Rmax < 10:
            ax[kk+1].set_yticks([Rmin,-3,0,3,Rmax])
            ax[kk+1].plot(x,np.ones(N)*-3,color='grey',linestyle='--')
            ax[kk+1].plot(x,np.ones(N)*3,color='grey',linestyle='--')
        else:
            ax[kk+1].set_yticks([Rmin,0,Rmax])
        if args.logx:
            ax[kk+1].set_xscale('log')
        ax[kk+1].set_ylabel('%s' % args.reslabel) #r'$\Delta I/\sigma$'

    kk_array = np.linspace(1,3,3)
    if Nk > 1:
         idx = np.argsort(chi2r_array)
         chi2r_sort = chi2r_array[idx]
         dof_sort = dof_array[idx]
         kk_sort = kk_array[idx]
        
         for kkk in range(Nk-1):
            number1 = -1-kkk
            for kkkk in range(1,Nk-kkk):
                number2 = number1-kkkk
                F0 = chi2r_sort[number1]/chi2r_sort[number2]
                p = 1-f.cdf(F0,dof_sort[number1],dof_sort[number2]) #find p-value of F test statistic (one sided)
                p *= 2 # two sided p-value 
                print('compare fit %d and fit %d, F-test p-value: %1.4f' % (kk_sort[number1],kk_sort[number2],p))

    if not args.nology:
            ax[0].set_yscale('log')
    if args.logx:
            ax[0].set_xscale('log')

    ax[kk+1].set_xlabel('%s' % args.xlabel)
    ax[0].set_ylabel('%s' % args.ylabel)
    ax[0].legend()

    if plot:
        plt.show()
