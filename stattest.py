#!/usr/bin/python3

"""
StatTest - statitiscal tests for SAXS or SANS data    
===================================================    

see description on the GitHub page: https://github.com/andreashlarsen/StatTest/tree/main          

"""

version = 'beta0.5' 

## importing python packages
import numpy as np
import argparse
import sys
import os

## import helper funcitons
try:
    from stattest_functions import *
except:
    printt("ERROR: stat-test tried to import functions from files stattest_functions.py")
    printt("this file should be in the same directory as stattest.py\n")
    from stattest_functions import *

if __name__ == "__main__":

    ## input values
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
    parser.add_argument("-p", "--path", help="Common path for data and fit file(s)", default='')
    parser.add_argument("-htest", "--htest", action="store_true", help="Apply h-tests", default=False)
    parser.add_argument("-no_plot", "--no_plot", action="store_true", help="Do not plot data and fit", default=False)
    parser.add_argument("-a", "--alpha", help="Significance criteria in F-test", default=0.01)
    parser.add_argument("-o", "--output_dir", help="Output directory", default="stattest_output")
    parser.add_argument("-kratky", "--kratky", action="store_true", help="Display Kratky plot",default=False)
    args = parser.parse_args()
    
    # import matplotlib.pyplot only if plotting enabled 
    if not args.no_plot:
        import matplotlib.pyplot as plt

    ## convert input 
    fit = convert_input(args.fit,False)
    k = convert_input(args.k,True)
    Nk = len(k)

    ## import data
    datafile = "%s/%s" % (args.path,args.data)
    if not os.path.exists(datafile):
         datafile = "%s%s" % (args.path,args.data)
    header,footer = get_header_footer(datafile)
    x,y,dy = np.genfromtxt(datafile,skip_header=header,skip_footer=footer,usecols=[0,1,2],unpack=True)
    N = len(x)

    ## read fit column input
    fit_column = int(args.fit_column)-1

    ## make output directory
    output_dir = args.output_dir
    CONTINUE,i = True,1
    while CONTINUE:
        if os.path.exists(output_dir):
            output_dir = '%s_%d' % (args.output_dir,i)
            i += 1
        else:
            CONTINUE = False 
    os.mkdir(output_dir)

    ## file for stdout using printt function
    f_out = open('%s/output.txt' % output_dir,'w')
    def printt(s):
        print(s)
        f_out.write('%s\n' %s)

    ## print welcome message
    printt('========================================================================================')
    printt('StatTest (version %s)' % version)
    printt('See description on the GitHub page: https://github.com/andreashlarsen/StatTest/tree/main')
    printt('========================================================================================')
    printt('data file: %s' % datafile)
    printt('Number of points: %d' % N)

    if not args.no_plot:
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
            printt("ERROR: currently, only up to 4 fits are supported at a time, contact developers if you need more.")
            sys.exit(-1)
        ax[0].errorbar(x,y,yerr=dy,marker='.',color='red',linestyle='none',label='data: %s' % args.data,zorder=0)

        if args.kratky:
            plt.figure(2)
            x2 = x**2
            plt.errorbar(x,x2*y,yerr=x2*dy,marker='.',color='red',linestyle='none',label='data: %s' % args.data,zorder=0)
            plt.figure(1)

    colors = ('black','green','blue','grey')
    chi2r_array,dof_array = np.zeros(Nk),np.zeros(Nk)

    for kk in range(Nk):

        fitfile =  "%s/%s" % (args.path,fit[kk])
        if not os.path.exists(fitfile):
            fitfile = "%s%s" % (args.path,fit[kk])

        K,color = k[kk],colors[kk]

        # import fit
        header,footer = get_header_footer(fitfile)
        yfit = np.genfromtxt(fitfile,skip_header=header,skip_footer=footer,usecols=[fit_column],unpack=True)

        # degrees of freedom
        DOF = N-K 
        dof_array[kk] = DOF

        ## normalized residuals
        R = (y-yfit)/dy 

        ## runs tests histogram
        h = get_runs_histogram(R)

        ## number of runs
        RN,RNr,RN_exp,RN_sigma,RN_p = get_RN(h,R,DOF)

        ## longest run
        RL,RLr,RL_exp,RL_sigma,RL_p = get_RL(h,DOF)

        ## chi2
        chi2,chi2r,chi2_sigma,chi2_p = get_chi2(R,DOF)
        chi2r_array[kk] = chi2r # used for model comparison with the F test

        printt('--------------------------------------------')
        printt('    ')
        printt('fit %d' % (kk+1))
        printt('    ')
        printt('    fitfile: %s' % fitfile)
        printt('    Number of free parameters: %d' % K)
        printt('    Degrees of freedom: %d' % DOF)
        printt('    ')
        printt('    chi2: %1.0f' % chi2)
        printt('    Expected chi2: %1.0f +- %1.0f' % (DOF,chi2_sigma))
        printt('    Reduced chi2: %1.1f' % chi2r)
        printt('    p-value chi2: %1.4f' % chi2_p)
        printt('    ')
        printt('    Longest run: %d' % RL)
        printt('    Expected Longest run: %1.1f +- %1.1f' % (RL_exp,RL_sigma))
        printt('    Reduced Longest run: %1.1f' % RLr)
        printt('    p-value Longest run: %1.4f' % RL_p)
        printt('    ')
        printt('    Number of runs: %d' % RN)
        printt('    Expected Number of runs: %1.0f +- %1.0f' % (RN_exp,RN_sigma))
        printt('    Reduced Number of runs: %1.1f' % RNr)
        printt('    p-value Number of run: %1.4f' % RN_p)
        printt('    ')

        if args.htest:
            printt('------------------------------ h test and hplusminus test ----------------------------------')
            from htest.evaluate import *
            from htest.io_hplus import *

            results = all_statistical_tests(R)
            print_pvalues_to_screen(results)
    
            printt('--------------------------------------------------------------------------------------------')

        if not args.no_plot:
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
            ax[kk+1].set_ylabel('%s' % args.reslabel)

            if args.kratky:
                plt.figure(2)
                plt.plot(x,x2*yfit,color=color,label=r'fit %d: %s, $\chi^2_r$: %1.1f, $R^L_r$: %1.1f, $R^N_r$: %1.1f' % (kk+1,fit[kk],chi2r,RLr,RNr),zorder=kk+1)
                plt.figure(1)

    kk_array = np.linspace(1,3,3)

    if Nk > 1:
        printt('------------------------------ F tests for model comparison ----------------------------------')
        printt('    ')
        idx = np.argsort(chi2r_array)
        chi2r_sort = chi2r_array[idx]
        dof_sort = dof_array[idx]
        kk_sort = kk_array[idx]

        for kkk in range(Nk-1):
            number1 = -1-kkk
            for kkkk in range(1,Nk-kkk):
                number2 = number1-kkkk              
                F0,F_p = get_F(chi2r_sort[number1],chi2r_sort[number2],dof_sort[number1],dof_sort[number2])
                printt('    compare fit %d and fit %d, F-test p-value %1.4f:' % (kk_sort[number1],kk_sort[number2],F_p))
                if F_p < args.alpha:
                    printt('        fit %d is significantly better than fit %d (p-value > significance level: %1.4f)' % (kk_sort[number2],kk_sort[number1],args.alpha))
                else:
                    printt('        fit %d is NOT significantly better than fit %d (p-value > significance level: %1.4f)' % (kk_sort[number2],kk_sort[number1],args.alpha))
                printt('    ')

    printt('Output send to output dir: %s' % output_dir)
    printt('    ')

    if not args.no_plot:
        if not args.nology:
                ax[0].set_yscale('log')
        if args.logx:
                ax[0].set_xscale('log')

        ax[kk+1].set_xlabel('%s' % args.xlabel)
        ax[0].set_ylabel('%s' % args.ylabel)
        ax[0].legend()
        fig.savefig('%s/stattest.pdf' % output_dir)

        if args.kratky:
            plt.figure(2)
            plt.ylabel('%s^2 %s' % (args.xlabel,args.ylabel))
            plt.xlabel('%s' % args.xlabel)
            plt.legend()
            plt.title('Kratky plot')
            plt.tight_layout()
            plt.savefig('%s/kratky.pdf' % output_dir)
            plt.figure(1)

        plt.show()

