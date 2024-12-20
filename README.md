# StatTest version beta0.5
statistical test for SAXS data (applicable to other data formats). 

## Description and motivation
StatTest is made to make statistical tests in SAXS easily acessible. The script can be included in other software packages, or run post modelling.   
StatTest takes a dataset and one or more fit files as input. As output, you get:
- the chi-square, the reduced chi-square (should be close to one), and the corresponding p-value.    
- the number of runs, the reduced number of runs (should be close to one), and the corresponding p-value.
- the longest run, the reduced longest run (should be close to one), and the corresponding p-value.
- plot of the data along with plot of the normalized residuals.

The chi-square assumeed normal-distributed (Gaussian) errors. The runs test assumed equal probability of data being higher than or lower than the model fit (due to statistical variation).   

If multiple fits are given, F-test statistics are also given. I.e., a measure for how much better one model is compared to another. Typically used to assess if a more complex model (i.e. with more free parameters) is significantly better than a simpler model.

## Installation
Download stattest.py and stattest_functions.py to the same folder. The program should be run from the download folder.    
If h-tests are opted for, the folder htest should be in the same folder.   

### Requirements
* python3
* standard python packages (can be installed with `pip` or `conda`):
  * numpy
  * matplotlib
  * argparse
  * sys
  * scipy
  * os

## Running the program

```
python stattest.py -d <DATAFILE> -f <FITFILE(S)> -k <NUMBER-OF-FREE-PARAMETERS> (<OPTIONS>)
```

### Input options

#### required options:    
-d (or --data): name of data file. Should contain x-data, y-data and error on y-data (sigma).    
-f (or --fit): name of fit file. Should contain the fits with the same x-array as the datafile. It is assumed that the fit is in the second column, but this can be adjusted with flag -col. More fit files can be provided, see examples. 
-k (or --k): number of free parameters. For example, for sphere model with radius, scale and background, the number of free parameters (K) is three.    

#### other options:    
-htest (or --htest): show h test statistics     
-p (or --path): common path for data and fit file(s)       
-o (or --output_dir): provide name for output directory (default is stattest_output/)    
-logx (or --logx): plot data and fit on log-log plot (default is log-lin)    

#### for all options, type:
```
python stattest.py -h
```

## Examples
example data and fits are provided in the examples folder. Output files also given.     

### Example, single fit
Data of tri-axial ellipsoids were simulated with Shape2SAS and exported (Isim_1.dat) and fitted with a tri-axial ellipsoid model (5 parameters: scale, axis a, axis b, axis c, background) in SasView. The fitfile was likewize exported (fit_ellips.txt):    
```
python stattest.py -d examples/Isim_1.dat -f examples/fit_ellips.txt -k 5 -o examples/stattest_output_example1
```

### Example, multiple alternative fits
The same data were also fitted with a simpler model of polydispere spheres (4 parameters: scale, background, radius, polydispersity) and a more complex model of tri-axial ellipsoids with polydispersity in the one axis (axis a):    
```
python stattest.py -p examples -d Isim_1.dat -f "fit_sph.txt fit_ellips.txt fit_ellips_poly.txt" -k "4 5 6" -o examples/stattest_output_example2
```
Note that multiple inputs should be surrounded by quotation marks. Moreover, the number of free parameters should match the number of fits (these can be the same).   
The option -p (or --path) gives the opportunity to provide a pathf for all fits and data, to avoid typing it multiple times.    

### Run examples with bash shell script
The examples can be run using one of the bash shell scripts. In unix this can be done by: 
```
bash run_example1.sh
```
```
bash run_example2.sh
```
```
bash run_all_examples.sh
```

## Notes on h-tests
The script for calculating h statistics (in the folder htests) are adapted from:    
https://github.com/bio-phys/hplusminus/tree/main    
with minimal changes.    

## Acknowledgements
The program was written by Andreas Haahr Larsen    

## Relevant litterature
#### Longest runs statitics:    
https://maa.org/sites/default/files/pdf/upload_library/22/Polya/07468342.di020742.02p0021g.pdf    
https://www.nature.com/articles/nmeth.3358   

#### Number of runs statistics:    
https://en.wikipedia.org/wiki/Wald%E2%80%93Wolfowitz_runs_test    

#### h statistics: 
https://doi.org/10.26434/chemrxiv-2021-mdt29-v3 

