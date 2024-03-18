# StatTest
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
Download stattest.py and stattest_functions.py to the same folder. The program should be run from this folder.

### Requirements
- python3
- packages: numpy, matplotlib, argparse, sys, scipy

## Running the program

```
python stattest.py <INPUT>
```

### Input options
required options:
-d: name of data. Should contain x-data, y-data and error on y-data (sigma).    
-f: name(s) of fit file(s). Should contain the fits with the same x-array as the datafile. It is assumed that the fit is in the second column, but this can be adjusted with flag -col.    
-k: number of free parameters (for each fit). For example, for sphere model with radius, scale and background, K is three.    
to see all options, type
```
python stattest.py -h
```

### Example, 1 fit: 
Data of tri-axial ellipsoids were simulated with Shape2SAS and exported (Isim_1.dat) and fitted with a tri-axial ellipsoid model (5 parameters: scale, axis a, axis b, axis c, background) in SasView. The fitfile was likewize exported (fit_ellips.txt):    
```
python stattest.py -d Isim_1.dat -f fit_ellips.txt -k 5
```

### Example, 1 fit: 
The same data were also fitted with a simpler model of polydispere spheres (4 parameters: scale, background, radius, polydispersity) and a more complex model of tri-axial ellipsoids with polydispersity in the one axis (axis a):    
```
python stattest.py -d Isim_1.dat -f "fit_sph.txt fit_ellips.txt fit_ellips_poly.txt" -k "4 5 6"
```

### Acknowledgements
The program was written by Andreas Haahr Larsen    

### Relevant litterature
to be written
