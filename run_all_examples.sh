

echo run example 1:
echo python stattest.py -d examples/Isim_1.dat -f examples/fit_ellips.txt -k 5
python stattest.py -d examples/Isim_1.dat -f examples/fit_ellips.txt -k 5 -o examples/stattest_output_example1

echo run example 2:
echo python stattest.py -p examples -d Isim_1.dat -f "fit_sph.txt fit_ellips.txt fit_ellips_poly.txt" -k "4 5 6"
python stattest.py -p examples -d Isim_1.dat -f "fit_sph.txt fit_ellips.txt fit_ellips_poly.txt" -k "4 5 6" -o examples/stattest_output_example2
