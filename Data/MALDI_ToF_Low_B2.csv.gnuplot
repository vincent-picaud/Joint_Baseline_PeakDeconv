set datafile separator ','
plot "Data/MALDI_ToF_Low_B2.csv" u 1 : 2 w l t "Y raw" lc "blue"
replot "Data/MALDI_ToF_Low_B2.csv" u 1 : 6 w l t "reconstructed" lw 2 lc "grey30"
replot "Data/MALDI_ToF_Low_B2.csv" u 1 : 4 w l t "baseline" lw 2 lc "grey70"
replot "Data/MALDI_ToF_Low_B2.csv" u 1 : 5 w l t "convolved peaks" lc "green"
replot "Data/MALDI_ToF_Low_B2.csv" u 1 : 3 w i t "peaks" lw 2 lc "red"
set terminal eps
set output "Data/MALDI_ToF_Low_B2.csv.eps"
replot
set terminal png
set output "Data/MALDI_ToF_Low_B2.csv.png"
replot
set terminal qt
set output
