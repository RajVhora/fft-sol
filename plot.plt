set ylabel "Normalised Stress"
set xlabel "Normalised Distance"
p "./numerical.dat" u 1:2 w p title"Numerical", "./analytical.dat" u 1:3 w l title"Analytical"
set term png
set output "11.png"
replot

p "./numerical.dat" u 1:3 w p title"Numerical", "./analytical.dat" u 1:4 w l title"Analytical"
set term png
set output "22.png"
replot

p "./numerical.dat" u 1:4 w p title"Numerical", "./analytical.dat" u 1:5 w l title"Analytical"
set term png
set output "12.png"
replot
