 set autoscale
 set title "Rossler bifurcation (a = -2.0--4.0, b = 0.2, c = 5.7)" 
 set key noautotitle
 
 plot "Rossler_a.d" using 1:2 with points pointtype 0 linecolor rgb "sea-green" title ''
 
 pause -1 "press OK" 
 #set output "Rossler_plot_a.pdf" 
 #set terminal postscript "Times-Roman" 20
 #set terminal pdf #monochrome
 replot
 reset 