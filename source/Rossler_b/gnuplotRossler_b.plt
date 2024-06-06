 set autoscale
 set title "Rossler bifurcation (a = 0.2, b = 0.0--2.0, c = 5.7)" 
 set key noautotitle
 
 plot "Rossler_b.d" using 1:2 with points pointtype 0 linecolor rgb "sea-green" title ''
 
 pause -1 "press OK" 
 #set output "Rossler_plot_b.pdf" 
 #set terminal postscript "Times-Roman" 20
 #set terminal pdf #monochrome
 replot
 reset 