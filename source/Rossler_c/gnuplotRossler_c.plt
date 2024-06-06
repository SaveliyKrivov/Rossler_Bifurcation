 set autoscale
 set title "Rossler bifurcation (a = 0.2, b = 0.2, c = 1.0--8.0)" 

 plot "Rossler_c.d" using 1:2 with points pointtype 0 linecolor rgb "sea-green" notitle
 
 pause -1 "press OK" 
 #set output "Rossler_plot_c.pdf" 
 #set terminal postscript "Times-Roman" 20
 #set terminal pdf #monochrome
 replot
 reset 