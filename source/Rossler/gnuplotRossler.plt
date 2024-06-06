set key  left  top  box 
set autoscale

HSVstart = 117/360.0
HSVend =  227/360.0  
# the HSV palette
set palette model HSV functions (1 - gray)*(HSVend - HSVstart) + HSVstart, 1, 0.68

set title "Rossler attractor (a = 0.2, b = 0.2, c = 5.7)" offset 0,1

splot "Rossler.d" every 1 using 2:3:4:1 with lines linecolor palette title ""

set ytics 1  # y scale marks and values will be at every ? unit
set grid      # grid turned on at each scale values

pause -1 "press OK" 
#set output "Rossler_plot.pdf" 
#set terminal pdf 
replot
reset 