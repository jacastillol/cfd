set isosamples 40
unset key
set title "Burgers"
set xrange [0:1]
set yrange [0:1]
set ztics 1
unset surface
set pm3d
t=0.0
Re = 80
set zrange [0.5:1] 
fu(x,y) = (0.75-1.0/(4.0*( 1.0+exp( Re*(-t-4.0*x+4.0*y)/32.0) ) ))
fv(x,y) = (0.75+1.0/(4.0*( 1.0+exp( Re*(-t-4.0*x+4.0*y)/32.0) ) ))
splot fu(x,y)
set view 29,53 #Done implicitly by mousing.
#set term pngcairo mono enhanced
#set out 'bessel.png'
#replot
pause -1

splot fu(x,y)
pause -1
