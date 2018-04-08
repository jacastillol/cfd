# Primer punto
set grid
samps = 10000
set samples samps
dx = 1./samps	
k1 = pi/(2*0.1)
k2 = 3./4.*k1
K1 = pi/(2*dx)
K2 = 3./4.*K1
plot [-pi:pi][] \
     exp(sin(x)) + 0.5*cos(K1*x) - 0.8*sin(K2*x) w l lw 0.2, \
     exp(sin(x)) + 0.5*cos(k1*x) - 0.8*sin(k2*x) w l lw 4
pause -1
set logscale y
plot [-pi:pi][] \
     cos(x)*exp(sin(x)) - 0.5*K1*sin(K1*x) - 0.8*K2*cos(K2*x) w l lw 0.2, \
     cos(x)*exp(sin(x)) - 0.5*k1*sin(k1*x) - 0.8*k2*cos(k2*x) w l lw 4
pause -1
plot [-pi:pi][] \
     (cos(x)**2-sin(x))*exp(sin(x)) - 0.5*K1**2*cos(K1*x) + 0.8*K2**2*sin(K2*x) w l lw 0.2, \
     (cos(x)**2-sin(x))*exp(sin(x)) - 0.5*k1**2*cos(k1*x) + 0.8*k2**2*sin(k2*x) w l lw 4
pause -1
#
set logscale xy
set format '10^{%T}'
set key right bottom
set title 'Error primera derivada vs $\Delta x$: Type real*8'
set linetype  1 lc rgb "dark-violet" lw 3 ps 1 pt 6
set linetype  2 lc rgb "#009e73" lw 3 ps 1 pt 12
set linetype  3 lc rgb "#56b4e9" lw 3 ps 1 pt 9
plot 'dat.dat' u (1/($1-1)):3 ls 1 t '2-CDS contorno', \
     'dat.dat' u (1/($1-1)):4 ls 2 t '4-CDS contorno', \
     'dat.dat' u (1/($1-1)):5 w l lw 3 t '6-CDS contorno', \
     'dat.dat' u (1/($1-1)):6 w l lw 3 t '2-CDS dominio', \
     'dat.dat' u (1/($1-1)):7 w l lw 3 t '4-CDS dominio', \
     'dat.dat' u (1/($1-1)):8 w l lw 3 t '6-CDS dominio'
pause -1
set title 'Error segunda derivada vs $\Delta x$: Type real*8'
plot 'dat.dat' u (1/($1-1)):9 w l lw 3 t '2-CDS contorno', \
     'dat.dat' u (1/($1-1)):10 w l lw 3 t '4-CDS contorno', \
     'dat.dat' u (1/($1-1)):11 w l lw 3 t '6-CDS contorno', \
     'dat.dat' u (1/($1-1)):12 w l lw 3 t '2-CDS dominio', \
     'dat.dat' u (1/($1-1)):13 w l lw 3 t '4-CDS dominio', \
     'dat.dat' u (1/($1-1)):14 w l lw 3 t '6-CDS dominio'
pause -1
