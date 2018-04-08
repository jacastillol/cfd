# Segundo punto
set grid
Pe_1_0 = -10.0
Pe0_0 = 0.001
Pe1_0 = 10.0
Pe = 60.0
phiL = 1.0
phi0 = 0.0
L = 1.0
plot [0:phiL][0:L] \
     phi0 + (exp(Pe_1_0*x/L)-1)/(exp(Pe_1_0)-1)*(phiL-phi0) lw 3 t 'Pe =-10.0', \
     phi0 + (exp(Pe0_0*x/L)-1)/(exp(Pe0_0)-1)*(phiL-phi0) lw 3 t 'Pe = 0.0', \
     phi0 + (exp(Pe1_0*x/L)-1)/(exp(Pe1_0)-1)*(phiL-phi0) lw 3 t 'Pe = 10.0', \
     phi0 + (exp(Pe*x/L)-1)/(exp(Pe)-1)*(phiL-phi0) lw 3 t 'Pe = 60.0'
pause -1