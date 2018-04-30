tarea01-a: tarea01-dif.f90 diferencias_finitas.o parametros.o
	gfortran -o tarea01-a tarea01-dif.f90 diferencias_finitas.o parametros.o -lblas -llapack

tarea01-b: tarea01-dif-conv.f90 diferencias_finitas.o parametros.o
	gfortran -o tarea01-b tarea01-dif-conv.f90 diferencias_finitas.o parametros.o -lblas -llapack

tarea01-c: tarea01-dif-conv-trans.f90 diferencias_finitas.o parametros.o utilidades.o
	gfortran -o tarea01-c tarea01-dif-conv-trans.f90 diferencias_finitas.o parametros.o utilidades.o -lblas -llapack

diferencias_finitas.o diferencias_finitas.mod: diferencias_finitas.f90 parametros.o 
	gfortran -c diferencias_finitas.f90

parametros.o parametros.mod: parametros.f90
	gfortran -c parametros.f90

utilidades.o utilidades.mod: utilidades.f90
	gfortran -c utilidades.f90

report.pdf: report.tex
	pdflatex report.tex
	pdflatex report.tex

clean:
	rm -f *.o *.mod *.aux *.log
