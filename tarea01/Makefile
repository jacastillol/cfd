tarea01: tarea01.f90 diferencias_finitas.o
	gfortran -o tarea01 tarea01.f90 diferencias_finitas.o -lblas -llapack

diferencias_finitas.o diferencias_finitas.mod: diferencias_finitas.f90
	gfortran -c diferencias_finitas.f90

report.pdf: report.tex
	pdflatex report.tex
	pdflatex report.tex

clean:
	rm -f *.o *.mod *.aux *.log
