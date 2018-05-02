OBJ=obj
LIBS=-lblas -llapack
SRC=.

all: tarea01-a tarea01-b tarea01-c pdf/report.pdf

tarea01-a: $(SRC)/tarea01-dif.f90 $(OBJ)/diferencias_finitas.o $(OBJ)/parametros.o
	gfortran -I $(OBJ)/ -o tarea01-a tarea01-dif.f90 obj/*.o $(LIBS)

tarea01-b: $(SRC)/tarea01-dif-conv.f90 $(OBJ)/diferencias_finitas.o $(OBJ)/parametros.o
	gfortran -I $(OBJ)/ -o tarea01-b tarea01-dif-conv.f90 obj/*.o $(LIBS)

tarea01-c: $(SRC)/tarea01-dif-conv-trans.f90 $(OBJ)/diferencias_finitas.o $(OBJ)/parametros.o $(OBJ)/utilidades.o
	gfortran -I $(OBJ)/ -o tarea01-c tarea01-dif-conv-trans.f90 obj/*.o $(LIBS)

$(OBJ)/diferencias_finitas.o $(OBJ)/diferencias_finitas.mod: $(SRC)/diferencias_finitas.f90 $(OBJ)/parametros.o 
	gfortran -J $(OBJ)/ -c diferencias_finitas.f90 -o $(OBJ)/diferencias_finitas.o

$(OBJ)/parametros.o $(OBJ)/parametros.mod: $(SRC)/parametros.f90
	gfortran -J $(OBJ)/ -c $(SRC)/parametros.f90 -o $(OBJ)/parametros.o

$(OBJ)/utilidades.o $(OBJ)/utilidades.mod: $(SRC)/utilidades.f90
	gfortran -J $(OBJ)/ -c $(SRC)/utilidades.f90 -o $(OBJ)/utilidades.o

pdf/report.pdf: report.tex images/001Funcion.pdf images/002Errores_df.pdf images/003Errores_d2f.pdf images/004Errores_df_10dx.pdf images/005Regions_dx1.pdf images/005Regions_dx01.pdf images/006Regions_df2_dx01.pdf images/007Performance.pdf images/008SolucionTeorica.pdf images/009Solu_Sch_Grid.pdf images/010Solu_Sch_non-uGrid.pdf images/011Solu_Sch_non-uGrid-Good.pdf images/012Solu_Trans_FE_Grid11pts_11_1201.pdf images/013Solu_Trans_FE_nonUGrid11pts_11_1201.pdf images/014Solu_Trans_RK4_Grid_51_1201.pdf images/015Solu_Trans_pts_t.pdf images/Figura_EFvsRK4dt0_05.pdf images/Figura_nuEFvsnuRK4.pdf images/Figura_nuGrid.pdf images/Figura_nuEF_dt_re_x0.pdf images/Figura_EF_dt_dx.pdf images/Figura_RK4_dt_dx.pdf images/Figura_EF_pt_on_t.pdf images/Figura_RK4_pt_on_t.pdf
	pdflatex -output-directory pdf report.tex
	pdflatex -output-directory pdf report.tex

images/001Funcion.pdf: scripts/postprocess_punto_1_01.gp
	gnuplot scripts/postprocess_punto_1_01.gp
	pdflatex -output-directory images images/001Funcion.tex

images/002Errores_df.pdf: scripts/postprocess_punto_1_02.gp
	gnuplot scripts/postprocess_punto_1_02.gp
	pdflatex -output-directory images images/002Errores_df.tex

images/003Errores_d2f.pdf: scripts/postprocess_punto_1_03.gp
	gnuplot scripts/postprocess_punto_1_03.gp
	pdflatex -output-directory images images/003Errores_d2f.tex

images/004Errores_df_10dx.pdf: scripts/postprocess_punto_1_04.gp
	gnuplot scripts/postprocess_punto_1_04.gp
	pdflatex -output-directory images images/004Errores_df_10dx.tex

images/005Regions_dx01.pdf: scripts/postprocess_punto_1_05.gp
	gnuplot scripts/postprocess_punto_1_05.gp
	pdflatex -output-directory images images/005Regions_dx01.tex

images/005Regions_dx1.pdf: scripts/postprocess_punto_1_05a.gp
	gnuplot scripts/postprocess_punto_1_05a.gp
	pdflatex -output-directory images images/005Regions_dx1.tex

images/006Regions_df2_dx01.pdf: scripts/postprocess_punto_1_06.gp
	gnuplot scripts/postprocess_punto_1_06.gp
	pdflatex -output-directory images images/006Regions_df2_dx01.tex

images/007Performance.pdf: scripts/postprocess_punto_1_07.gp
	gnuplot scripts/postprocess_punto_1_07.gp
	pdflatex -output-directory images images/007Performance.tex

images/008SolucionTeorica.pdf: scripts/postprocess_punto_2_01.gp
	gnuplot scripts/postprocess_punto_2_01.gp
	pdflatex -output-directory images images/008SolucionTeorica.tex

images/009Solu_Sch_Grid.pdf: scripts/postprocess_punto_2_02.gp
	gnuplot scripts/postprocess_punto_2_02.gp
	pdflatex -output-directory images images/009Solu_Sch_Grid.tex

images/010Solu_Sch_non-uGrid.pdf: scripts/postprocess_punto_2_03.gp
	gnuplot scripts/postprocess_punto_2_03.gp
	pdflatex -output-directory images images/010Solu_Sch_non-uGrid.tex

images/011Solu_Sch_non-uGrid-Good.pdf: scripts/postprocess_punto_2_04.gp
	gnuplot scripts/postprocess_punto_2_04.gp
	pdflatex -output-directory images images/011Solu_Sch_non-uGrid-Good.tex

images/012Solu_Trans_FE_Grid11pts_11_1201.pdf: scripts/postprocess_punto_3_01.gp
	gnuplot scripts/postprocess_punto_3_01.gp
	pdflatex -output-directory images images/012Solu_Trans_FE_Grid11pts_11_1201.tex

images/013Solu_Trans_FE_nonUGrid11pts_11_1201.pdf: scripts/postprocess_punto_3_02.gp
	gnuplot scripts/postprocess_punto_3_02.gp
	pdflatex -output-directory images images/013Solu_Trans_FE_nonUGrid11pts_11_1201.tex

images/014Solu_Trans_RK4_Grid_51_1201.pdf: scripts/postprocess_punto_3_03.gp
	gnuplot scripts/postprocess_punto_3_03.gp
	pdflatex -output-directory images images/014Solu_Trans_RK4_Grid_51_1201.tex

images/015Solu_Trans_pts_t.pdf: scripts/postprocess_punto_3_04.gp
	gnuplot scripts/postprocess_punto_3_04.gp
	pdflatex -output-directory images images/015Solu_Trans_pts_t.tex

images/Figura_EFvsRK4dt0_05.pdf: scripts/postprocess_tmp01.gp
	gnuplot scripts/postprocess_tmp01.gp
	pdflatex -output-directory images images/Figura_EFvsRK4dt0_05.tex

images/Figura_nuEFvsnuRK4.pdf: scripts/postprocess_tmp02.gp
	gnuplot scripts/postprocess_tmp02.gp
	pdflatex -output-directory images images/Figura_nuEFvsnuRK4.tex

images/Figura_nuGrid.pdf: scripts/postprocess_tmp03.gp
	gnuplot scripts/postprocess_tmp03.gp
	pdflatex -output-directory images images/Figura_nuGrid.tex

images/Figura_nuEF_dt_re_x0.pdf: scripts/postprocess_tmp04.gp
	gnuplot scripts/postprocess_tmp04.gp
	pdflatex -output-directory images images/Figura_nuEF_dt_re_x0.tex

images/Figura_EF_dt_dx.pdf: scripts/postprocess_tmp05.gp
	gnuplot scripts/postprocess_tmp05.gp
	pdflatex -output-directory images images/Figura_EF_dt_dx.tex

images/Figura_RK4_dt_dx.pdf: scripts/postprocess_tmp06.gp
	gnuplot scripts/postprocess_tmp06.gp
	pdflatex -output-directory images images/Figura_RK4_dt_dx.tex

images/Figura_EF_pt_on_t.pdf: scripts/postprocess_tmp07.gp
	gnuplot scripts/postprocess_tmp07.gp
	pdflatex -output-directory images images/Figura_EF_pt_on_t.tex

images/Figura_RK4_pt_on_t.pdf: scripts/postprocess_tmp08.gp
	gnuplot scripts/postprocess_tmp08.gp
	pdflatex -output-directory images images/Figura_RK4_pt_on_t.tex

dirs:
	mkdir -p data
	mkdir -p images
	mkdir -p obj
	mkdir -p pdf
	mkdir -p scripts

clean:
	rm -f obj/*.o obj/*.mod pdf/*.aux pdf/*.log
	rm -f images/*.tex images/*.aux images/*.log
