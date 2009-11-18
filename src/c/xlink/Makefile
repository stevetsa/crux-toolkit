#Makefile
#

crux_source = "../../../../crux/src/c/"
#crux_source = "../c/"
crux_lib = "../../../../crux/src/c/.libs/libcrux.a"
#crux_lib = "../c/.libs/libcrux.a"

mst_lib="../../../../crux/src/c/MSToolkit/.libs/libmstoolkit.a"

LIBS=$(crux_lib) $(mst_lib)

objects = xhhc.o xhhc_scorer.o xhhc_ion_series.o crux_ion_series.o

program_objects = xhhc_predict_peptide_ions.o xhhc_score_peptide_spectrum.o xhhc_make_histogram.o xhhc_search.o

progs = xhhc-predict-peptide-ions xhhc-score-peptide-spectrum \
	xhhc-make-histogram xhhc-search


options = "-Wall"

# default target
all: $(program_objects) $(objects) $(progs)
	mkdir -p crux-output
#all: xhhc-predict-peptide-ions xhhc-score-peptide-spectrum \
#	xhhc-make-histogram xhhc-search


xhhc-predict-peptide-ions: xhhc_predict_peptide_ions.o
	g++ -I $(crux_source) $(options) -o xhhc-predict-peptide-ions \
	$(objects) xhhc_predict_peptide_ions.o $(LIBS)

xhhc-score-peptide-spectrum: xhhc_score_peptide_spectrum.o
	g++ -I $(crux_source) $(options) -o xhhc-score-peptide-spectrum \
	$(objects) xhhc_score_peptide_spectrum.o $(LIBS)

xhhc-make-histogram: xhhc_make_histogram.o
	g++ -I $(crux_source) $(options) -o xhhc-make-histogram \
	$(objects) xhhc_make_histogram.o $(LIBS)

xhhc-search: xhhc_search.o
	g++ -I $(crux_source) $(options) -o xhhc-search \
	$(objects) xhhc_search.o $(LIBS)

crux_ion_series.o:crux_ion_series.c crux_ion_series.h
	gcc -I $(crux_source) $(options) -c crux_ion_series.c 

xhhc.o: xhhc.cpp xhhc.h
	g++ -I $(crux_source) $(options) -c xhhc.cpp $(LIBS)

xhhc_scorer.o: xhhc_scorer.cpp xhhc_scorer.h xhhc.h
	g++ -I $(crux_source) $(options) -c xhhc_scorer.cpp $(LIBS)

xhhc_ion_series.o: xhhc_ion_series.cpp xhhc_ion_series.h xhhc.h xhhc_scorer.h
	g++ -I $(crux_source) $(options) -c xhhc_ion_series.cpp $(LIBS)

xhhc_predict_peptide_ions.o: xhhc_predict_peptide_ions.cpp $(objects)
	g++ -I $(crux_source) $(options) -c xhhc_predict_peptide_ions.cpp $(LIBS)

xhhc_score_peptide_spectrum.o: xhhc_score_peptide_spectrum.cpp $(objects)
	g++ -I $(crux_source) $(options) -c xhhc_score_peptide_spectrum.cpp $(LIBS)

xhhc_make_histogram.o: xhhc_make_histogram.cpp $(objects)
	g++ -I $(crux_source) $(options) -c xhhc_make_histogram.cpp $(LIBS)

xhhc_search.o: xhhc_search.cpp $(objects)
	g++ -I $(crux_source) $(options) -c xhhc_search.cpp $(LIBS)


clean:
	rm -rf $(progs) *.o 
