OBJ = gzstream.o self_util.o main.o stat_soap_coverage.o
CC = g++
LIBRARY = -O3 -lz -pipe -march=core2 -mtune=generic -pedantic -Wfloat-equal \
#-Wall

gc_coverage_bias:${OBJ}
	${CC} -o gc_coverage_bias ${LIBRARY} ${OBJ}
gzstream.o:gzstream.h gzstream.cpp
	${CC} -c ${LIBRARY} gzstream.h gzstream.cpp -static
self_util.o:self_util.cpp self_util.h
	${CC} -c ${LIBRARY} self_util.cpp self_util.h -static
main.o:main.cpp
	${CC} -c ${LIBRARY} main.cpp -static
stat_soap_coverage.o:stat_soap_coverage.h stat_soap_coverage.cpp
	${CC} -c ${LIBRARY} stat_soap_coverage.h stat_soap_coverage.cpp -static
test: gc_coverage_bias
	./gc_coverage_bias -r gcdeptestref.fa.gz -o test -w 100,200 --gcdump --depwindump gcdeptest.depth.gz
	-./gc_coverage_bias_plot.pl test_100.dat
	-./gc_coverage_bias_plot.pl test_200.dat
	ls -l test*
clean: justclean
	-rm -v test_*
justclean:
	-rm -v ${OBJ} *.gch

