all: simulator stator symlinks

PERL_LIST = stator/baseCallingMatrix/baseCalling_Matrix_analyzer \
	stator/baseCallingMatrix/baseCalling_Matrix_calculator \
	stator/baseCallingMatrix/baseCalling_Matrix_merger \
	stator/alignment_stator \
	stator/gcContCvgBias/gc_coverage_bias_plot

simulator:
	cd ./src/pirs && ${MAKE}

stator:
	cd ./src/stator/gcContCvgBias && ${MAKE}

symlinks:
	-@ln -s ./src/pirs/pirs 2> /dev/null
	-@ln -s ./src/stator/gcContCvgBias/gc_coverage_bias 2> /dev/null
	-@ln -s ./src/stator/alignment_stator.pl alignment_stator 2> /dev/null
	-@for P in ${PERL_LIST}; do \
	    ln -s ./src/$${P}.pl `basename $${P}` 2> /dev/null; \
	done

clean:
	cd ./src/pirs && ${MAKE} clean
	cd ./src/stator/gcContCvgBias && ${MAKE} clean
	-rm pirs gc_coverage_bias
	-@for P in ${PERL_LIST}; do \
	    rm `basename $${P}`; \
	done
