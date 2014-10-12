all: simulator stator symlinks

PERL_LIST = stator/baseCallingMatrix/baseCalling_Matrix_analyzer \
	stator/baseCallingMatrix/baseCalling_Matrix_calculator \
	stator/baseCallingMatrix/baseCalling_Matrix_merger \
	stator/alignment_stator \
	stator/indelstat_sam_bam \
	stator/gcContCvgBias/gc_coverage_bias_plot

simulator:
	cd ./src/pirs && ${MAKE} -f gccMakefile

stator:
	cd ./src/stator/gcContCvgBias && ${MAKE}

symlinks:
	-@ln -s ./src/pirs/pirs 2> /dev/null
	-@ln -s ./src/stator/gcContCvgBias/gc_coverage_bias 2> /dev/null
	-@ln -s ./src/stator/alignment_stator.pl alignment_stator 2> /dev/null
	-@for P in ${PERL_LIST}; do \
	    ln -s ./src/$${P}.pl `basename $${P}` 2> /dev/null; \
	done

test: all
	cd ./src/pirs && ${MAKE} test
	cd ./src/stator/gcContCvgBias && ${MAKE} test

distclean:
	cd ./src/pirs && ${MAKE} distclean
	cd ./src/stator/gcContCvgBias && ${MAKE} distclean
	-rm pIRS_*.tgz

tDATE := $(shell date +%Y%m%d)
tTIME := $(shell date +%H%M%S)
dist: all distclean
	@echo "Packing pIRS_$(tDATE).tgz ..."
	@tar -czf /var/tmp/pIRS_$(tDATE)_$(tTIME).tgz --exclude '.git*' .
	@mv /var/tmp/pIRS_$(tDATE)_$(tTIME).tgz ./pIRS_$(tDATE).tgz

clean: distclean
	cd ./src/pirs && ${MAKE} clean
	cd ./src/stator/gcContCvgBias && ${MAKE} clean
	-rm pirs gc_coverage_bias
	-@for P in ${PERL_LIST}; do \
	    rm `basename $${P}`; \
	done

.PHONY push:
	git push github master
	git push
	git push google master

