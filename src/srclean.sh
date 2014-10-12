#!/bin/sh

rm *~ */*~
rm config.log config.status Makefile.in ./pirs/Makefile.in

rm Makefile depcomp configure missing aclocal.m4 INSTALL install-sh stamp-h1

rm -fr ./pirs/.deps ./autom4te.cache
