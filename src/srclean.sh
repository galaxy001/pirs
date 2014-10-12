#!/bin/sh

rm *~ */*~
rm config.log config.status Makefile.in ./pirs/Makefile.in

rm Makefile ./pirs/Makefile depcomp configure missing aclocal.m4 install-sh stamp-h1

rm -fr ./pirs/.deps ./autom4te.cache

find .|sort > afterclean.lst

diff -u cleansrc.lst afterclean.lst
