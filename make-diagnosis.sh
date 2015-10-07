#!/bin/bash

FC=gfortran
OUTFILE=diagnose

SRCFLD=./src/diagnose
LIBFLD=./xtt-lib-fortran
mkdir ./bin

$FC ${SRCFLD}/main.f90 ${LIBFLD}/* -o ./bin/$OUTFILE
$FC ${SRCFLD}/main.f90 ${LIBFLD}/* -o ./bin/$OUTFILE

rm *.mod
