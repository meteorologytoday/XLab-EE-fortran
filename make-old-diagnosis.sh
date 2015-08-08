#!/bin/bash

FC=gfortran
OUTFILE=old-diagnose

SRCFLD=./src/old-diagnose

mkdir ./bin

$FC ${SRCFLD}/diagnose.f90 ${SRCFLD}/xtt-lib/* -o ./bin/$OUTFILE
rm *.mod

