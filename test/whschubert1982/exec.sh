#!/bin/bash

set +x

BINFLD=../../bin
DIAGPRO=old-diagnose.exe

python mkfield.py

${BINFLD}/${DIAGPRO} < input-150000_10/diagnose.txt
