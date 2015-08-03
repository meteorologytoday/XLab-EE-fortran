"""
 Programmer : Hsu, Tien-Yiao
 Last Update: 2015-08-03
 
 This program is the configuration file to generate test files for
 Eliassen-Sawyer operator solver.
	
 Test files include : gaussian field, sine x cosine field. They are
 all with constant A, C and with B = 0

"""
import sys;
import os;
sys.path.append(os.path.realpath('../../xtt-lib-python'));

import writeDiagConfig as WDC;

WDC.writeDiagnoseFile('./diag.txt', WDC.emptySetting());