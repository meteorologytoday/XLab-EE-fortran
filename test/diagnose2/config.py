"""
 Programmer : Hsu, Tien-Yiao
 Last Update: 2015-08-03
 
 This program is the configuration file to generate test files for
 Eliassen-Sawyer operator solver.
	
 Test files include : gaussian field, sine x cosine field. They are
 all with constant A, C and with B = 0

"""
sys.path.append(os.path.realpath('../../xtt-lib-python'));

import numpy as np;
import eeenum;
import writeDiagConfig as WDC;

Lr = np.array([0.0, 1.0]); nr = 200;  Lr_span = Lr[1] - Lr[0];
Lz = np.array([0.0, 1.0]);  nz = 200; Lz_span = Lz[1] - Lz[0];

r_vec   = np.linspace(Lr[0], Lr[1], nr);
z_vec   = np.linspace(Lz[0], Lz[1], nz);

setting_opt = {
		'GEOMETRY' : eeenum.GEOMETRY.CYLINDRICAL,
		'DENSITY'  : eeenum.DENSITY.NORMAL,
		'OPERATOR_COMPLEXITY' : eeenum.OPERATOR_COMPLEXITY.BAROTROPIC,
		'POINTS' : {'horizontal': nr, 'vertical': nz},
		'DOMAIN_RANGE' : {'horizontal': Lr, 'vertical': Lz, 'planet_radius': 6371000.0},
		'INPUT_FOLDER':'.',
		'OUTPUT_FOLDER':'.',
		'A_FILE':'A.bin',
		'B_FILE':'B.bin',
		'C_FILE':'C.bin',
		'RESIDUE_CRITERIA' : {'absolute': 5e-3, 'relative': 5e-3, 'alpha': 1.0, 'max_iteration': 100000}
};

A_fld = np.zeros((nz, nr), dtype=np.float32);
B_fld = np.zeros((nz, nr), dtype=np.float32);
C_fld = np.zeros((nz, nr), dtype=np.float32);

A_fld[:,:] = 1.0;
C_fld[:,:] = 1.0;

for i in range(0, nr):
	for j in range(0, nz):			
		_r = r_vec[i];
		_z = z_vec[j];
		
		B_fld[j,i] = 1e-2 * np.sin( 2.0 * np.pi * ( _r - r_vec[0] ) / Lr_span ) * np.sin( 3.0 * np.pi * ( _z - z_vec[0] ) / Lz_span );

A_fld.tofile('%s/%s' % (setting_opt['INPUT_FOLDER'], setting_opt['A_FILE']));
B_fld.tofile('%s/%s' % (setting_opt['INPUT_FOLDER'], setting_opt['B_FILE']));
C_fld.tofile('%s/%s' % (setting_opt['INPUT_FOLDER'], setting_opt['C_FILE']));

WDC.writeDiagnoseFile('diag.txt', setting_opt);