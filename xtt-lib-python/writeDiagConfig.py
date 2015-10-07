import eeenum;

def emptySetting():
	return {
	    'TYPE': eeenum.TYPE.DYNAMIC_EFFICIENCY,
		'GEOMETRY' : eeenum.GEOMETRY.CYLINDRICAL,
		'DENSITY'  : eeenum.DENSITY.NORMAL,
		'OPERATOR_COMPLEXITY' : eeenum.OPERATOR_COMPLEXITY.BAROTROPIC,
		'POINTS' : {'horizontal': 100, 'vertical': 100},
		'DOMAIN_RANGE' : {'horizontal': [0.0, 100.0], 'vertical': [0.0, 100.0], 'planet_radius': 6371000.0},
		'INPUT_FOLDER':'.',
		'OUTPUT_FOLDER':'.',
		'A_FILE':'A.bin',
		'B_FILE':'B.bin',
		'C_FILE':'C.bin',
		'FORCING_FILE': 'forcing.bin',
		'BC_INIT_FILE': 'bc_init.bin',
		'RESIDUE_CRITERIA' : {'absolute': 1e-5, 'relative':1e-3, 'alpha':1.0, 'max_iteration': 100000}
	};
	
	


def writeDiagnoseFile(file_path, setting):
	try:
		diag_file = open(file_path, 'w');
		diag_file.write('%s-%s-%s-%s   // geometry-density-operator_complexity\n' % (setting['TYPE'], setting['GEOMETRY'], setting['DENSITY'], setting['OPERATOR_COMPLEXITY']));
		if setting['GEOMETRY'] == eeenum.GEOMETRY.CYLINDRICAL:
			diag_file.write('%f %f %f %f // domain size\n'  % (setting['DOMAIN_RANGE']['horizontal'][0],
													 setting['DOMAIN_RANGE']['horizontal'][1],
													 setting['DOMAIN_RANGE']['vertical'][0],
													 setting['DOMAIN_RANGE']['vertical'][1]));
		elif setting['GEOMETRY'] == eeenum.GEOMETRY.SPHERICAL:	
			diag_file.write('%f %f %f %f // domain size\n'  % (setting['DOMAIN_RANGE']['planet_radius'],
													 setting['DOMAIN_RANGE']['vertical'][0],
													 setting['DOMAIN_RANGE']['vertical'][1]));
		diag_file.write('%d %d // grid points\n'  % (setting['POINTS']['horizontal'], setting['POINTS']['vertical']));
		diag_file.write('%s    // input folder\n' % setting['INPUT_FOLDER']);
		diag_file.write('%s    // output folder\n'% setting['OUTPUT_FOLDER']);
		diag_file.write('%s    // file: A\n' % setting['A_FILE']);
		diag_file.write('%s    // file: B\n' % setting['B_FILE']);
		diag_file.write('%s    // file: C\n' % setting['C_FILE']);
		diag_file.write('%s    // file: forcing\n' % setting['FORCING_FILE']);
		diag_file.write('%s    // file: boundary condition and initial guess\n' % setting['BC_INIT_FILE']);
		diag_file.write('%f %f %d %f // rchi solver residue absolute, residue relative, max iteration time, and alpha.\n' % (setting['RESIDUE_CRITERIA']['absolute'],setting['RESIDUE_CRITERIA']['relative'],setting['RESIDUE_CRITERIA']['max_iteration'],setting['RESIDUE_CRITERIA']['alpha']));
	except IOError:
		print('Cannot open diagnose file.');
	finally:
		diag_file.close();