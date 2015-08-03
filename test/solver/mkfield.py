"""
 Programmer : Hsu, Tien-Yiao
 Last Update: 2014-11-18
 
 This program generate fields in the paper of Schubert and Hack (1982)
 "Inertial Stability and Tropical Cyclone Development".

 5 fields will be generated:
 
    Q, F, A, B, C

 Notice : This program uses row-major arrangement.
"""

import numpy as np;
import config as cf;
import os;

Q_fld = np.zeros((cf.nz-1, cf.nr-1), dtype=np.float32);
F_fld = np.zeros((cf.nz-1, cf.nr-1), dtype=np.float32);
A_fld = np.zeros((cf.nz, cf.nr), dtype=np.float32);
B_fld = np.zeros((cf.nz, cf.nr), dtype=np.float32);
C_fld = np.zeros((cf.nz, cf.nr), dtype=np.float32);

A_fld[:,:] = cf.N_freq ** 2.0;
B_fld[:,:] = 0.0;	
F_fld[:,:] = 0.0;

try:
	os.mkdir(cf.in_fld);
except OSError:
	print('Folder [%s] already exists.' % (cf.in_fld,));

try:
	os.mkdir(cf.out_fld);
except OSError:
	print('Folder [%s] already exists.' % (cf.out_fld,));


for (i,j), exp in np.ndenumerate(cf.exps):
	
	print("Producing (%d,%d)" % (i,j));
	
	in_fld  = exp['in_fld_case'];
	out_fld = exp['out_fld_case'];

	intl_stb_now = exp['region_intl_stb'];
	heating_now = exp['region_heating'];
	r_now = exp['region_r'];

	# Q field
	for i in range(0,cf.nr-1):
		for j in range(0,cf.nz-1):
			_r = (cf.r_vec[i] + cf.r_vec[i+1])/2.0;
			_z = (cf.z_vec[j] + cf.z_vec[j+1])/2.0;
			_pi = (cf.pi_vec[j] + cf.pi_vec[j+1])/2.0;
			
			if _r < r_now[0]:
				_J = heating_now[0];
			elif _r < r_now[1]:
				_J = heating_now[1];
			elif _r < r_now[2]:
				_J = heating_now[2];
			elif _r < r_now[3]:
				_J = heating_now[3];
			else:
				_J = heating_now[4];
			
			"""for _i in range(0, len(r_now)):
				if(_i < len(r_now) - 1 and _r < r_now[_i]):
					_J = heating_now[_i];
					break;
				elif(_i == len(r_now) - 1):
					_J = heating_now[_i+1];
			if(j==0 and _i != 2 and _i!=3):
				print(_i)"""
			
			Q_fld[j,i] = _J * cf.Cp * np.sin(np.pi * _z / cf.zt) * _pi;

	# A, C field (notice the grid size is different)
	for i in range(0,cf.nr):
		for j in range(0,cf.nz):			
			_r = cf.r_vec[i];
			_z = cf.z_vec[j];
			_rho = cf.rho_vec[j];

			if _r < r_now[0]:
				_c = intl_stb_now[0];
			elif _r < r_now[1]:
				_c = intl_stb_now[1];
			elif _r < r_now[2]:
				_c = intl_stb_now[2];
			elif _r < r_now[3]:
				_c = intl_stb_now[3];
			else:
				_c = intl_stb_now[4];

			C_fld[j,i] = _c;

	try:
		os.mkdir(in_fld);
	except OSError:
		print('Folder [%s] already exists.' % (in_fld,));
	
	try:
		os.mkdir(out_fld);
	except OSError:
		print('Folder [%s] already exists.' % (out_fld,));


	A_fld.tofile(in_fld + '/A.bin');
	B_fld.tofile(in_fld + '/B.bin');
	C_fld.tofile(in_fld + '/C.bin');
	Q_fld.tofile(in_fld + '/Q.bin');
	F_fld.tofile(in_fld + '/F.bin');

	try:
		### w   heating w/o pumping test
		diag_file = open('%s/%s' % (in_fld, 'diagnose.txt'), 'w');
		diag_file.write('%s          // mode\n'         % cf.mode);
		diag_file.write('%f          // delta time\n'   % cf.dt);
		diag_file.write('%f %f %f %f // domain size\n'  % (cf.Lr[0], cf.Lr[1], cf.Lz[0], cf.Lz[1]));
		diag_file.write('%d %d      // grid points\n'  % (cf.nr, cf.nz));
		diag_file.write('%s         // input folder\n' % in_fld);
		diag_file.write('%s         // output folder\n'% out_fld);
		diag_file.write('%s         // file: A\n' % 'A.bin');
		diag_file.write('%s         // file: B\n' % 'B.bin');
		diag_file.write('%s         // file: C\n' % 'C.bin');
		diag_file.write('%s         // file: Q\n' % 'Q.bin');
		diag_file.write('%s         // file: F\n' % 'F.bin');
		diag_file.write('%d %f %d %f // psi solver strategy and residue.\n' % (cf.solver_strategy[0], cf.solver_strategy_residue[0], cf.solver_max_iteration[0], cf.solver_alpha[0]));
		diag_file.write('%d %f %d %f // chi solver strategy and residue.\n' % (cf.solver_strategy[0], cf.solver_strategy_residue[0], cf.solver_max_iteration[0], cf.solver_alpha[0]));
		diag_file.write('no\nno\n');
	except IOError:
		print('Cannot open diagnose file.');
	finally:
		diag_file.close();
	
