import os;

def readEfficiency(eff_file_path):
	eff = {};
	eff['semi_internal'] = 0.0;
	eff['semi_cb1'] = 0.0;
	eff['internal'] = 0.0;
	eff['semi_cb1'] = 0.0;
	eff['semi_cb1'] = 0.0;
	eff['wtheta'] = 0.0;
	eff['local_response'] = 0.0;
	if os.path.exists(eff_file_path):
		with open(eff_file_path, 'r') as eff_file:
			for line in eff_file:
				if   line.startswith(r' eta [L(B=0)    = 0]      w/  boundary'):
					eff['semi_internal'] += float((line.split(':')[1]).split(',')[1]);
				elif line.startswith(r' eta [L(B=0)    = dB]     wo/ boundary'):
					eff['semi_internal'] += float((line.split(':')[1]).split(',')[1]);
					eff['internal']+= float((line.split(':')[1]).split(',')[1]);
				elif line.startswith(r' eta [L(B=0)    = B0]     wo/ boundary'):
					eff['semi_internal'] += float((line.split(':')[1]).split(',')[1]);
					eff['internal'] += float((line.split(':')[1]).split(',')[1]);
				elif line.startswith(r' bndconv [L(B=0) = B0dB]   w/ boundary'):
					eff['semi_cb1'] += float((line.split(':')[1]).split(',')[1]);
				elif line.startswith(r' wtheta [L(B=0)    = J F] w/  boundary'):
					eff['wtheta'] += float((line.split(':')[1]).split(',')[1]);
				elif line.startswith(r' Local heat response (sum Q / sum dtheta_dt)'):
					eff['local_response'] += float(line.split(':')[1]);
		eff['semi_total'] = eff['semi_internal'] + eff['semi_cb1'];
		return eff;
	else:
		raise IOError();