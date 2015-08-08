"""
 Programmer : Hsu, Tien-Yiao
 Last Update: 2014-11-18
 
 This program is the configuration file of experiments in Schubert and Hack (1982)
 "Inertial Stability and Tropical Cyclone Development".

"""

import numpy as np;

mode = "CYLINDRICAL-TENDENCY-DENSITY_BOUSSINESQ-BAROTROPIC";
in_fld = "input";
out_fld = "output";

solver_strategy = np.array([2,2]);
solver_strategy_residue = np.array([3e-3,3e-3]);
solver_max_iteration = np.array([1000000,1000000]);
solver_alpha = np.array([1.0,1.0]);

theta0 = 298.0;
g0 = 9.8;
p0 = 100000.0;
R0 = 287.0;
kappa = 2.0/7.0;
Cp = R0 / kappa;
H0 = Cp * theta0 / g0;
zt = 5000.0 * np.pi;
dt = 3600.0;
N_freq = 1.2e-2;

radius_a  = 150000.0;
radius_b1 = 0;
radius_b2 = radius_a;

mu_variation_ratio = 10;
mu0 = 1000000.0**(-1.0);

f0 = N_freq * zt * mu0 / np.pi;
f_hat = (1.0 + mu_variation_ratio) * f0;

Q_hat = ((250000.0)**2.0 * 10.0 / 86400.0) / (radius_b2**2.0 - radius_b1**2.0);

Lr = np.array([0, 2000000.0]); nr = 1001;
Lz = np.array([0,zt]);  nz = 41;


r_vec   = np.linspace(Lr[0], Lr[1], nr);
z_vec   = np.linspace(Lz[0], Lz[1], nz);

konst = 1.0/4.0 * radius_a**4.0 * (f_hat**2.0 - f0**2.0);
def vfun(r):
	if( r < radius_a ):
		return r * (f_hat - f0) / 2.0;
	else:
		return  (1.0/4.0 * r**2.0 * f0**2.0 + konst / r**2.0)**(0.5) - f0 * r / 2.0;

wind_profile = [vfun(_r) for _r in r_vec];

in_fld += ("-%d_%d" % (radius_a, mu_variation_ratio));
out_fld += ("-%d_%d" % (radius_a, mu_variation_ratio));