1;

in_fld = "input";
out_fld = "output";
mode = 'TENDENCY-DENSITY_BOUSSINESQ-INTEGRAL_CHECK';
solver_strategy = 2;
solver_strategy_residue = 1e-5;
solver_max_iteration = 100000;


cases = {  'A',   'B',   'C',   'D',   'E', 'fig3ab', 'fig3cd'};
delta = [  0.0,   2.0,   4.0,  10.0,  24.0,      0.0,      5.0];
a =     [300.0, 250.0, 200.0, 150.0, 100.0,    200.0,    200.0] * 1000.0;
b1 =    [  0.0,   0.0,   0.0,   0.0,   0.0,      0.0,      0.0] * 1000.0;
b2 =    [300.0, 250.0, 200.0, 150.0, 100.0,    250.0,    250.0] * 1000.0;



theta0 = 298.0;
g0 = 9.8;
p0 = 100000.0;
R0 = 287.0;
rho0 = p0 / (R0*theta0);
kappa = 2.0/7.0;
Cp = R0 / kappa;
H0 = Cp * theta0 / g0;
zt = H0;

N_freq = 1.2e-2;
mu = (1000000.0)^(-1.0);

f_far  = mu * N_freq * zt / pi;
f_core = (1.0 + delta) * f_far;

dt = 10800.0;

QMax = (ones(1,length(cases)) * 10.0 / 86400.0 * (250000.0)^2.0) ./ (b2.^2.0 - b1.^2.0) * Cp;

Lr = [0 1000000.0]; nr = 1001;
Lz = [0        zt]; nz =  41;


r_vec   = linspace(Lr(1), Lr(2), nr);
z_vec   = linspace(Lz(1), Lz(2), nz);
rho_vec = z_vec .* 0 + rho0;

