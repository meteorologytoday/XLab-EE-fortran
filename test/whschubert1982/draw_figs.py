# -*- coding: UTF-8 -*-
import sys, os;
sys.path.append(os.path.realpath('../../xtt-lib-python'));

import numpy as np;

from matplotlib import rc;
import matplotlib.pyplot as pplt;
import matplotlib.font_manager as fm;
import matplotlib as mplt;

import config as cf;
import XWindProfile as wp;
from XContourExt import manualCLabelIfNotExists;



rc('font', **{'family':'sans-serif', 'serif': 'Bitstream Vera Serif', 'sans-serif': 'MS Reference Sans Serif', 'size': 20.0, 'weight' : 100})
rc('axes', **{'labelsize': 20.0, 'labelweight': 100})
rc('mathtext', **{'fontset':'stixsans'});

print("%s: %d"%(fm.FontProperties().get_name(),fm.FontProperties().get_weight()));

class FormattedFloat(float):
	def __str__(self):
		return '%.1f' % self.__float__();


r_vec = cf.r_vec.copy() / 1000.0;
z_vec = cf.z_vec.copy() / 1000.0;

rB_vec = (r_vec[0:len(r_vec)-1] + r_vec[1:len(r_vec)])/2; zB_vec = (z_vec[0:len(z_vec)-1]+z_vec[1:len(z_vec)])/2;
rA_vec = rB_vec; zA_vec = z_vec;
rC_vec = r_vec; zC_vec = zB_vec;

r_range = np.array([0, 500.0]);
z_range = cf.Lz / 1000;

# position of wind vector
wind_r_vec = np.linspace(5, 96, 10); 
wind_z_vec = np.linspace(1, 15, 7);
wind_threshold = 0.5;



# figure size in inches
figsize = [18,12];
figdpi=300;
border_linewidth = 2;
pplot_linewidth = 2;
avg_efficiency = 0.1;



try:
	Q = np.fromfile(cf.out_fld + '/J-B.bin',dtype='<f4', count=(cf.nr-1)*(cf.nz-1)).reshape((cf.nz-1,cf.nr-1));
	rpsi = np.fromfile(cf.out_fld + '/rpsi_before-O.bin',dtype='<f4', count=cf.nr*cf.nz).reshape((cf.nz,cf.nr));
	dtheta_dt = np.fromfile(cf.out_fld + '/dtheta_dt-B.bin',dtype='<f4', count=(cf.nr-1)*(cf.nz-1)).reshape((cf.nz-1,cf.nr-1));
	dtheta_dr = np.fromfile(cf.out_fld + '/RHS_rchi-O.bin',dtype='<f4', count=cf.nr*cf.nz).reshape((cf.nz,cf.nr)) * cf.theta0/cf.g0;
	
	rchi_0_B0 = np.fromfile(cf.out_fld + '/rchi-[0_B0]-O.bin',dtype='<f4', count=cf.nr*cf.nz).reshape((cf.nz,cf.nr));
	eta_0_B0 = np.fromfile(cf.out_fld + '/eta-[0_B0]-A.bin',dtype='<f4', count=(cf.nr-1)*cf.nz).reshape((cf.nz,cf.nr-1));
	rchi_0_dB = np.fromfile(cf.out_fld + '/rchi-[0_dB]-O.bin',dtype='<f4', count=cf.nr*cf.nz).reshape((cf.nz,cf.nr));
	eta_0_dB = np.fromfile(cf.out_fld + '/eta-[0_dB]-A.bin',dtype='<f4', count=(cf.nr-1)*cf.nz).reshape((cf.nz,cf.nr-1));

	rchi_0_B0dB = rchi_0_B0 + rchi_0_dB;
	eta_0_B0dB = eta_0_B0 + eta_0_dB;
	

except IOError as e:
	print(e);
	# Interpolate wind field to specified position	

### First graph: heating and efficiency comparison
fig, ax = pplt.subplots(2, 3, sharex = True, sharey = True, figsize=figsize);
fig.subplots_adjust(left=0.1, right=0.85, bottom=0.1, top=0.9, hspace=0.15, wspace=0.1);

for i in range(0,len(ax)):
	for j in range(0,len(ax[i])):
		for border in ax[i][j].spines:
			ax[i][j].spines[border].set_linewidth(border_linewidth);
		
		ax[i][j].set_xlim(r_range);
		ax[i][j].set_ylim(z_range);
		if(i==(len(ax)-1)):
			ax[i][j].set_xlabel(r'Radius [$\mathrm{km}$]');
		if(j==0):
			ax[i][j].set_ylabel(r'Pseudo Height [$\mathrm{km}$]');



for (i, text) in enumerate([
	r'(a) $v$ $[\mathrm{m}\,\mathrm{s}^{-1}]$',
	r'(b) $r \psi$ $[\times 10^{13}\,\mathrm{kg}\,\mathrm{day}^{-1}]$',
	r'(c) $\partial \theta / \partial t$ $[\mathrm{K}\,\mathrm{day}^{-1}]$',
	r'(d) $\partial \theta / \partial r$ $[\times 10^{-3}\,\mathrm{K}\,\mathrm{km}^{-1}]$',
	r'(e) $r \chi$ $[\times 10^{10}\,\mathrm{kg}]$',
	r'(f) $\eta$ $[\%]$']):
	idx = [int(i/3), i%3];
	ax[idx[0]][idx[1]].text(0.5, 1.01, text, transform=ax[idx[0]][idx[1]].transAxes, verticalalignment='bottom', horizontalalignment='center', color='black', fontsize=25);


# Heatin region

for i in range(0,len(ax)):
	for j in range(0,len(ax[i])):
		cmap = mplt.cm.hot_r;cmap.set_under(color='#ffffff', alpha=None);
		CS_heating = ax[i][j].contourf(rB_vec, zB_vec, Q * 86400.0, np.arange(0, 62, 6.0), cmap=cmap, antialiased=False, extend='both');

### Fig1 ###
# r'Heating region (red box) and heating efficiency ($10^{-1}\%$, contour)'
# Add contour

# Draw tangential wind
second = ax[0][0].twinx();
second.plot(r_vec, cf.wind_profile, '-k', linewidth=border_linewidth);
second.set_xlim(r_range);
second.set_ylim([0,60]);
second.set_ylabel(r'Tangential Wind $[\mathrm{m}\,\mathrm{s}^{-1}]$')
second.yaxis.set_label_coords(0.8, 0.5);
second.set_yticks([10, 20, 30, 40, 50]);
second.yaxis.set_tick_params(pad=-35);


### Fig2 ###
CS = ax[0][1].contour(r_vec, z_vec, rpsi * 86400.0 / 1e13, np.arange(0,10,1.0), colors='k', zorder=2, linewidths=border_linewidth);
CS.levels = [('%.1f' % val.__float__()) for val in CS.levels ];
#manualCLabelIfNotExists(ax[0][1], CS, 'rpsi');
ax[0][1].clabel(CS, inline=1, zorder=3);

### Fig3 ###
CS = ax[0][2].contour(rB_vec, zB_vec, dtheta_dt * 86400.0, np.arange(-500,500,2), colors='k', zorder=2, linewidths=border_linewidth);
CS.levels = [('%.0f' % val.__float__()) for val in CS.levels ];
#manualCLabelIfNotExists(ax[0][2], CS, 'dthetadt');
ax[0][2].clabel(CS, inline=1, zorder=3);

### Fig4 ###
CS = ax[1][0].contour(r_vec, z_vec, dtheta_dr * 1e6 , np.arange(-10,22,2), colors='k', zorder=2, linewidths=border_linewidth);
CS.levels = [('%.1f' % val.__float__()) for val in CS.levels ];
#manualCLabelIfNotExists(ax[1][0], CS, 'dthetadr');
ax[1][0].clabel(CS, inline=1, zorder=3);

### Fig5 ###
CS = ax[1][1].contour(r_vec, z_vec, rchi_0_B0dB / 1e10 , np.arange(0,200,10), colors='k', zorder=2, linewidths=border_linewidth);
CS.levels = [('%.0f' % val.__float__()) for val in CS.levels ];
#manualCLabelIfNotExists(ax[1][1], CS, 'rchi');
ax[1][1].clabel(CS, inline=1, zorder=3);

### Fig6 ###
CS = ax[1][2].contour(rA_vec, zA_vec, eta_0_B0dB * 100.0, np.arange(-10,20,0.1), colors='k', zorder=2, linewidths=border_linewidth);
CS.levels = [('%.2f' % val.__float__()) for val in CS.levels ];
#manualCLabelIfNotExists(ax[1][2], CS, 'eta');
ax[1][2].clabel(CS, inline=1, zorder=3);
ax[1][2].text(0.5, 0.5, r'$\overline{\eta} = ' + ('%.1f'% (avg_efficiency,)) + r'\%$', transform=ax[1][2].transAxes, fontsize=25);


cb = fig.colorbar(CS_heating, ticks=np.arange(0,66,6), cax=fig.add_axes([0.90, 0.1, 0.03, 0.8]), orientation='vertical');
cb.outline.set_linewidth(border_linewidth);
cb.ax.set_ylabel(r'$Q$ $[\mathrm{K}\,\mathrm{day}^{-1}]$', fontsize=25);

fig.savefig('exp-whs1982.png', dpi=figdpi);
