
1;

source('./config.m');

for i=1:length(cases)
	printf('Case %d : %s \n', i, cases{i});
endfor

draw_case = input('# Choose case: ');
printf('Case chosen %d, %s \n',draw_case, cases{draw_case});
out_fld_case = [out_fld '_' cases{draw_case}];
[rr,zz] = meshgrid(r_vec, z_vec);
[rrB,zzB] = meshgrid(r_vec(1:end-1),z_vec(1:end-1));
[rrA,zzA] = meshgrid(r_vec(1:end-1),z_vec(1:end));

figure(1);
subplot(2,2,1);
% draw figab figcd

data = readField([out_fld_case '/psi_before.bin'],nr,nz) ./ (1e13 / 86400.0);
[c,h] = contour(rr,zz,data,[0:1:100],'k-');
clabel_h = clabel(c,h,'rotation',0);
xlim([0 500000.0]);

subplot(2,2,2);
data = readField([out_fld_case '/dtheta_dt-B.bin'],nr-1,nz-1) * 86400.0;
[c,h] = contour(rrB,zzB,data);
clabel_h = clabel(c,h,'rotation',0);
xlim([0 500000.0]);


subplot(2,2,3);
data = readField([out_fld_case '/w-A.bin'],nr-1,nz);
[c,h] = contour(rrA,zzA,data);
clabel_h = clabel(c,h,'rotation',0);
xlim([0 500000.0]);

subplot(2,2,4);
data = readField([out_fld_case '/CHI-B.bin'],nr-1,nz);
[c,h] = contour(rrA,zzA,data);
clabel_h = clabel(c,h,'rotation',0);
xlim([0 500000.0]);

figure(2);%subplot(2,2,4);

data = readField([out_fld_case '/rhs_of_chi.bin'],nr,nz);
[c,h] = contour(rr,zz,data);
clabel_h = clabel(c,h,'rotation',0);
xlim([0 500000.0]);

% draw w
%{
hold on;
xlim([0 500000.0]);
ylim([0 z_vec(end)]);
rng=[0:0.01:0.15];
[c,h] = contourf(rrA,zzA,data, rng ,'linecolor','none');
colormap(jet(length(rng)-1));
caxis([rng(1) rng(end)]);
colorbar('EastOutside');

data = readField([out_fld_case '/wtheta-B.bin'],nr-1,nz-1);
[c,h] = contour(rrB,zzB,data);
clabel_h = clabel(c,h,'rotation',0);

%}


