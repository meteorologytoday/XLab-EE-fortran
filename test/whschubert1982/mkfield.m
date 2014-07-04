1;

% This program generate fields in the paper of W. H. 
% Schubert et al.(1982) "Inertial stability and Tropical
% Cyclone Developement."
%
% Four thermal fields will be generated:
% 
%    Q, A, B, C   of case A, B, C, D, E, fig3ab, fig3cd
%
% Notice : This program uses row-major arrangement.
%


source('./config.m')

Q_fld = zeros(nr-1, nz-1);
A_fld = zeros(nr, nz);
B_fld = zeros(nr, nz);
C_fld = zeros(nr, nz);

A_fld(:,:) = N_freq^2.0;
B_fld(:,:) = 0.0;	

[status, msg, msgid] = mkdir(in_fld);
if(status != 1)
	printf("Error: (%d)%s\n", msgid, msg);
endif

for k = 1:length(cases)
	fprintf(stderr,"Now processing case %s...", cases{k});
	fflush(stderr);
	% Q field
	for i = 1:nr-1
		for j = 1:nz-1
			_r = (r_vec(i) + r_vec(i+1)) / 2.0;
			_z = (z_vec(j) + z_vec(j+1)) / 2.0;
			if ( _r > b1(k) && _r < b2(k) )
				Q_fld(i,j) = QMax(k) * sin(pi * (_z - Lz(1)) / (Lz(2) - Lz(1)));
			else
				Q_fld(i,j) = 0.0;
			endif
		endfor
	endfor	

	% C field
	for i = 1:nr
		for j = 1:nz
			_r = r_vec(i);
			if ( _r < a(k) )
				C_fld(i,j) = f_core(k)^2.0; 
			else
				C_fld(i,j) = f_far^2.0;
			endif
		endfor
	endfor	

	prefix = [cases{k} '_'];

	files.diag = [prefix 'diagnose.txt'];
	files.A = [prefix 'A.bin'];
	files.B = [prefix 'B.bin'];
	files.C = [prefix 'C.bin'];
	files.Q = [prefix 'Q.bin'];

	writeField([in_fld '/' files.A], A_fld);
	writeField([in_fld '/' files.B], B_fld);
	writeField([in_fld '/' files.C], C_fld);
	writeField([in_fld '/' files.Q], Q_fld);

	[fd msg] = fopen(files.diag,'w');
	if(fd == -1) 
		fprintf(stdout, "Error : %s\n", msg);
	endif

	out_fld_case = [out_fld '_' cases{k}];
	fprintf(fd,"%s\n%f\n", mode, dt);
	fprintf(fd,"%f %f\n", Lr(2) - Lr(1), Lz(2) - Lz(1));
	fprintf(fd,"%d %d\n", nr, nz);
	fprintf(fd,"%s\n", in_fld);
	fprintf(fd,"%s\n", out_fld_case);
	fprintf(fd,"%s\n", files.A);
	fprintf(fd,"%s\n", files.B);
	fprintf(fd,"%s\n", files.C);
	fprintf(fd,"%s\n", files.Q);
	fprintf(fd,"%d %f %d %f\n", solver_strategy(1), solver_strategy_residue(1), solver_max_iteration(1), solver_alpha(1));
	fprintf(fd,"%d %f %d %f\n", solver_strategy(2), solver_strategy_residue(2), solver_max_iteration(2), solver_alpha(2));
	fprintf(fd,"no\nno\n");
	fclose(fd);
	fprintf(stderr,"done.\n");
	fflush(stderr);

	[status, msg, msgid] = mkdir(out_fld_case);
	if(status != 1)
		printf("Error: (%d)%s\n", msgid, msg);
	endif
endfor


