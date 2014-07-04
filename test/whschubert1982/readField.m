function success = readField(fname, nr, nz)
	num = nr*nz;
	[fd, msg]  = fopen(fname, "rb");
	if(fd == -1)
		error(['Open file error: ' msg]);
	endif

	[val, cnt] = fread(fd, num, "float32", 0 , "ieee-le");
	if(cnt != num )
		error('Read binary data error!');
	endif
	fclose(fd);

	success = reshape(val, nr, nz)';
endfunction
