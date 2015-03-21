function success = writeField(fname, data)
	[fd, msg]  = fopen(fname, "wb");
	if(fd == -1)
		error(['Open file error: ' msg]);
	endif

	cnt = fwrite(fd, data, "float32", 0 , "ieee-le");
	if(cnt != numel(data))
		error('Output binary data error!');
	endif
	fclose(fd);
endfunction
