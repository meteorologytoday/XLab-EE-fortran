open (unit=15,file=trim(output_folder)//"/result.txt",action="write",status="replace")
write (15,*) "Time elapsed (sec) : ", (time_end - time_beg)
close(15)
