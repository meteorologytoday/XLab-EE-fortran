module DIAGNOSE_DATA
implicit none

type DIAG_DATA
    integer :: nr,nz,nrA,nzA,nrB,nzB,nrC,nzC,nrQ,nzQ
    character(256) :: A_file, B_file, C_file, Q_file
end type

contains

subroutine setPoints(diag, nr, nz)
implicit none
type(DIAG_DATA), intent(inout) :: diag
integer, intent(in) :: nr, nz

diag%nr  = nr  ; diag%nz  = nz  ;
diag%nrA = nr-1; diag%nzA = nz-2;
diag%nrB = nr-1; diag%nzB = nz-1;
diag%nrC = nr-2; diag%nzC = nz-1;
diag%nrQ = nr-1; diag%nzQ = nz  ;

end subroutine

end module
