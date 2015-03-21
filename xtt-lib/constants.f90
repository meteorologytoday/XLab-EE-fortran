Module constants
implicit none

real(4), parameter :: g0=9.8, theta0=298.0, Rd=287.0,  Cv=5.0/2.0*Rd, Cp=Cv + Rd, kappa=Rd/Cp, &
                      h0=Cp*theta0/g0, p0 = 101300.0

real(4), parameter :: PI = acos(-1.0), PII = 2*PI


contains

subroutine print_constants()
implicit none

print *, "Constants in xtt-lib/constants.f90 : "
print *, "Rd      (J/K/kg) : ", Rd
print *, "Cv      (J/K/kg) : ", Cv
print *, "Cp      (J/K/kg) : ", Cp
print *, "kappa            : ", kappa
print *, "g0      (m/s^2)  : ", g0
print *, "p0      (pa)     : ", p0
print *, "theta0  (K)      : ", theta0
print *, "h0      (m)      : ", h0


end subroutine





end module
