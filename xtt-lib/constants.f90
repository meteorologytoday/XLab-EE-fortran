Module constants
implicit none

real(4), parameter :: g0=9.8, theta0=298.0, Rd=287.0,  Cv=5.0/2.0*Rd, Cp=Cv + Rd, kappa=Rd/Cp, &
                      h0=Cp*theta0/g0, p0 = 101300.0

real(4), parameter :: PI = acos(-1.0), PII = 2*PI

end module
