program simanhosc; implicit none

real(8), parameter :: twopi = 6.28318530717958647692528676655900577d0
real(8) kphys, kcom
logical break, change
integer j
real(8) dt1, y(7), yinit(7)

yinit(1) = 30.0
yinit(2) = 0.0
yinit(3) = hubb(yinit)
yinit(4) = 0.0
kphys = 1.e3*hubb(yinit)
dt1 = twopi*(1.0/(1.e3*yinit(3)))/60.0
change = .true.
break = .true.
y = yinit

do while (change)
    call gl8(y,dt1)
    if (y(4)>= 5.0) then
        yinit = y
        kcom = kphys*exp(y(4))
        yinit(5) = 1.0/sqrt(2.0*kcom)
        yinit(6) = 0.0
        yinit(7) = sqrt(kcom*kcom-zpp_over_z(yinit))
        change = .false.
    end if
end do

y = yinit
j = 0
do while (break)
    call gl8(y,dt1/20.0)
    j = j+1
    if (eps_inflation(y)>= 2.0) then
        break = .false.
    end if
    if (mod(j,1000)==0) then
        write(*,*) log(1/y(3)), y(4) 
    end if
end do

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! implicit Gauss-Legendre methods; symplectic with arbitrary Hamiltonian, A-stable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! equations of motion simple anharmonic oscillator

subroutine evalf(y, dydx)
        real(8) y(7), dydx(7)
! dx/dt = pi
        dydx(1) = y(2)
        dydx(2) = -3.d0*y(3)*y(2)-Vprime(y)
        dydx(3) = -0.5*y(2)**2
        dydx(4) = y(3)
        dydx(5) = y(6)*exp(-y(4))
        dydx(6) = -(kcom*kcom-zpp_over_z(y)-y(7)*y(7))*y(5)*exp(-y(4))
        dydx(7) = -2.0*(y(7)*y(6)/y(5))*exp(-y(4))
end subroutine evalf


! 10th order implicit Gauss-Legendre integrator
! only change n -> the number of equations
subroutine gl8(y, dt)
        integer, parameter :: s = 4, n = 7
        real(8) y(n), g(n,s), dt; integer i, k
        
        ! Butcher tableau for 8th order Gauss-Legendre method
        real(8), parameter :: a(s,s) = reshape((/ &
                 0.869637112843634643432659873054998518d-1, -0.266041800849987933133851304769531093d-1, &
                 0.126274626894047245150568805746180936d-1, -0.355514968579568315691098184956958860d-2, &
                 0.188118117499868071650685545087171160d0,   0.163036288715636535656734012694500148d0,  &
                -0.278804286024708952241511064189974107d-1,  0.673550059453815551539866908570375889d-2, &
                 0.167191921974188773171133305525295945d0,   0.353953006033743966537619131807997707d0,  &
                 0.163036288715636535656734012694500148d0,  -0.141906949311411429641535704761714564d-1, &
                 0.177482572254522611843442956460569292d0,   0.313445114741868346798411144814382203d0,  &
                 0.352676757516271864626853155865953406d0,   0.869637112843634643432659873054998518d-1 /), (/s,s/))
        real(8), parameter ::   b(s) = (/ &
                 0.173927422568726928686531974610999704d0,   0.326072577431273071313468025389000296d0,  &
                 0.326072577431273071313468025389000296d0,   0.173927422568726928686531974610999704d0  /)
        
        ! iterate trial steps
        g = 0.0; do k = 1,16
                g = matmul(g,a)
                do i = 1,s
                        call evalf(y + g(:,i)*dt, g(:,i))
                end do
        end do
        
        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl8


function Vphi(y)
real(8) y(7), Vphi
real(8), parameter :: lambda = 1.d-14
Vphi = lambda*y(1)**4/4.0
end function Vphi

function Vprime(y)
real(8) y(7), Vprime
real(8), parameter :: lambda = 1.d-14
Vprime = lambda*y(1)**3
end function Vprime

function hubb(y)
real(8) y(7), hubb
hubb = sqrt(y(2)**2/6.0+Vphi(y)/3.0)
end function hubb

function eps_inflation(y)
real(8) y(7), eps_inflation
eps_inflation = y(2)**2/(2.0*y(3)**2)
end function eps_inflation

function zpp_over_z(y)
real(8) y(7), epsilon1, epsilon2, epsilon3, zpp_over_z, hubbp, hubbpp, hubbppp
real(8) Vpp, ddot_phi, tdot_phi, eps2dot
Vpp = 3.0*Vprime(y)/y(1)
ddot_phi = -3.d0*y(3)*y(2)-Vprime(y)
epsilon1 = epsilon(y)
hubbp = -0.5*y(2)**2
tdot_phi = -3.0*hubbp*y(2)-3.0*y(3)*ddot_phi-Vpp*y(2)
hubbpp = -y(2)*ddot_phi
hubbppp = -ddot_phi**2-y(2)*tdot_phi
epsilon2 = hubbpp/(y(3)*hubbp)-2.0*hubbp/y(3)**2
eps2dot = hubbppp/(y(3)*hubbp)-(hubbpp*(hubbp**2+hubbpp*y(3)))/(y(3)*hubbp)**2-2.0*hubbpp/y(3)**2+4.0*hubbp**2/y(3)**3
epsilon3 = eps2dot/(y(3)*epsilon2)

zpp_over_z = exp(2.0*y(4))*y(3)**2*(2.d0-epsilon1+1.5*epsilon2+0.25*epsilon2**2-0.5*epsilon2*epsilon1+0.5*epsilon2*epsilon3)
end function zpp_over_z

end program simanhosc
