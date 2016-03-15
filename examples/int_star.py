# execute on interactive command line
# %pylab
import sys
sys.path.append('../mypylib')
from stars import star_utils as su

# print su.grav_const
# G in units of 'cm^3 g^-1 s^-2'
# introduce code units:
# one unit of density: g_cu * cm_cu**-3 = 0.1691355
# one length unit:  1R_sun = 6.955e+10 cm
# one unit of mass: 1M_sun = 1.9891e+33 g
# one unit of time: 1yr    = 3.14e7s
L_cu = su.rsun_cm    # length in code unit
cm_cu = 1./L_cu 
M_cu = su.msun_g
g_cu = 1./M_cu
T_cu = 3.14e7
s_cu = 1/T_cu
G_cu = su.grav_const * cm_cu**3 * g_cu**-1 * s_cu**-2
# code unit of G: R_sun^3 Msun^-1 yr^-2
# one unit of pressure: [cgs: g cm^-1 s^-2] -> Msun Rsun^-1 yr^-2

rho_c = 0.0010887436239364599    # in cgs
rho_c = rho_c * g_cu * cm_cu**-3 # cu
P_c = 89386992764.462601         # in cgs
P_c = P_c * g_cu * cm_cu**-1 * s_cu**-2 # cu
print(rho_c, P_c)

# p2.get('radius')[-1],p2.get('mass')[-1]
r_c = 0.058674868915132176   # in L_cu = R_sun
m_c = 1.5612859321774683e-07 # in M_cu = M_sun

# EOS:
gamma_ad = 5./3.
K_ad = P_c/(rho_c**gamma_ad)
K    = K_ad
print(K)

def rho(p):   # polytropic EOS
        return (p/K)**(1./gamma_ad)

# unit tester EOS:
print(rho(P_c), rho_c)

# RHS of system of ODEs
def f_rhs(y,r):
    dm_dr = 4.*pi*(r**2)*rho(y[1])
    dp_dr = -rho(y[1])*G_cu*y[0]/(r**2)
    return [dm_dr,dp_dr]


# initial conditions
V0  = 4./3*pi*r_c**3
m_0 = rho_c*V0
print(m_0,m_c)

def int_eul(nsteps,r_c,dr,y0):
    r  = r_c
    dr = dr
    y  = array(y0)
    pp = []; mp = []; rp = [] 
    for i in range(nsteps):
        rhs = f_rhs(y,r)
        y += array(rhs)*dr
        r += dr
        mp.append(y[0])
        pp.append(y[1])
        rp.append(r)
    return mp, pp, rp

y0=[m_0,P_c]
dr = 0.01
steps = 2500
mass,pressure,radius = int_eul(steps,r_c,dr,y0)
plot(mass,radius)
