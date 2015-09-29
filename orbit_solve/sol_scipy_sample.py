from numpy import *
from scipy.integrate import *
from scipy.optimize import *
from scipy.interpolate import *

# Parameters 
# sed editor run to replace parameters in run 

sR   = 0.08               #solenoid radius
sl   = SOL_LENGTH         #solenoid length
Lz   = 20.0               #simulation length of longitudinal direction
cl   = 2.99792458e8       #speed of light
ms   = 1.66053892e-27     #mass of proton
qc   = 1.60217657e-19     #charge of proton
en   = KINE_ENE            #kinetic energy of a particle
gmm  = qc*en/(ms*cl**2.0)+1 # Lorentz factor gamma
vz0   = cl*sqrt(1-1/gmm**2.0) # speed of a particle
bz0  = 1.7647058823529411 #scaled factor to set the center of magnetic field to 1 Tesla when sol_length=30cm
#bz0  = 1.41421 #scaled factor to set the center of magnetic field to 1 Tesla when sol_length=16cm
#bz0  = 0.894427 #scaled factor to set the center of magnetic field to 1 Tesla when sol_length=8cm
#bz0  = 0.485071 #scaled factor to set the center of magnetic field to 1 Tesla when sol_length=4cm
bfac = 1.0                # Center of magnetic field[T]
twpi = 2.0*pi             #two pi

#Integrand for calculation of Br(r,z)
def brf(th,r,z):
    global twpi,sR,sL,bz0,bfac
    return(bfac/twpi*(
        (sR*cos(th))/sqrt((z-sl/2.0)**2.0 + sR**2.0 + r**2.0 - 2.0*sR*r*cos(th))
       -(sR*cos(th))/sqrt((z+sl/2.0)**2.0 + sR**2.0 + r**2.0 - 2.0*sR*r*cos(th))
               )/bz0)

#Br(r,z)
def br(r,z):  
    return(quad(brf,0.0,2.0*pi,args=(r,z))[0])

#Integrand for calculation of Bz(r,z)
def bzf(th,r,z):
    global twpi,sR,sL,bz0
    return(bfac/twpi*
            (sR**2.0 - sR*r*cos(th))/(sR**2.0 + r**2.0 - 2.0*sR*r*cos(th))*
            ((z+sl/2.0)/sqrt((z+sl/2.0)**2.0 + sR**2.0 + r**2.0 - 2.0*sR*r*cos(th)) -
             (z-sl/2.0)/sqrt((z-sl/2.0)**2.0 + sR**2.0 + r**2.0 - 2.0*sR*r*cos(th))
            )/bz0)

#Bz(r,z)
def bz(r,z):  
    return(quad(bzf,0.0,2.0*pi,args=(r,z))[0])
    
#Integrand for calculation of A_theta = Atheta(r,z)
def atf(th,r,z):
    global twpi,sR,sL,bz0
    return(bfac*(r*sR**2.0)/twpi*
            (sin(th)**2.0)/(sR**2.0 + r**2.0 - 2.0*sR*r*cos(th))*
            ((z+sl/2.0)/sqrt((z+sl/2.0)**2.0 + sR**2.0 + r**2.0 - 2.0*sR*r*cos(th)) -
             (z-sl/2.0)/sqrt((z-sl/2.0)**2.0 + sR**2.0 + r**2.0 - 2.0*sR*r*cos(th))
            )/bz0)

#Atheta(r,z)
def at(r,z):  
    return(quad(atf,0.0,2.0*pi,args=(r,z))[0])

#function of calculating physical parameters of each step
def physpara(yvec):
    global ms,cl
    xx,xp,yy,yp,zz,zp = yvec                   #set particle condition
    rr = sqrt(xx**2.0 + yy**2.0)               #radius of beam
    beta = sqrt(xp**2.0+yp**2.0+zp**2.0)/cl    #Lorentz factor beta
    gamma = 1.0/sqrt(1.0-beta**2.0)            #Lorentz factor gamma
    pkin = ms*gamma*(xx*yp-yy*xp)              #kinetic term of ptheta
    ppot = qc*rr*at(rr,zz)                     #potential term of ptheta
    pthn = pkin + ppot                         #ptheta
    prr3 = ((pthn/ms/cl/gamma/beta)**2.0)/rr**3.0                     #ptheta term of radial kick
    fterm = (qc/ms/gamma)*(bz(rr,zz)/(rr * zp**2.0))*(xx*yp - yy*xp)  #first term of radial kick
    sterm = (qc/ms/gamma)*(br(rr,zz)/(rr**2 * zp**3.0))*(xx*xp + yy*yp)*(xx*yp - yy*xp) #second term of radial kick
    tterm = (xp**2 + yp**2)/(rr * zp**2.0) - (xx*xp + yy*yp)**2/(rr**3 * zp**2.0)       #third term of radial kick
    toradkick = fterm + sterm + tterm - prr3   #Radial kick for each step
    return(pthn,toradkick)
    
nt = NSTEP                        #total step
tvec = linspace(0,Lz/vz0,nt)      #time step of simulation
init = [INIT_AMP,0,0,INIT_VELO/gmm,-Lz/2.0,vz0] #initial values  [x, vx, y, vy, z, vz0]

#Coupled ordinary differential equations describing particle orbit to be integraded 
#  fr[6] = vector length 6 defining rhs of equations to integrate
#     fr[0] = d x/dt 
#     fr[1] = d^2 x/dt^2 
#     fr[2] = d y/dt 
#     fr[3] = d^2 y/dt^2 
#     fr[4] = d z/dt 
#     fr[5] = d^2 z/dt^2 
def deriv(t,yvec):
    global qc,gmm,ms
    fr = [0]*6
    xx,xp,yy,yp,zz,zp = yvec
    rr = sqrt(xx**2.0 + yy**2.0)
    fr[0] = xp
    fr[1] = (qc/ms/gmm)*( bz(rr,zz)*yp-br(rr,zz)*yy*zp/rr)
    fr[2] = yp
    fr[3] = (qc/ms/gmm)*(-bz(rr,zz)*xp+br(rr,zz)*xx*zp/rr)
    fr[4] = zp
    fr[5] = (qc/ms/gmm)*br(rr,zz)*(yy*xp - xx*yp)/rr
    return fr

#set integration conditions
ans=ode(deriv)
ans.set_integrator('dopri5') 
"""
 “dopri5”
 This is an explicit runge-kutta method of order (4)5 
 due to Dormand & Prince (with stepsize control and dense output).
 See also,
 http://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html
"""
ans.set_initial_value(init)
hst=empty((nt,6))
hst[0] = init #set initial condition of integration

#do the numerical integration
k=1
while ans.successful() and ans.t <  Lz/vz0 :
    ans.integrate(tvec[k])
    hst[k]=ans.y
    k += 1

#list of orbit from integration result
xx = hst[:,0] 
vx = hst[:,1]
yy = hst[:,2]
vy = hst[:,3]
zz = hst[:,4]
vz = hst[:,5]

#list of physical parameter
pt = array( [ physpara(hst[i]) for i in range(nt)])

#output data list 
# 0:x   
# 1:vx   
# 2:y   
# 3:vy   
# 4:z   
# 5:vz  
# 6:p_theta = mass*gamma*(xx*vy - yy*vx)+q*r*A_theta    Canonical angular momentum    
# 7:Radial kick for each simulation step   r' integral argument as a function of z 
savetxt("alldata.dat",transpose((xx,vx,yy,vy,zz,vz,pt[:,0],pt[:,1])))
