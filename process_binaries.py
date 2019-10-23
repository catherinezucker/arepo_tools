import numpy as np

##AREPO code units##
ulength=3.0856e+20 # 100 pc
umass=1.9910000e+33 # 1 solar mass
uvel=1.0e+5 #converting km/s to cm/s
udensity=umass/ulength**3 #needed to convert from native density units (Msun/100 pc^3) to cgs units
x_He=0.1 #Assumed helium abundance
mp=1.67e-24 # mass of proton in grams
n=500 #grid size

gridnum='x916y400' #grid to unpack
inpath = './' #assuming files are in current directory

#################Read in density#####################
fdensity = open(inpath+"density_grid_"+gridnum, "r")

d = np.fromfile(fdensity, dtype=np.float32)
d = np.delete(d, [0,1,2])
d = d.reshape([n,n,n])

#total gas density (in units of g/cm^3)
rhogas= d * udensity

print("Computed the total gas density, in units of g/cm^3")

fdensity.close()
#################Read in dust temperature#####################
fdusttemp = open(inpath+"tdust_grid_"+gridnum, "r")

d = np.fromfile(fdusttemp, dtype=np.float32)
d = np.delete(d, [0,1,2])
d = d.reshape([n,n,n])

#dust temperature, in units of K
dusttemp = d

print("Computed the dust temperature, in units of K")
fdusttemp.close()

#################Read in gas temperature#####################
fgastemp = open(inpath+"temp_grid_"+gridnum, "r")

d = np.fromfile(fgastemp, dtype=np.float32)
d = np.delete(d, [0,1,2])
d = d.reshape([n,n,n])

gastemp = d

print("Computed the gas temperature, in units of K")
fgastemp.close()

#################Read in CO Abundance#####################

faco = open(inpath+"xCO_grid_"+gridnum, "r")

d = np.fromfile(faco, dtype=np.float32)
d = np.delete(d, [0,1,2])
d = d.reshape([n,n,n])
aco = d

#number density of H-nucleons
# nHtot = nHI + nHp + 2*nH2
nHtot=rhogas/((1 + 4 * x_He) * mp)

#CO number density is the CO abundance times the number density of H-nucleons
nco_nden=aco*nHtot

print("Computed the CO number density, in units of cm^-3")

faco.close()

#################Read in H2 Abundance#####################

fah2 = open(inpath+"xH2_grid_"+gridnum, "r")

d = np.fromfile(fah2, dtype=np.float32)
d = np.delete(d, [0,1,2])
d = d.reshape([n,n,n])

ah2 = d

#number density of H-nucleons
nHtot=rhogas/((1 + 4 * x_He) * mp)

#H2 number density is the H2 abundance times the number density of H-nucleons
nh2_nden=ah2*nHtot

print("Computed the H2 number density, in units of cm^-3")

fah2.close()

#################Read in Velocities#####################

fvel = open(inpath+"velocity_grid_"+gridnum, "r")

d = np.fromfile(fvel, dtype=np.float32)
d = np.delete(d, [0,1,2])

velx=d[0::3]
vely=d[1::3]
velz=d[2::3]

#compute the x gas velocity, in units of cm/s
velx = velx.reshape([n,n,n]) * uvel

#compute the y gas velocity, in units of cm/s
vely = vely.reshape([n,n,n]) * uvel 
            
#compute the z gas velocity, in units of cm/s
velz = velz.reshape([n,n,n]) * uvel

print("Computed the vx, vy, vz velocities, in units of cm/s")
            
fvel.close()

print("All Done!")

