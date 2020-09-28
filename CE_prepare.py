import numpy as np

# Read in header of GGchem 2D file
file   = 'Static_Conc_2D.dat'
data   = open(file)
dummy  = data.readline()
dimens = data.readline()
dimens = np.array(dimens.split())
NELEM  = int(dimens[0])
NMOLE  = int(dimens[1])
NDUST  = int(dimens[2])
NPOINT = int(dimens[3])
header = data.readline()
data.close()
dat = np.loadtxt(file,skiprows=3)
keyword = np.array(header.split())
#NPOINT = len(dat[0:])

# Conversion factors
bar   = 1.E+6                    # 1 bar in dyn/cm2
Tg    = dat[:,0]                 # T [K]
nHtot = dat[:,1]                 # n<H> [cm-3]
lognH = np.log10(nHtot)
press = dat[:,2]/bar             # p [bar]

# Note, in my edited version of GGChem, I output mu on the 4th column!!
mu = dat[:,3]

# Remove pressures that are the same from the grid
press = np.unique(press)

# Here, make a list of molecules and elements you want,
# NOTE: molecules are ALL CAPS e.g. (TIO2), while elements are lowercase after the first letter e.g (He, Ti, Vo)
)
#mol_list = ['H2','He','H2O','CO','CH4','NH3','Na','K']
mol_list = ['H2', 'H', 'He']
nmol = len(mol_list)
x_mol = np.zeros((nmol,NPOINT**2))

# Calcuate the volume mixing ratio
ntot  = 0.0*nHtot
for i in range(4,5+NELEM+NMOLE): # electrons, all atoms, ions and cations
  ntot = ntot + 10.0**dat[:,i]
lntot = np.log10(ntot)
count = 0
for j in range(nmol):
  for i in range(4,5+NELEM+NMOLE):
    if (mol_list[j] == keyword[i]):
      yy = dat[:,i]-lntot            # log10 nmol/ntot
      x_mol[j,:] = yy
      #print(i,keyword[i])
      #ist = NPOINT*150
      #ien = ist+NPOINT
      #print(press[ist],press[ien-1])
      #plt.scatter(Tg[ist:ien],yy[ist:ien])
      #plt.xscale('log')

# Output VMR and mu to a file
output = open('interp_table.txt','w')
output.write(str(NPOINT) + ' ' + str(nmol) + '\n')
for i in range(nmol):
  output.write(mol_list[i] + '\n')
for i in range(NPOINT):
  output.write(str(Tg[i]) + '\n')
for i in range(NPOINT):
  output.write(str(press[i]) + '\n')

e = 0
for j in range(NPOINT):
  print(j)
  for i in range(NPOINT):
    output.write(str(mu[e]) + ' ' + ' '.join(str(g) for g in x_mol[:,e]) + '\n')
    e = e + 1
