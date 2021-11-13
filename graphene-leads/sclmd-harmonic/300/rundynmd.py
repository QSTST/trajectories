import time
import numpy as np
from sclmd.baths import ebath
from sclmd.tools import calHF, calTC
from sclmd.lammpsdriver import lammpsdriver
from sclmd.md import md
lammpsinfile = [
    #"log none",
    "units metal ",
    "dimension 3 ",
    "boundary f p p",
    "atom_style full",
    "read_data C10.data",
    "pair_style rebo ",
    "pair_coeff * * CH.rebo C H",
    "region reg_0   block   0.0     4.261       INF INF INF INF units   box",
    "region reg_1   block   59.2725 63.5335007  INF INF INF INF units   box",
    "group  g_0 region  reg_0",
    "group  g_1 region  reg_1",
    "fix    g_0 g_0 setforce    0   0   0",
    "fix    g_1 g_1 setforce    0   0   0",
    "min_style  cg",
    "minimize   1e-25   1e-25   5000    10000",
    "unfix  g_0",
    "unfix  g_1",
    "dump   1   all     xyz     1   dump.xyz",
    "dynamical_matrix all eskm 0.000001 file dynmat.dat",
    "run    0",
]
# temperature
T = 300
delta = 0.1
nstart = 0
nstop = 20
# time = 0.658fs #time unit
dt = 0.5/0.658
# number of md steps
nmd = 2*10**5
# initialise lammps run
lmp = lammpsdriver(infile=lammpsinfile)
time_start = time.time()
print("initialise md")
fixatoms = [range(0*3, (23+1)*3), range(304*3, (327+1)*3)]
# Molecular Junction atom indices
slist = list(range(149*3, (178+1)*3))
# atom indices that are connecting to debyge bath
ecatsl = list(range(24*3, (59+1)*3))
ecatsr = list(range(268*3, (303+1)*3))
# if slist is not given, md will initialize it using xyz
#mdrun = md(dt, nmd, T, syslist=None, axyz=lmp.axyz, writepq=False, rmnc=False,
#           nstart=nstart, nstop=nstop, npie=1, constr=fixatoms, nstep=1000)
## attache lammps driver to md
#mdrun.AddPotential(lmp)

dynmatdat = np.loadtxt("dynmat.dat")# Dynmat units in ps^-2, THz^2
dynlen = int(3*np.sqrt(len(dynmatdat)/3))
rpc = 6.582119569e-4
dynmat = dynmatdat.reshape((dynlen, dynlen))*rpc*rpc# Dynmat units in eV^2

mdrun = md(dt, nmd, T, axyz=lmp.axyz, dyn=dynmat,nstart=nstart,nstop=nstop)

# unit in 0.658211814201041 fs
damp = 100/0.658211814201041
etal = (1.0/damp)*np.identity(len(ecatsl), np.float)
etar = (1.0/damp)*np.identity(len(ecatsr), np.float)
# atom indices that are connecting to bath
ebl = ebath(ecatsl, T*(1+delta/2), mdrun.dt, mdrun.nmd,
            wmax=1., nw=500, bias=0.0, efric=etal, classical=False, zpmotion=False)
mdrun.AddBath(ebl)
ebr = ebath(ecatsr, T*(1-delta/2), mdrun.dt, mdrun.nmd,
            wmax=1., nw=500, bias=0.0, efric=etar, classical=False, zpmotion=False)
mdrun.AddBath(ebr)

mdrun.AddConstr(fixatoms)
mdrun.noranvel()
mdrun.CalPowerSpec()
atomlist = [ecatsl, slist, ecatsr]
mdrun.AddPowerSection(atomlist)
#mdrun.CalAveStruct()
mdrun.SaveTraj(1000)
mdrun.RemoveNC()
mdrun.Run()
lmp.quit()
calHF()
calTC(delta=delta, dlist=1)
time_end = time.time()
print('time cost', time_end-time_start, 's.')
