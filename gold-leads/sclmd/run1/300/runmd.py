import time
import numpy as np
from sclmd.baths import ebath
from sclmd.tools import calHF, calTC
from sclmd.lammpsdriver import lammpsdriver
from sclmd.md import md

lammpsinfile = [
    "log    none",
    "units  real",
    "dimension 3 ",
    "boundary   f p p",
    "atom_style full",
    "read_data  avestructure.data",
    "pair_style reax/c NULL",
    "pair_coeff * * AuSCH_2013.ff Au C S H",
    #"compute reax all pair reax/c",
    "fix qep all qeq/reax 1 0.0 10.0 1.0e-6 reax/c",
    "region reg_0   block   0   3.5         INF INF INF INF units   box",
    "region reg_1   block   87  93.0567017  INF INF INF INF units   box",
    "group  g_0 region  reg_0",
    "group  g_1 region  reg_1",
    "fix    g_0 g_0 setforce    0   0   0",
    "fix    g_1 g_1 setforce    0   0   0",
    "min_style  cg",
    "minimize   1e-25   1e-25   5000    10000",
    "unfix  g_0",
    "unfix  g_1",
    "dump   1   all     xyz     1   dump.xyz",
    #"dynamical_matrix all eskm 0.000001 file dynmat.dat",
    "run    0",
]

# temperature
T = 300
delta = 0.1
nstart = 20
nstop = 21
# time = 0.658fs #time unit
dt = 1/0.658
# number of md steps
nmd = 2*10**5
# initialise lammps run
lmp = lammpsdriver(infile=lammpsinfile,eunit="Kcal/mole")

time_start = time.time()

print("initialise md")
fixatoms = [range(0*3, (71+1)*3), range(1016*3, (1090+1)*3)]

# Molecular Junction atom indices
slist = range(180*3, (907+1)*3)
#cutslist = [range(70*3, (89+1)*3),
#            range(90*3, (109+1)*3), range(110*3, (130+1)*3)]
# atom indices that are connecting to debyge bath
ecatsl = range(72*3, (179+1)*3)
ecatsr = range(908*3, (1015+1)*3)
lsection=range(180*3,540*3)
rsection=range(583*3,908*3)
csection=range(540*3,583*3)
#dynmatdat = np.loadtxt("dynmat.dat")# Dynmat units in ps^-2, THz^2
#dynlen = int(3*np.sqrt(len(dynmatdat)/3))
#rpc = 6.582119569e-4
#dynmat = dynmatdat.reshape((dynlen, dynlen))*rpc*rpc# Dynmat units in eV^2

mdrun = md(dt, nmd, T, axyz=lmp.axyz, nstart=nstart,nstop=nstop)
# attache lammps driver to md
mdrun.AddPotential(lmp)
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
#mdrun.noranvel()
mdrun.CalPowerSpec()
atomlist = [ecatsl, slist, ecatsr,lsection,rsection,csection]
mdrun.AddPowerSection(atomlist)
#mdrun.CalAveStruct()
mdrun.SaveTraj(1000)
mdrun.SaveAll()
#mdrun.RemoveNC()
mdrun.Run()

lmp.quit()

calHF()
calTC(delta=delta)
time_end = time.time()
print("time cost", time_end-time_start, "s.")
