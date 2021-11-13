from sclmd.tools import dumpdisp
lammps = "C10.data"
trajectories = ["trajectories.300.run0.ani", "trajectories.300.run1.ani",
                "trajectories.300.run2.ani", "trajectories.300.run3.ani", 
                "trajectories.300.run4.ani", "trajectories.300.run5.ani",
                "trajectories.300.run6.ani", "trajectories.300.run7.ani",
                "trajectories.300.run8.ani", "trajectories.300.run9.ani",
                "trajectories.300.run10.ani", "trajectories.300.run11.ani",
                "trajectories.300.run12.ani", "trajectories.300.run13.ani",
                "trajectories.300.run14.ani", "trajectories.300.run15.ani",
                "trajectories.300.run16.ani", "trajectories.300.run17.ani",
                "trajectories.300.run18.ani", "trajectories.300.run19.ani",
                ]
dumpdisp(lammps, trajectories, index=[1,10,20,30,40,50,60,70,80,90,100,1000], outputname="dispstructure")
