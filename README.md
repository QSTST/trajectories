# Molecular structures and MD trajectories of manuscript submitted to PRB

This project lists molecular structures used in moldcular dynmics simulation and output trajectiries files.

```bash
.
├── graphene-leads # Alkane chain placed between graphene electrode
│   ├── graphene-lead-c10.data # Input structure for molecular junctions composed of graphene electrodes
│   ├── clmd # Classical Langevin molecular dynamics
│   │   └── run1 #Simulation times
│   │       └── 300 # Temperature
│   │           ├── trajectories.300.run0.ani # Trajectories file
│   │           └── ...
│   ├── sclmd # Semi-classical Langevin molecular dynamics
│   │   └── run1
│   │       └── 300
│   │           ├── trajectories.300.run0.ani
│   │           └── ...
│   ├── sclmd-harmonic # Semi-classical Langevin molecular dynamics with harmonic force
│   │   └── 300
│   │       ├── trajectories.300.run0.ani
│   │       └── ...
│   └── sclmd-zpe # Semi-classical Langevin molecular dynamics with zero-point energy
│       └── run1
│           └── quantum
│               └── 300
│                   ├── trajectories.300.run0.ani
│                   └── ...
├── gold-leads # Alkane chain placed between gold electrode
│   ├── gold-lead-c10.data # Input structure for molecular junctions composed of gold electrodes
│   ├── clmd
│   │   └── run1
│   │       └── 300
│   │           ├── trajectories.300.run0.ani
│   │           └── ...
│   └── sclmd
│       └── run1
│           └── 300
│               ├── trajectories.300.run0.ani
│               └── ...
└── README.md
```
