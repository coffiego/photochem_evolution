# photochem_evolution

## Description

This photochemistry model is originally developed by Chaffin+2017.
Original code can be found from the link below:
https://github.com/planetarymike/chaffin_natgeo_mars_photochemistry

This code is modified by Koyama to investigate atmospheric evolution of Mars.
This code solves 1-D photochemistry and vertical tranport with time-dependent EUV flux and volcanic outgassing.

## Requirements
- Julia 1.5.3
- Python3 (only for plot)

Packages
- PyPlot
- HDF5
- JLD
- DelimitedFiles
- SparseArrays
- LinearAlgebra
- specialfunctions

## Installation

Cloning in any directory with the command below
```
$ git clone git@github.com:coffiego/photochem_evolution.git
```

## Usage
Please change the directory path to yours in the code first.

```
folder_directory = "your/directory/"
```

a) Change paramters in photochemistry_c.jl

b) Prepare initial profile
Please create the initial density profile, and write the path into variable readfile in photochemistry_c.jl. 

```
readfile = "path/to/your/inputdata.h5"
```

c) run run_c.jl

```
$ julia run_c.jl
```


#### Boundary Condtition

Upper BD
- C, O, H, H2:
modified jeans escape -> effusion velocity is assigned
- C,O:
non-thermal escape for each EUV flux

Lower BD
- CO2
time-dependent volcanic flux
- Deposition velocity:
H2O2,O3
You can modify dictionary object of "speciesbclist" in the code.