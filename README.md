################### DLCA ###############################

This repo contains code that simulates particulate aggregation, building on the Diffusion-Limited Cluster Aggregation (DLCA) model by adding rotational Brownian motion, realistic size-depedent diffusivities and settling under gravity.
System is hard-coded to be three-dimensional, and periodic boundary conditions (PBC) are implemented in all spatial directions.

### Basic Structure ###
    'main100.cpp' is the driver
    'funcs100.cpp' contains the relevant functions to perform DLCA, following PBC routine
    'common100.h' contains global parameters and functions declarations
    
### TO COMPILE ###
    Type './compile100.sh'

### TO RUN ###
    Type './run100.sh'
    
Note: make sure to have made the shell scripts into executables 'chmod +x <filename>'
