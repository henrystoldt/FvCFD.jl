# JuliaCFD

![](https://github.com/henrystoldt/JuliaCFD/workflows/Tests/badge.svg)
# Nomenclature
All units SI  
For the simple variables:

Var | Meaning
--- | ---
e | internal energy  
P | Pressure  
rho | Density  
T | Temperature  
U | Velocity   

And the state & flux variables:  

Var | Meaning | Definition  
--- | --- | ---  
xMom | x-Momentum | rho*U  
eV2 | total energy | rho*(e + U^2/2)  
rhoU2p | flux of x-momentum | rho*U^2 + P   
rhoUeV2PU | x-direction flux of total energy | U*eV2 + P*U  

# References
These are referred to throughout the code:  
**Moukalled et al.**: The Finite Volume Method in Computational Fluid Dynamics: An Advanced Introduction with OpenFOAM and Matlab  
**Hoffman**: Numerical Methods for Engineers and Scientists  
**Versteeg et al.**: An Introduction to Computational Fluid Dynamics: The Finite Volume Method (2nd Edition)  
**Anderson**: Computational Fluid Dynamics: The Basics with Applications  
