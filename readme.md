# Thermodynamic library: humidityRhoThermo
This library takes humidity effects into account and is mainly built for HVAC analysis. The validity is limited to around 1 bar and -50 degC to 100 degC. The library is built for different OpenFOAM versions. For more complex analysis (different pressure than atmospheric one or higher temperatures) one should use the reactingFoam for taking into account the humidity.

# How to use it

Clone the repository to any place you want using the following command:
```console
@-: git clone https://github.com/shor-ty/humidityRhoThermo.git
```

After that load your OpenFOAM environment (if not already happend) and move into the repository. Here checkout your version you want:
```console
@~: git checkout <TAB><TAB>
OpenFOAM-6.x
OpenFOAM-7.x
OpenFOAM-v8
OpenFOAM-v9
OpenFOAM-v1712
OpenFOAM-v1806
OpenFOAM-v1812
OpenFOAM-v1906
OpenFOAM-v1912
OpenFOAM-v2012
OpenFOAM-v2112
@~: git checkout OpenFOAM-v9
```
After you switched to your OpenFOAM version, you have to compile the library first and after that the solver:
```console
@~: cd src/thermodynamic/basic/
@~: wmake libso
@~: cd ../../../applications/solvers/heatTransfer/buoyantHumidityPimpleFoam/
@~: wmake
```
You are done. After that you can use the buoyantHumidityPimpleFoam solver.

# Sponsored
This project was sponsered by Tian Building Engineering

# Rebuilt
This version was initially created by Tobias Holzmann for OpenFOAM-v9 but needed to be changed and modified as the Foundation version was changed in addition. In the updated implementation, humidity values are read from a time-variant table, allowing for the specification of humidity levels that change over time. This approach employs the Function1<scalar> class, enabling the simulation to update the boundary condition at each time step based on predefined values in an external file or table. By incorporating time-dependent boundary conditions, the simulation can now more accurately reflect varying environmental conditions, enhancing the realism and applicability of the results.
The rebuilt for the latest OpenFOAM-v9 version was done by Dr. Robert Castilla. Thanks for that.
