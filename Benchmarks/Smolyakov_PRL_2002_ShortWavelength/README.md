Benchmark using local, sheareless slab geometry
===============================================

  Note : Needs cleanup


This geometry is a shearless, slab geometry. The parallel wavenumber accounting
for the Landau damping is set directly by $k_\parallel = k_z$. As we 
basically fix kx,ky and kz, only velocity space (v,m) and kinetic
species are evolved. 

We reproduce Smolyakov et al. results, who found the short-wavelength ITG mode,
a further destabilization of the ITG mode beyond ky > 1.


Benchmark
____________

Using this test, we benchmark the gyro-averaging modules, as well as the inclusion of
kinetic species.

Additionally, we use the eigenvalue solver to find that the short-wavelength branch is
indeed the ITG-sw branch and does not account from the ETG.


Results
_________


## Benchmark results (with adiabatic electrons)
![Benchmark results with adiabatic electrons](Smolyakov_Adiabatic.png)

## Benchmark results (with kinetic electrons)
![Caption text](Smolyakov_Kinetic.png "")

## Benchmark results (Eigenvalues solver resolved ETG and ITG mode)
![Caption text](Smolyakov_Eigenvalue.png "")
