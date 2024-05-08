# KerrSpinningSecondary

A Mathematica package for the trajectory and related quantities for a spinning body orbiting a Kerr black hole.

**This package is in heavy development and not ready for general use**

Initial probably easiest to include linear-in-spin results for spin-aligned (anti-aligned) orbits. That said, when designing the API we want to allow it to describe:

- Non-aligned spins, e.g., `KerrSpinEnergy[a,p,e,x,{σ1,σ2,σ3}]`
- Different spin supplementary conditions, e.g., `KerrGeoEnergy[a,p,e,x,s, "SpinSupplementaryCondition"->"linear/Tulczyjew/Pirani/OKS/etc"]`

As much as is possible let's try to follow the API and code structure of the [KerrGeodesics](https://bhptoolkit.org/KerrGeodesics/) package.



