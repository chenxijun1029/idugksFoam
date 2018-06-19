# idugksFoam

An OpenFOAM inplementation for  incompressible isothermal fluid. While based on [dugksFoam][1] \[1\], idugksFoam is much more efficient when applied for incompressible flow such as porous media flow.

## Key features

- A standard OpenFOAM solver for Boltzmann model equation using discrete unified gas kinetic scheme \[2\]

- Maintain major feautres of dugksFoam including 1D & 2D & 3D in a single solver, various boundary condition types and arbitrary unstructured meshes \[3\]

- Based on BGK collision model, which is suitable for incompressible isothermal fluid

- Bounce-back scheme is used for velocity boundary condition, which is much more effiecient when applied for porous media flow

- Easyly applied to REV scale for engineer computing

## Installation

OpenFOAM-2.4.0 is supported only at present, with Intel/gcc compilers.

    of240 # load the environment of OpenFOAM-2.4.0
    cd idugksFoam/src
    ./Allwmake

Done!

## Usage

There is a typical example in demo/testPressureBoundary for simulation of poiseuille flow, from which you could know how to set boundary conditon and parameters. A full documentation will be written if necessary. As for other module such as meshes, control and postprocession, just find answer in [OpenFOAM User Guide][5]. 

If there is any question, leave it here.

## Reference

- \[1\] L. Zhu, S. Chen, Z. Guo, dugksFoam: An open source OpenFOAM solver for the Boltzmann model equation, [Comp. Phys. Commun., 213(2017) 155-164][2]

- \[2\] Z. Guo, K. Xu, R. Wang, Discrete unified gas kinetic scheme for all Knudsen number flows: low-speed isothermal case, [Phys. Rev. E, 88 (2013) 033305][3]

- \[3\] L. Zhu, Z. Guo, K. Xu, Discrete unified gas kinetic scheme on unstructured meshes, [Comp. Fluids, 127(2016) 211-225][4]

[1]: https://github.com/zhulianhua/dugksFoam "dugksFoam"

[2]: https://www.sciencedirect.com/science/article/pii/S0010465516303642 "dugksFoam article"

[3]: https://journals.aps.org/pre/abstract/10.1103/PhysRevE.88.033305 "DUGKS for isothermal fluid"

[4]: https://www.sciencedirect.com/science/article/pii/S0045793016000177 "DUGKS on unstructed mesh"

[5]: https://www.openfoam.com/documentation/user-guide/ "OpenFOAM User Guide"
