# idugksFoam

An OpenFOAM inplementation for  incompressible isothermal fluid. While modified from [dugksFoam][1] \[1\], idugksFoam is much more efficient  when applied for incompressible flow such as porous media flow.

## key feature

- A standard OpenFOAM solver using discrete unified gas kinetic scheme \[2\].

- Based on BGK collision model, which is suitable for incompressible isothermal fluid

- Bounce-back scheme is used for velocity boundary condition, which is much more effiecient when applied for porous media flow

- Arbitrary unstructured meshes \[3\].

- Easyly applied to REV scale for engineer computing

## Reference

- \[1\] L. Zhu, S. Chen, Z. Guo, dugksFoam: An open source OpenFOAM solver for the Boltzmann model equation, [Comp. Phys. Commun., 213(2017) 155-164][2]

- \[2\] Z. Guo, K. Xu, R. Wang, Discrete unified gas kinetic scheme for all Knudsen number flows: low-speed isothermal case, [Phys. Rev. E, 88 (2013) 033305][3]

- \[3\] L. Zhu, Z. Guo, K. Xu, Discrete unified gas kinetic scheme on unstructured meshes, [Comp. Fluids, 127(2016) 211-225][4]

[1]: https://github.com/zhulianhua/dugksFoam "dugksFoam"

[2]: https://www.sciencedirect.com/science/article/pii/S0010465516303642 "dugksFoam article"

[3]: https://journals.aps.org/pre/abstract/10.1103/PhysRevE.88.033305 "DUGKS for isothermal fluid"

[4]: https://www.sciencedirect.com/science/article/pii/S0045793016000177 "DUGKS on unstructed mesh"
