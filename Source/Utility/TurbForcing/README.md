# Turbulent Forcing

Adds support for the capability for maintained homogeneous isotropic turbulence. A source term in the momentum equation injects energy at large scales, which naturally cascades to form homogeneous isotropic turbulence.

## Acknowledgments

Original implementation: Andrew Aspden, Nikos Nikiforakis, Stuart Dalziel, John Bell (see citation below).

Integration into IAMR: Candace Gilet.

Integration into PelePhysics: Andrew Aspden, Edward Hunt, Thomas Howarth.

## Citation
To cite the forcing scheme, please cite the [CAMCOS article](http://dx.doi.org/10.2140/camcos.2008.3.103)
```
@article{aspden2008analysis,
  title={{Analysis of implicit LES methods}},
  author={Aspden, Andrew and Nikiforakis, Nikos and Dalziel, Stuart and Bell, John},
  journal={Communications in Applied Mathematics and Computational Science},
  volume={3},
  number={1},
  pages={103--126},
  year={2008},
  publisher={Mathematical Sciences Publishers}
}
```