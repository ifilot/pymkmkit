PyMKMKit documentation
======================

Purpose
-------

This toolkit standardizes how microkinetic modeling data from plane-wave DFT
(VASP) calculations are collected, stored, and reproduced. It extracts energies
and vibrational frequencies to compute activation energies and pre-exponential
factors for reaction steps. Data are stored in human-readable YAML files, with
separate files for each thermodynamic state and a master file defining the full
reaction network.

Features
--------

- **VASP output parsing**: Convert frequency and geometry optimization
  ``OUTCAR`` files into structured YAML data.
- **Standardized data generation**: Automatically produce consistent, reusable
  YAML representations of thermodynamic states.
- **Reaction network evaluation**: Compute activation barriers and adsorption
  energies from a network definition.
- **Reaction pathway analysis**: Calculate overall reaction energies for defined
  pathways.
- **Potential energy diagrams**: Generate potential energy profiles for selected
  reaction paths.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   tutorial
   api
   cli