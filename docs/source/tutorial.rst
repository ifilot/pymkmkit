Tutorials
=========

Converting OUTCAR to YAML file
##############################

This tutorial explains how one can parse a single VASP frequency calculation and
constructs a single :code:`.yaml`-based data entry from this calculation.

Prerequisites
-------------

First, we need to prepare an :code:`OUTCAR` file that we wish to parse. Several
example :code:`OUTCAR` are stored in :code:`tests/data` which we can unpack. For
this tutorial, we assume that we are working under :code:`examples/individual`.

Go to this folder

.. code-block:: bash

   cd examples/individual

and run

.. code-block:: bash

   unzip ../../tests/data/OUTCAR_Ru1121_CH.zip

This should give the following output::

   Archive:  ../../tests/data/OUTCAR_Ru1121_CH.zip
      inflating: OUTCAR_Ru1121_CH

Building YAML files
-------------------

The :code:`OUTCAR_Ru1121_CH` file contains a frequency calculation of a CH
molecule adsorbed on a Ru(11-21) surface. To build the :code:`yaml` file, we run

.. code-block:: python

   pymkmkit freq2yaml OUTCAR_Ru1121_CH -o ru1121_c.yaml --average-pairs

The option :code:`--average-pairs` is needed in this case because the
:code:`OUTCAR` file contains a metallic slab on which CH has been adsorbed on
both sides. Thus, when sampling the frequencies for a single CH adsorbate, we
need to average out the frequencies of both the top and bottom CH fragment.

.. note::

   The ``--average-pairs`` option is required for some older slab calculations.
   Historically, adsorbates were placed on **both sides of a metal slab** to
   cancel dipole moments before dipole corrections became standard in VASP.
   Modern calculations typically use dipole corrections instead and only place
   the adsorbate on one side.

   In these legacy setups, two identical adsorbates (top and bottom) produce
   duplicate vibrational modes. Because frequencies are computed using a
   finite-difference approach, small numerical differences appear between the
   paired modes. Averaging the pairs recovers the vibrational frequencies
   corresponding to a single adsorbate while reducing numerical noise.

.. tip::

   To visualize the frequency calculation and inspect vibrational modes, you can
   use `Atom Architect <https://github.com/ifilot/atom-architect>`_. Visualizing
   the displacements can help identify paired modes and confirm the averaging.

Understanding the generated YAML file
-------------------------------------

The generated YAML file stores all information required to reproduce and reuse
the thermodynamic state in a microkinetic model.

Below we highlight the key sections and why they are stored.

``pymkmkit``
~~~~~~~~~~~~

Metadata describing how the file was created.

.. code-block:: yaml

   pymkmkit:
     version: 0.1.0
     generated: 2026-02-26T07:21:02Z

.. list-table::
   :widths: 30 70

   * - Field
     - Purpose
   * - ``version``
     - pymkmkit version used to generate the file
   * - ``generated``
     - Timestamp for provenance and reproducibility

This metadata ensures traceability and helps diagnose differences between datasets.

``structure``
~~~~~~~~~~~~~

Describes the atomic structure used in the calculation.

.. code-block:: yaml

   structure:
     formula: Ru48C2H2
     lattice_vectors: [...]
     coordinates_direct: [...]
     pbc: [true, true, true]

.. list-table::
   :widths: 30 70

   * - Field
     - Purpose
   * - ``formula`` / ``n_atoms``
     - Quick identification of system composition and size
   * - ``lattice_vectors``
     - Defines the periodic simulation cell
   * - ``coordinates_direct``
     - Atomic positions in fractional coordinates
   * - ``pbc``
     - Periodic boundary conditions

Storing the full structure allows the system to be reconstructed and verified.

``calculation``
~~~~~~~~~~~~~~~

Records how the calculation was performed.

.. code-block:: yaml

   calculation:
     code: VASP
     type: frequency
     incar: {...}
     potcar:
       - PAW_PBE Ru
       - PAW_PBE C
       - PAW_PBE H

.. list-table::
   :widths: 30 70

   * - Field
     - Purpose
   * - ``code`` / ``type``
     - Simulation software and calculation type
   * - ``incar``
     - Important settings affecting accuracy and results
   * - ``potcar``
     - Pseudopotentials used

These parameters are essential for reproducibility and validation.

``energy``
~~~~~~~~~~

.. code-block:: yaml

   energy:
     electronic: -433.134185

The electronic energy forms the basis for adsorption energies, reaction energies,
and activation barriers.

``vibrations``
~~~~~~~~~~~~~~

Contains vibrational data required to compute thermodynamic properties.

.. code-block:: yaml

   vibrations:
     frequencies_cm-1:
       - 2886.2
       - 714.1
     partial_hessian:
       dof_labels:
         - 49X
         - 49Y
         - 49Z
         - 50X
         - 50Y
         - 50Z
         - 51X
         - 51Y
         - 51Z
         - 52X
         - 52Y
         - 52Z
       matrix: [...]
     paired_modes_averaged: true

.. list-table::
   :widths: 30 70

   * - Field
     - Purpose
   * - ``frequencies_cm-1``
     - Vibrational modes used in partition functions
   * - ``partial_hessian``
     - Hessian submatrix for adsorbate degrees of freedom (stored for reuse)
   * - ``partial_hessian.dof_labels``
     - Maps each Hessian row/column to an *atom index* and *Cartesian direction* (e.g. ``49X``),
       i.e., which coordinates were displaced in the finite-difference frequency calculation
   * - ``paired_modes_averaged``
     - Indicates symmetric slab modes were averaged

Vibrational data enable computation of:

- zero-point energies
- entropy and heat capacity
- pre-exponential factors in rate constants

The ``dof_labels`` make the partial Hessian interpretable: they encode exactly
which atoms were perturbed and in which Cartesian direction when constructing
the finite-difference Hessian, so the stored matrix can be validated,
visualized, or reprocessed later.

.. tip::

   YAML files generated by **pymkmkit** can be visualized using
   `Atom Architect <https://github.com/ifilot/atom-architect>`_. The program can
   display structures, vibrational modes, and displaced coordinates, making it
   easier to inspect the finite-difference perturbations and verify the results.