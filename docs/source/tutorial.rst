Tutorial: from VASP outputs to a reaction network
=================================================

This tutorial walks through a complete PyMKMKit workflow using the bundled
``examples/Ru1121`` dataset:

1. parse individual ``OUTCAR`` files into YAML state files,
2. define a reaction network,
3. evaluate barriers and pathway energies,
4. generate a potential energy diagram (PED).

Prerequisites
-------------

Install PyMKMKit and verify the CLI is available:

.. code-block:: bash

   pip install pymkmkit
   pymkmkit --help

Example project layout
----------------------

The Ru(11\ :sub:`2`\ 1) example uses separate folders for stable states,
transition states, gas-phase species, and one network definition file.

.. code-block:: text

   examples/Ru1121/
   ├── GAS/
   ├── ISFS/
   ├── TS/
   └── network.yaml

Step 1: create YAML files from OUTCAR files
-------------------------------------------

Use ``opt2yaml`` for optimized geometries and ``freq2yaml`` for frequency jobs.

.. code-block:: bash

   pymkmkit opt2yaml path/to/OUTCAR -o ISFS/co.yaml
   pymkmkit freq2yaml path/to/OUTCAR -o TS/co_diss.yaml

If your system contains paired adsorbates and should be averaged in sequential
mode pairs, add ``--average-pairs``:

.. code-block:: bash

   pymkmkit freq2yaml path/to/OUTCAR -o TS/ch2_hydr.yaml --average-pairs

The repository includes a convenience script showing this harvesting process for
all states in a dataset.

.. literalinclude:: ../../examples/Ru1121/harvest.sh
   :language: bash
   :caption: Example batch harvesting script used to generate state YAML files.

Step 2: inspect a generated state file
--------------------------------------

A state YAML captures structure, calculation metadata, and energies. Frequency
states also include vibrational modes and (when present) imaginary modes and
partial Hessian data.

.. literalinclude:: ../../examples/Ru1121/TS/co_diss.yaml
   :language: yaml
   :caption: Transition-state YAML example.

For optimized minima, the same schema is used without the ``vibrations`` block:

.. literalinclude:: ../../examples/Ru1121/ISFS/empty.yaml
   :language: yaml
   :caption: Stable-state YAML example.

Step 3: define and understand the network file
----------------------------------------------

The network file links named states to elementary steps and optional paths.

.. literalinclude:: ../../examples/Ru1121/network.yaml
   :language: yaml
   :caption: Full network definition for the Ru(11\ :sub:`2`\ 1) methanation example.

Key sections in ``network.yaml``:

- ``stable_states`` and ``transition_states`` map symbolic names to state files.
- ``network`` defines each elementary step:

  - ``type: surf`` uses ``forward``/``backward`` blocks with ``ts`` and ``is`` terms.
  - ``type: ads`` uses ``is`` and ``fs`` terms for adsorption heat.

- ``normalization`` rescales per-site/per-event energetics.
- ``paths`` defines cumulative pathways as weighted sums of step reaction heats.

Step 4: evaluate elementary-step energetics
-------------------------------------------

Print forward/reverse barriers (or adsorption heat) for every step:

.. code-block:: bash

   pymkmkit read_network examples/Ru1121/network.yaml

This command resolves each referenced state file, computes electronic and ZPE
contributions, and reports total values in eV.

Step 5: evaluate pathway energies
---------------------------------

Compute net pathway energies from the ``paths`` section:

.. code-block:: bash

   pymkmkit evaluate_paths examples/Ru1121/network.yaml

Each path accumulates ``factor * reaction_heat_total`` for all listed steps.

Step 6: generate a potential energy diagram
-------------------------------------------

Build a PED for a selected path and save it as an image:

.. code-block:: bash

   pymkmkit build_ped examples/Ru1121/network.yaml methanation ped_methanation.png

If no output filename is provided, the figure is displayed interactively.

What PyMKMKit can do
--------------------

Using this workflow, PyMKMKit supports:

- extracting structured, reusable YAML state data from VASP outputs,
- tracking vibrational information and ZPE corrections,
- evaluating forward/reverse barriers and adsorption heats from a unified
  network schema,
- aggregating elementary steps into named mechanistic paths,
- producing publication-ready potential-energy diagrams.
