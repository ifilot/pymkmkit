CLI reference
=============

The ``pymkmkit`` command-line interface is implemented with Click.

Command overview
----------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Command
     - Purpose
   * - ``freq2yaml``
     - Parse a VASP frequency ``OUTCAR`` and write YAML output.
   * - ``opt2yaml``
     - Parse a VASP optimization ``OUTCAR`` and write YAML output.
   * - ``read_network``
     - Evaluate elementary-step barriers/adsorption heats from a network YAML.
   * - ``evaluate_paths``
     - Evaluate total reaction energy for each named reaction path.
   * - ``build_ped``
     - Build a potential energy diagram for a selected path.

CLI module API
--------------

.. automodule:: pymkmkit.cli
   :members:
   :undoc-members:
   :show-inheritance:
