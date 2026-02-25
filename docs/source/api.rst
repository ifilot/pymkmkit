API reference
=============

This section provides an overview of all core classes and functions and then
expands each module with full per-member documentation.

Classes overview
----------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Class
     - Purpose
   * - :class:`pymkmkit.network_reader.State`
     - Immutable thermodynamic-state container used for energy and ZPE lookup.
   * - :class:`pymkmkit.network_reader.ElementaryStep`
     - Immutable container for evaluated forward/reverse barriers and reaction heats.
   * - :class:`pymkmkit.network_reader.ReactionPath`
     - Immutable container for net reaction energy of a named pathway.
   * - :class:`pymkmkit.yaml_writer.InlineList`
     - List subclass that forces flow-style YAML rendering.

Functions overview
------------------

.. autosummary::

   pymkmkit._version.get_version
   pymkmkit.vasp_freq.extract_incar_settings
   pymkmkit.vasp_freq.extract_potcar_info
   pymkmkit.vasp_freq.extract_ionic_energies
   pymkmkit.vasp_freq.extract_total_energies
   pymkmkit.vasp_freq.extract_frequencies
   pymkmkit.vasp_freq.extract_perturbed_hessian
   pymkmkit.vasp_freq.frequencies_from_partial_hessian
   pymkmkit.vasp_freq.formula_from_atom_order
   pymkmkit.vasp_freq.lattice_vectors
   pymkmkit.vasp_freq.geometry_direct_strings
   pymkmkit.vasp_freq.parse_vasp_frequency
   pymkmkit.vasp_freq.parse_vasp_optimization
   pymkmkit.vasp_freq.average_mode_pairs
   pymkmkit.yaml_writer.represent_inline_list
   pymkmkit.yaml_writer.clean_none
   pymkmkit.yaml_writer.write_yaml
   pymkmkit.network_reader.read_network
   pymkmkit.network_reader.evaluate_paths
   pymkmkit.network_reader.build_ped

Module details
--------------

pymkmkit.vasp_freq
~~~~~~~~~~~~~~~~~~

.. automodule:: pymkmkit.vasp_freq
   :members:
   :undoc-members:
   :show-inheritance:

pymkmkit.network_reader
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: pymkmkit.network_reader
   :members:
   :undoc-members:
   :show-inheritance:

pymkmkit.yaml_writer
~~~~~~~~~~~~~~~~~~~~

.. automodule:: pymkmkit.yaml_writer
   :members:
   :undoc-members:
   :show-inheritance:

pymkmkit._version
~~~~~~~~~~~~~~~~~

.. automodule:: pymkmkit._version
   :members:
   :undoc-members:
   :show-inheritance:
