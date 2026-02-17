PyMKMKit documentation
======================

Welcome to the PyMKMKit documentation.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation

Within ab initio microkinetic modeling, chemokinetic networks are constructed
using electronic structure calculations, typically periodic density functional
theory (DFT). To date, there is no standardized format for storing the research
data underlying these models. As a result, authors in the field employ a wide
range of approaches, ranging from large tables of rate constants at fixed
temperatures accompanied by a few renderings of geometric structures, to highly
detailed datasets with interactive scripts that allow input files to be
regenerated under different process conditions.

The purpose of this toolkit is to provide a reproducible and, hopefully, facile
procedure for collecting and storing microkinetic research data. Here, we focus
on the use case most commonly employed in our research group: the collection of
kinetic parameters from plane-wave DFT calculations performed with VASP. In
constructing the microkinetic model, energies and vibrational frequencies are
extracted from the DFT calculations and used, via standard statistical
thermodynamics approaches, to derive the input parameters underlying reaction
rate constants—namely, activation energies and pre-exponential factors for
elementary reaction steps.

All data are stored in YAML format, ensuring cross-platform interoperability
while retaining a high degree of human readability. In addition to storing
information for individual thermodynamic states (namely initial, transition, and
final states) we define a unified YAML schema that encapsulates the complete
microkinetic network by referencing these states. Specifically, each
thermodynamic state is represented by a dedicated YAML file containing the
electronic energy, vibrational frequency data, structural geometry, and
calculation parameters. An overarching YAML file then references these state
files and defines the microkinetic network through the specification of the
elementary reaction steps.