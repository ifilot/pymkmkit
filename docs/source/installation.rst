Installation
============

We recommend installing **pymkmkit** inside a Python virtual environment.

Why use a virtual environment?
------------------------------

A virtual environment creates an isolated Python workspace for your project.
This prevents dependency conflicts between different projects and avoids
modifying your system-wide Python installation.

Using virtual environments helps you:

- keep project dependencies separate
- ensure reproducibility
- avoid permission issues
- safely test different package versions

Create and activate a virtual environment
-----------------------------------------

First, create a virtual environment:

.. code-block:: bash

   python -m venv .venv

Activate it:

**Linux / macOS**

.. code-block:: bash

   source .venv/bin/activate

**Windows (PowerShell)**

.. code-block:: powershell

   .venv\Scripts\Activate.ps1

Once activated, your terminal prompt should show ``(.venv)``.

Install pymkmkit
----------------

With the virtual environment active, install the package:

.. code-block:: bash

   pip install pymkmkit

You can verify the installation with:

.. code-block:: bash

   pip show pymkmkit