Installation 
============

1. Install Intel oneAPI
-----------------------
Download and install Intel oneAPI Base Toolkit and HPC Toolkit from `Intel oneAPI <https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html>`_.


2. LAMMPS Requirements
----------------------
We test the code with the LAMMPS version : LAMMPS/2024.8.29. 
We use the above Intel MPI to build LAMMPS with following packages enabled:

`EXTRA-PAIR <https://docs.lammps.org/Packages_details.html#pkg-extra-pair>`_,
`EXTRA-FIX <https://docs.lammps.org/Packages_details.html#pkg-extra-fix>`_,
`MANYBODY <https://docs.lammps.org/Packages_details.html#pkg-manybody>`_,

.. _prepare-python-environment:

3. Prepare python environment
-----------------------------
Install following python packages

.. code-block:: bash
   
   conda install numpy
   conda install scipy
   conda install matplotlib
   conda install scikit-learn
   conda install jq