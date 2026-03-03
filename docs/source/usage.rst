Usage
=====

1. Pareparing the environment
-----------------------------

1.1. Add src/tools/pylib to PYTHONPATH
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   export PYTHONPATH=$PYTHONPATH:/path/to/Gfit/src/tools/pylib

1.2. Add /src/tools/scripts to PATH
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   export PATH=$PATH:/path/to/Gfit/src/tools/scripts

1.3. Compile cpp tools
~~~~~~~~~~~~~~~~~~~~~~
Switch to the src/tools/scripts directory and compile C++ tools:

.. code-block:: bash

   cd /path/to/Gfit/src/tools/scripts
   g++ -std=c++11 -O3  meanforce-alloy-para-all.cpp -o  mf-alloy.x
   cd -

1.4. Give execution permission to fit.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Switch to the src directory and give execution permission to fit.sh:

.. code-block:: bash

   cd /path/to/Gfit/src
   chmod +x fit.sh
   cd -

1.5. Activate python environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Switch to the python environment prepared in the installation step:

See :ref:`prepare-python-environment` for details.


2. Modify bash script
---------------------

2.1 Modify src/fit.sh
~~~~~~~~~~~~~~~~~~~~~

In the file ``fit.sh``, modify the variable ``LMP`` to point to your LAMMPS executable:

.. code-block:: bash

   LMP=/your/lammps/path/lmp

2.2 Modify src/init/rerun/init/gen_property.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In the file ``gen_property.sh``, modify the variable ``LMP`` to point to your LAMMPS executable:

.. code-block:: bash 

   LMP=/your/lammps/path/lmp


3. Run Gfit
---------------------
If you are using batch job system, prepare your job script accordingly like src/fit.pbs. 

.. note::

   Please request at least 64 CPU cores to handle all the tasks.

.. note::
      
   Required (Intel MPI): export I_MPI_HYDRA_BOOTSTRAP=fork
   This enforces fork-based launching, allowing multiple tasks to run concurrently on one node by distributing them across cores.


Then run the fit script:

.. code-block:: bash

   cd /path/to/Gfit/src
   sbatch fit.pbs

else, you can run it directly in terminal:

.. code-block:: bash

   cd /path/to/Gfit/src
   ./fit.sh

4. Output files
----------------

After the fitting is done, you will find potential files of the iteration `i` in the ``/path/to/Gfit/src/pot-ite`` directory,
which is named as ``pot_combine-v{i}.fs``.
. The free energy of each iteration is recorded in the file ``/path/to/Gfit/src/G-v{i+1}/Gmix.dat``.
