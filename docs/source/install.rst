Installation
============

.. _dependencies:

Install dependencies and build source code
------------------------------------------

The source code of IPK depends on libraries you need to install from your distro's repository: boost_ (v1.67+) and zlib_. 
Other depencies are provided as submodules and should be compiled within IPK.

.. _boost: https://www.boost.org/

.. _zlib: https://www.zlib.net/

On a Debian-based system, it can be installed with the following commands:

.. code-block:: console

   $ sudo apt-get update
   $ sudo apt-get install -yq libboost-dev libboost-serialization-dev libboost-filesystem-dev libboost-iostreams-dev libboost-program-options-dev zlib1g-dev
   $ pip3 install click

Clone the repository including the submodules:

.. code-block:: console

    $ git clone --recursive https://github.com/phylo42/IPK.git ipk


Compile the code:

.. code-block:: console

    $ cd ipk && mkdir bin && cd bin
    $ cmake -DHASH_MAP=USE_TSL_ROBIN_MAP -DCMAKE_CXX_FLAGS="-O3" ..
    $ make -j4

Install system-wide:

.. code-block:: console

   $ sudo cmake --install .

Or for the current user (replace ``DIRECTORY`` with any directory you like):

.. code-block:: console

   $ cmake --install . --prefix DIRECTORY
   $ export PATH=DIRECTORY/bin:$PATH


Other dependencies
------------------

To run IPK, you will also need to install either ``raxml-ng`` or ``phyml``. 
Both of them are available in Bioconda_:

.. _Bioconda: https://bioconda.github.io/

.. code-block:: console

    $ conda install raxml-ng

Or:

.. code-block:: console

    $ conda install phyml


.. _test

Check installation
------------------

To check that installation is complete, run the following command to see the help message:

.. code-block:: console

    $ ipk.py --help
    $ ipk.py build --help


