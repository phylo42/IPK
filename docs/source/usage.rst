Usage
=====

.. _installation:

Installation
------------

XPAS depends on a few libraries you need to install from your distro repository: `boost` (v1.67+) and `zlib`. 
Other depencies are provided as submodules and should be compiled within XPAS.

On a Debian-based system, it can be installed with the following command:

.. code-block:: console

   $ sudo apt-get update
   $ sudo apt-get install -yq libboost-dev libboost-serialization-dev libboost-filesystem-dev libboost-iostreams-dev libboost-program-options-dev zlib1g-dev

When you clone the repository, do not forget to clone the submodules:

.. code-block:: console

    git clone --recursive https://github.com/phylo42/xpas.git


To compile the library:

.. code-block:: console

    cd xpas && mkdir bin && cd bin
    cmake -DHASH_MAP=USE_TSL_ROBIN_MAP -DCMAKE_CXX_FLAGS="-O3" ..
    make -j4

.. _dependencies:

Other dependencies
------------

To compute phylo-*k*-mers, you need to install `raxml-ng` or `phyml`.

.. _computing:

Computing a database
------------

Use the python wrapper `xpas.py` to create a new database:

.. code-block:: console

    python xpas.py build -s [nucl|amino] -b `which raxml-ng` -w workdir -r alignment.fasta -t tree.newick -k 10
