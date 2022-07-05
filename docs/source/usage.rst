User Guide
============


.. _oneliner:

Give me a one-liner
--------------------

.. warning::

   This section gives a lazy copy-paste for experienced users who understand the meaning of all parameters. 
   If you are not one of them, you may want to :ref:`read the docs<readthedocs>` instead.


.. note::
   We assume that you want to use RAxML-ng as the external ancestral reconstruction tool. If you want to use PhyML instead, change the command accordingly.



.. tabs::

   .. tab:: DNA

      .. code-block:: console

         $ python xpas.py build -w workdir \
                                -r alignment.fasta \
                                -t tree.newick \
                                -m GTR \
                                -a 1.0 \
                                -k 10


   .. tab:: Proteins

      .. code-block:: console

         Temporarily disabled in this version of XPAS.

.. _readthedocs:


I want to *read the docs*
--------------------------

To create a phylo-k-mer database, you need to provide the following input:

- Multiple alignment of reference sequences
- Reference phylogenetic tree, tips of which are in 1-to-1 correspondance with the sequences of the alignment. The tree should be inferred with the maximum likelihood method
- Evolutionary model used to infer the tree


.. _parameters:


Besides that, you need to define a number of important parameters discussed below.

- Getting started
- Parameter reference


