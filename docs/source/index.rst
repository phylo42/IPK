IPK: Inference of Phylo-K-mers
===================================

IPK_ is a program for inference of phylo-k-mers. You will need to use :program:`IPK`
if you want to use any of the following:

1. Alignment-free phylogenetic placement with EPIK_.
2. Detection of novel virus recombinants with SHERPAS_.
3. Protein gene family classification with CLAPPAS_.
4. Or do your own research experimenting with phylo-k-mers.

.. _IPK: https://github.com/phylo42/IPK

Check out the :doc:`install` and :doc:`usage` for general information. 
See :doc:`cli` for detailed information about program options.


.. warning::
   Although I2L (IPK's Interface Library) is currently stable, it is still under development. 
   Future updates of the library may not be backward compatible. 
   This means that databases built with the current version of IPK
   may not be compatible with the future versions of phylo-k-mer-based tools 
   (e.g., :program:`EPIK` and :program:`SHERPAS`). Check out the changelog_ for breaking changes.

   .. _changelog: https://github.com/phylo42/I2L/blob/master/CHANGELOG.txt 
   .. _EPIK: https://github.com/phylo42/EPIK
   .. _SHERPAS: https://github.com/phylo42/sherpas
   .. _CLAPPAS: http://blinard.net/


Contents
--------

.. toctree::

   install
   usage
   cli
