XPAS: computation of phylo-k-mers
===================================

XPAS_ (a recursive acronym for XPAS: Phylo-k-mers of Ancestral Sequences) is a
tool for creating phylo-k-mer databases. You may want to use :program:`XPAS`
if you do one of the following:

1. Alignment-free phylogenetic placement with RAPPAS2_.
2. Detection of novel virus recombinants with SHERPAS_.
3. Protein gene family classification with CLAPPAS_.
4. Your own research that experiments with phylo-k-mers.

.. _XPAS: https://github.com/phylo42/xpas

Check out the :doc:`install` and :doc:`usage` for general information. 
See :doc:`cli` for detailed information about program options.


.. warning::
   Although XCL (XPAS Core Library) is currently stable, it is still under development. 
   Therefore, future updates of the library may not be backward compatible. 
   Unfortunately, this means that databases built with the current version of XPAS 
   may not be compatible with the future versions of phylo-k-mer-based tools 
   (e.g., :program:`RAPPAS2` and :program:`SHERPAS`). Check out the changelog_ for breaking changes.

   .. _changelog: https://github.com/phylo42/xpas/blob/master/CHANGELOG.txt 
   .. _RAPPAS2: https://github.com/phylo42/rappas2
   .. _SHERPAS: https://github.com/phylo42/sherpas
   .. _CLAPPAS: http://blinard.net/


Contents
--------

.. toctree::

   install
   usage
   cli
