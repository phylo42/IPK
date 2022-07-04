Usage
=====


.. _computing:

Computing a database
------------

Use the python wrapper ``xpas.py`` to create a new database:

.. code-block:: console

    python xpas.py build -s [nucl|amino] -b `which raxml-ng` -w workdir -r alignment.fasta -t tree.newick -k 10
