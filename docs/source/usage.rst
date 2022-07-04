Usage
======


.. _oneliner:

Give me a one-liner
--------------------

Use the python wrapper ``xpas.py`` to create a new database:


.. tabs::

   .. tab:: DNA

    .. code-block:: console

        $ python xpas.py build -s nucl \
                               -b `which raxml-ng` \
                               -w workdir \
                               -r alignment.fasta \
                               -t tree.newick \
                               -k 10 \


   .. tab:: Proteins

    .. code-block:: console

        $ python xpas.py build -s amino \
                               -b `which raxml-ng` \
                               -w workdir \
                               -r alignment.fasta \
                               -t tree.newick \
                               -k 4 \

.. note::
    We assume that you want to use RAxML-ng as the external ancestral reconstruction tool. If you want to use PhyML instead, change the command accordingly.


.. _readthedocs:

I want to *read the docs*
--------------------------




