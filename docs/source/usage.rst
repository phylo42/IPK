***************
Getting started
***************


.. _oneliner:

I want a quick example
===============

.. warning::

   This section gives a copy-paste for experienced users who understand the meaning of all parameters. 
   If you are not one of them, you may want to :ref:`read the docs<readthedocs>` instead.


.. note::
   We assume that you want to use RAxML-ng as the external ancestral reconstruction tool. If you want to use PhyML instead, change the command accordingly.



.. tabs::

   .. tab:: DNA

      .. code-block:: console

         $ python ipk.py build -w workdir \
                               -r alignment.fasta \
                               -t tree.newick \
                               -m GTR \
                               -a 0.42 \
                               -k 10


   .. tab:: Proteins

         Temporarily disabled in this version of IPK.


.. _readthedocs:

I want to *read the docs*
=========================


First, you must ensure that your data satisfies the assumptions that :program:`IPK` makes. 
Then, you need to ensure that specific default parameters fit your case. 

Here is a checklist for any phylo-k-mer-based analysis. This section will talk about all
of these points in detail below.

- :option:`-r`: Reference sequences are short and aligned in a .fasta file
- :option:`-t`: Reference tree is inferred with ML, rooted, and in accordance with the alignment
- :option:`-m` or :option:`--ar-parameters`: Evolutionary model set matches the one used to infer the tree 
- :option:`-k` and :option:`-o`: The values fit the data type and size
- (optional) :option:`-u`: desired k-mer filtering rate is set.

.. danger::

   Failing to meet any of these requirements (but the last one) will lead to
   meaningless results or incorrect operation.


Reference multiple alignment
----------------------------

Phylo-k-mers are applied for phylogenies inferred from either short genetic markers 
common in metabarcoding (such as 16S and others) or short viral genomes. 
Therefore, we do not know how phylo-k-mers in general and :program:`IPK` in particular 
will work with data representing much longer sequences, e.g., hundreds of kilobases.

:option:`-r` or :option:`--refalign` define the path to the reference sequences. Your reference sequences should be aligned 
and stored in ``fasta`` format. 

.. note::

   Depending on the gap percentage of your data, you may want to change the default value 
   of :option:`--reduction-ratio` (see command line reference for detail).


Reference phylogenetic tree
---------------------------


The reference phylogenetic tree provided with :option:`-t` or :option:`--reftree` should satisfy 
the following criteria:

1. The tips of the tree are in one-to-one correspondance with the reference sequences defined
by :option:`-r`.

2. It was inferred with a maximum likelihood method. You can use `PhyML`_, `RAxML-ng`_, 
or `IQ-TREE`_ for these purposes. 


.. tip::

   If the tree is too large for ML inference, you can infer it with any other method. 
   Then, keep its topology but reoptimize the branch lengths with maximum likelihood.

3. The tree should be rooted.



.. _PhyML: http://www.atgc-montpellier.fr/phyml/
.. _RAxML-ng: https://github.com/amkozlov/raxml-ng
.. _IQ-TREE: http://www.iqtree.org/




Ancestral reconstruction
------------------------

Ancestral reconstruction is a required step of phylo-k-mer computation. 
To complete it, :program:`IPK` runs either :program:`PhyML` or :program:`RAxML-ng`. 
By default, it searches for ``raxml-ng`` in your PATH. Should it be absent 
(or you want to use a specific version of the software), 
you must provide its path via :option:`-b` or :option:`--ar`.

Besides that, you must provide the evolutionary model used to infer the tree, as long as
necessary parameters for that model. There are two ways of doing it.


Passing parameters via command line
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pass the evolutionary model via :option:`-m` or :option:`--model`. 
We follow the :ref:`RAxML-ng's naming notation<models>`, i.e., it must be one of the 
following ones:

+------------+-----------------------------------------------------------------+
| Data type  | Model name                                                      |
+------------+-----------------------------------------------------------------+
| DNA        | JC, K80, F81, HKY, TN93ef, TN93, K81, K81uf,                    |
|            |                                                                 |
|            | TPM2, TPM2uf, TPM3, TPM3uf, TIM1, TIM1uf, TIM2,                 |
|            |                                                                 |
|            | TIM2uf, TIM3, TIM3uf, TVMef, TVM, SYM, GTR                      |
+------------+-----------------------------------------------------------------+
| Proteins   | Blosum62, cpREV, Dayhoff, DCMut, DEN, FLU, HIVb,                |
|            |                                                                 |
|            | HIVw, JTT, JTT-DCMut, LG, mtART, mtMAM, mtREV,                  |
|            |                                                                 |
|            | mtZOA, PMB, rtREV, stmtREV, VT, WAG, LG4M,                      |
|            |                                                                 |
|            | LG4X, PROTGTR                                                   |
+------------+-----------------------------------------------------------------+


.. _models: https://github.com/amkozlov/raxml-ng/wiki/Input-data#evolutionary-model


Command line supports among-site rate heterogeneity parameters for the GAMMA model: 
:option:`-a` (:option:`--alpha`) and :option:`-c` (:option:`--categories`). Those should be
set accordingly to their values. In case another model was used, see the next section
to pass the parameters via config file.


Passing parameters via config file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Alternatively, you can create a .json-formatted config file to pass arbitrary arguments to 
the ancestral reconstruction tool. For example, let us say the tree was inferred 
with :program:`RAxML-ng` under the GTR model as follows:

.. code-block:: console

   $ raxml-ng --all --msa msa.fasta --model GTR+G4+FC --tree pars{20} --bs-trees 200

resulting in the best model of 

.. code-block:: console

   $ cat msa.fasta.raxml.bestModel

   ``GTR{1.229596/5.824854/2.865999/2.632363/9.525945/1.000000}+FC+G4m{0.492429}, noname = 1-10262``

Then, the following config file should be created and passed via :option:`--ar-config` to
:program:`IPK`:


.. code-block:: json
   :caption: An example of a .json-formatted config for :program:`RAxML-ng`

   {
      "arguments": {
         "data-type": "DNA",
         "model": "GTR{1.229596/5.824854/2.865999/2.632363/9.525945/1.000000}+FC+G4m{0.492429}",
         "opt-branches": "on",
         "opt-model": "on",
         "blopt": "nr_safe"
      }
   }


.. note::
   The evolutionary model should be reoptimized after the extention of the tree done by :program:`IPK`. This
   is why we set ``opt-model: on`` in this example.



K-mer size and score threshold
------------------------------

K-mer size and score threshold parameter (:option:`-k` and :option:`--omega`)
have **tremendous** performance impact on :program:`IPK` and tools using 
resulting phylo-k-mers, both in running time and memory. It is important to 
understand what those parameter mean.

K-mer size
~~~~~~~~~~
K-mer size determines how long considered k-mers are, or, in other words, 
**how many hypothetical strings we consider**. The number of considered k-mers, as well as
running time of :program:`IPK`, **grows exponentially** with ``k``. For these reasons, 
working values of ``k`` should be small (see :ref:`discussion<recommendations>` below).

Score threshold
~~~~~~~~~~~~~~~
Score threshold determines how low probability of a considered k-mer can be, or 
**how picky we are considering k-mers**. All k-mers that get a score below the specific threshold value
are excluded from consideration. The threshold is computed according to the following formula:


.. math::
    \text{minimal score} = {\Big(\frac{\omega}{\sigma}\Big)}^k


where :math:`\omega` is ``omega``, and :math:`\sigma` is the alphabet size (4 for DNA).

Thus, increasing ``omega`` leads to excluding low-scored k-mers from consideration. 

.. note::

   Note that the number of k-mers considered decreases non-linearly with increasing ``omega``.

.. _recommendations:

Recommendations
~~~~~~~~~~~~~~~

You may need to experiment with tweaking these parameters for your data. 
The instruction below may be helpful.



.. tabs::

   .. tab:: DNA

      :option:`-k`: 
         1. For short genetic markers (<2Kbp), leave it to be default 10. 
         If the tree is large (e.g., >10K taxa), the analysis may be too slow. Then, the value of 8
         can be an acceptable trade-off between accuracy and speed.

         2. For long viral genomes (>10Kb), use the value of 12. If it takes too long, consider
         using higher values of :option:`--omega`.

         3. For any kind of data, do not use values lower than 8.

      :option:`-o` (:option:`--omega`): 
         1. For large trees, use the default value of 1.5.

         2. For small trees (<1K taxa), consider using values 1.25 and 1.0 for higher accuracy.


   .. tab:: Proteins

         Temporarily disabled in this version of IPK.