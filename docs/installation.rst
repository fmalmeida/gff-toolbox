.. _installation:

Installation
============

This package can be used without installation or installed via pip or conda. Remember, when not installing it via conda you must make sure that the dependencies are properly installed.

Dependencies
------------

* Python >= 3.6, and its packages:

  * `pandas <https://pandas.pydata.org/>`_
  * `biopython <https://biopython.org/>`_
  * `bcbiogff <https://github.com/chapmanb/bcbb/tree/master/gff>`_ 
  * `dna features viewer <https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer>`_
  * `matplotlib <https://matplotlib.org/>`_
  * `docopt <http://docopt.org/>`_
  * `pymongo <https://pypi.org/project/pymongo/>`_
* `mongo shell <https://docs.mongodb.com/manual/>`_

From source
-----------

.. code-block:: bash

   # Download
    git clone https://github.com/fmalmeida/gff-toolbox.git
    cd gff-toolbox

    # Run without installing
    python3 gfftoolbox-runner.py -h

    # Install and run in any place
    python3 setup.py install
    gff-toolbox -h

Via conda
---------

.. code-block:: bash

    # Get the conda package
    conda create -n gff-toolbox -c conda-forge -c bioconda -c anaconda -c defaults -c falmeida gff-toolbox

.. tip:: Users are advised to use mamba.
