Installation guide
==================

The ``scmmib`` package can be installed via following 2 steps:

1. Preparing the envrionment. 
-  **Option 1:**/  install minimum dependencies by pip.

.. code-block:: bash
    :linenos:

    pip install scib scglue scanpy


-  **Option 2 (recomended):**  install all dependencies with conda.
  
.. code-block:: bash
    :linenos:

    conda env create -f scmmib_env.yml
    conda activate scmmib


2. Install scmmib package.
   
.. code-block:: bash
    :linenos:

    # download SCMMIB
    git clone https://github.com/bm2-lab/SCMMI_Benchmark
    # set dir to folder
    cd SCMMI_benchmark
    pip install .


3. Test the installation in python
 
.. code-block:: python
    :linenos:

    import scmmib
    scmmib.__version__

.. note:: 
 - A bug may occur for graph LISI metrics as follows:

    `"FileNotFoundError, [Errno 2] No such file or directory: '/tmp/lisi_svo3el2i/graph_lisi_indices_0.txt'""`

    The related GitHub issue in scib project is `here <https://github.com/theislab/scib/issues/333>`__ and a posssible `solution <https://github.com/theislab/scib/blob/main/scib/knn_graph/README.md>`__ .
