.. _docs-quickstart:

Try it out now!
============================================================

*diffTF* runs on Linux and macOS and is even independent on the operating system if combined with ``Singularity``. The following quick start briefly summarizes the necessary steps to install and use it.

Principally, there are two ways of installing *diffTF* and the proper tools:

1a. **The "easy" way**: Using ``Singularity`` and our preconfigured *diffTF* containers that contain all necessary tools, R, and R libraries

  You only need to install Snakemake (see below for details) and ``Singularity``. *Snakemake* supports Singularity in Versions >=2.4. You can check whether you already have ``Singularity`` installed by simply typing

  .. code-block:: Bash

    singularity --version

  Snakemake requires at least version 2.4. If your version is below, please update to the latest ``Singularity`` version.

  .. note:: Make to read the section :ref:`docs-singularityNotes` properly!

1b. **The "more complicated" way**:  Install the necessary tools (*Snakemake*, *samtools*, *bedtools*, *Subread*, and *R* along with various packages).

  .. note:: Note that all tools require Python 3.

  We recommend installing all tools except R via conda, in which case the installation then becomes as easy as

  .. code-block:: Bash

    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda install snakemake bedtools samtools subread

  If conda is not yet installed, follow the `installation instructions <https://conda.io/docs/user-guide/install/index.html>`_. Installation is quick and easy. Make sure to open a new terminal after installation, so that *conda* is available.

  .. note:: You do not need to uninstall other Python installations or packages in order to use conda. Even if you already have a system Python, another Python installation from a source such as the macOS Homebrew package manager and globally installed packages from pip such as pandas and NumPy, you do not need to uninstall, remove, or change any of them before using conda.

  If you want to install the tools manually and outside of the conda framework, see the following instructions for each of the tools: `snakemake  <http://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_, `samtools <http://www.htslib.org/download>`_, `bedtools <http://bedtools.readthedocs.io/en/latest/content/installation.html>`_, `Subread <http://subread.sourceforge.net>`_.

  In addition, *R* is needed along with various packages (see below for details).

2. **Clone the Git repository:**

    .. code-block:: Bash

      git clone https://git.embl.de/grp-zaugg/diffTF.git

    If you receive an error, *Git* may not be installed on your system. If you run Ubuntu, try the following command:

    .. code-block:: Bash

      sudo apt-get install git

    For macOS, there are multiple ways of installing it. If you already have *Homebrew* (http://brew.sh) installed, simply type:

    .. code-block:: Bash

      brew install git

    Otherwise, consult the internet on how to best install Git for your system.

3. **To run diffTF with an example ATAC-Seq / RNA-seq dataset for 50 TF, simply perform the following steps (see section**  :ref:`exampleDataset` **for dataset details)**:

  * Change into the ``example/input`` directory within the Git repository

      .. code-block:: Bash

        cd diffTF/example/input

  * Download the data via the download script

        .. code-block:: Bash

          sh downloadAllData.sh

  * To test if the setup is correct, start a dryrun via the first helper script

        .. code-block:: Bash

          sh startAnalysisDryRun.sh

  * Once the dryrun is successful, start the analysis via the second helper script.

    .. code-block:: Bash

      sh startAnalysis.sh

    If you want to include ``Singularity`` (which we strongly recommend), simply edit the file and add the ``--use-singularity`` and ``--singularity-args`` command line arguments in addition to the other arguments (see the Snakemake documentation and the section :ref:`docs-singularityNotes` for more details).

    Thus, the command you execute should look like this:

        .. code-block:: Bash

          snakemake --snakefile ../../src/Snakefile --cores 2 --configfile config.json \
           --use-singularity --singularity-args "--bind /your/diffTF/path"

    Read in section :ref:`docs-singularityNotes` about the ``--bind`` option and what ``/your/diffTF/path`` means here , it is actually very easy!

    You can also run the example analysis with all TF instead of only 50. For this, simply modify the ``TF`` parameter and set it to the special word ``all`` that tells *diffTF* to use all recognized TFs instead of a specific list only (see section :ref:`parameter_TFs` for details).

4. **To run your own analysis**, modify the files ``config.json`` and ``sampleData.tsv``. See the instructions in the section `Run your own analysis`_ for more details.
5. **If your analysis finished successfully**, take a look into the ``FINAL_OUTPUT`` folder within your specified output directory, which contains the summary tables and visualization of your analysis. If you received an error, take a look in Section :ref:`docs-errors` to troubleshoot.

.. _docs-prerequisites:

Prerequisites for the "easy" way
==================================

The only prerequisite here is that Snakemake and ``Singularity`` must be installed on the system you want to run *diffTF*. See above for details with respect to the supported versions etc. For details how to install Snakemake, see below.


Prerequisites for the "manual" way
=====================================

Note that most of this section is only relevant if you use Snakemake without ``Singularity``. This section lists the required software and how to install them. As outlined in Section :ref:`docs-quickstart`, the easiest way is to install all of them via ``conda``. However, it is of course also possible to install the tools separately.

Snakemake
--------------------------

Please ensure that you have at least version 5.3 installed. Principally, there are `multiple ways to install Snakemake <http://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_. We recommend installing it, along with all the other required software, via conda.

*samtools*, *bedtools*, *Subread*
----------------------------------

In addition, `samtools <http://www.htslib.org/download>`_, `bedtools <http://bedtools.readthedocs.io>`_ and `Subread <http://subread.sourceforge.net>`_ are needed to run *diffTF*. We recommend installing them, along with all the other required software, via conda.


R and R packages
--------------------------

A working ``R`` installation is needed and a number of packages from either CRAN or Bioconductor have to be installed.  Type the following in ``R`` to install them:

.. code-block:: R

  install.packages(c("checkmate", "futile.logger", "tidyverse", "reshape2", "RColorBrewer", "ggrepel", "lsr", "modeest", "boot", "grDevices", "pheatmap", "matrixStats", "locfdr"))

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install(c("limma", "vsn", "csaw", "DESeq2", "DiffBind", "geneplotter", "Rsamtools", "preprocessCore", "apeglm"))


.. _docs-runOwnAnalysis:

Run your own analysis
============================================================

Running your own analysis is almost as easy as running the example analysis (see section :ref:`exampleDataset`). Carefully read and follow the following steps and notes:

1. Copy the files ``config.json`` and ``startAnalysis.sh`` to a directory of your choice.
2. Modify the file ``config.json`` accordingly. For example, we strongly recommend running the analysis for all TF instead of just 50 as for the example analysis. For this, simply change the parameter “TFs” to “all”. See Section :ref:`configurationFile` for details about the meaning of the parameters. Do not delete or rename any parameters or sections.
3. Create a **tab-separated** file that defines the input data, in analogy to the file ``sampleData.tsv`` from the example analysis, and refer to that in the file ``config.json`` (parameter ``summaryFile``)
4. Adapt the file ``startAnalysis.sh`` if necessary (the exact command line call to Snakemake and the various Snakemake-related parameters). If you run with Singularity, see the section below for modifications.
5. Since running the pipeline is often computationally demanding, read Section :ref:`timeMemoryRequirements` and decide on which machine to run the pipeline. In most cases, we recommend running *diffTF* in a cluster environment (see Section :ref:`clusterEnvironment` for details). The pipeline is written in Snakemake, and we strongly suggest to also read Section :ref:`workingWithPipeline` to get a basic understanding of how the pipeline works.


.. _docs-singularityNotes:

Adaptations and notes when running with Singularity
============================================================
 With ``Singularity``, each rule will be executed in pre-configured isolated containers that contain all necessary tools.  To enable it, you only have to add the following arguments when you execute Snakemake:

1. ``--use-singularity``: Just type it like this!

2. ``--singularity-args``: You need to make all directories that contain files that are referenced in the *diffTF* configuration file available within the container also. By default, only the directory and subdirectories from which you start the analysis are automatically mounted inside the container. Since the *diffTF* source code is outside the ``input`` folder for the example analysis, however, at least the root directory of the Git repository has to be mounted. This is actually quite simple! Just use ``--singularity-args "--bind /your/diffTF/path"`` and replace ``/your/diffTF/path`` with the root path in which you cloned the *diffTF* Git repository (the one that has the subfolders ``example``, ``src`` etc.). If you reference additional files, simply add one or multiple directories to the bind path (use the comma to separate them). For example, if you reference the files ``/g/group1/user1/mm10.fa`` and ``/g/group2/user1/files/bla.txt`` in the configuration file file, you may add ``/g/group1/user1,/g/group2/user1/files`` or even just ``/g`` to the bind path (as all files you reference are within ``/g``).

  .. note:: We note again that within a Singularity container, you cannot access paths outside of the directory from where you started executing Snakemake. If you receive errors in the ``checkParameterValidity`` rule that a directory does not exist even though you can cd into it, you most likely forgot to include the path this folder or a parent path as part of the ``bind`` option.

3. ``--singularity-prefix /your/directory`` (optional): You do not have to, but you may want to add the ``--singularity-prefix`` argument to store all ``Singularity`` containers in a central place (here: ``/your/directory``) instead of the local ``.snakemake`` directory. If you intend to run multiple *diffTF* analyses in different folders, you can save space and time because the containers won't have to be downloaded each time and stored in multiple locations.

Please read the following additional notes and warnings related to ``Singularity``:

- .. warning:: If you use ``Singularity`` version 3, make sure you have at least version 3.0.3 installed, as there was an issue with Snakemake and particular ``Singularity`` versions. For more details, see `here <https://bitbucket.org/snakemake/snakemake/issues/1017/snakemake-process-suspended-upon-execution>`_.
