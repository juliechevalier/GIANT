GIANT
=====

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3611057.svg
   :target: https://doi.org/10.5281/zenodo.3611057

Galaxy-based Interactive tools for ANalysis of Transcriptomic data


User-friendly tools suite for micro-arrays analyses and for exploring RNA-seq & Micro-Arrays differential results
=================================================================================================================

Overview of GIANT
-----------------
GIANT is a modular Galaxy tools suite interfacing R packages (as Limma) allowing users to perform transcriptomic analyses in a user-friendly environment. With GIANT, users can consecutively perform all required steps for microarray-based analyses ranging from data QC to differential analyses and complex visualization and unsupervised explorations. Furthermore, specific modules such as QC-plot, Volcano-plot, GSEA-format or Heatmap modules can be used with any type of pre-processed data (as RNA-seq normalized data or RNA-seq differential results). Finally, we use the power of Plotly and Heatmaply packages to provide interactive graphics and interactive requestable result tables in each module, a very innovative feature for Galaxy tools. 


Details of the modules of GIANT
-------------------------------
Our tools-suite wraps functionalities for: 1) QC plot generation (Array images, boxplots, signal density plots, MA plots, PCA), data normalization and annotation (with APT tool from Affymetrix, with probe-set and gene-level analysis), 2) complex differential analyses with no restriction on number of conditions and factors and taking into account blocking effects (as batch effects or paired-analysis), 3) data visualization with raw p-values histograms, volcanos and circular plots and 4) heatmap and hierarchical clustering generation based on normalized expressions and on differential analysis results.
Interactive graphics integrated in our modules allow users, for example : to explore data within 3D PCA plots, to target particular points in the volcano plot to easily identify genes of interest or to zoom in heatmaps to explore the content of clusters


Citation
========

Zenodo. http://doi.org/10.5281/zenodo.1477870

Each tool wraps R/Bioconductor packages or programs, cite them as indicated in the citation section of the corresponding wrapper

Documentation
=============

See help sections in wrappers & R packages documentation

Tutorial slides and Troubleshooting are available on Documentation directory


Galaxy Installation
===================
see testtoolshed : correct version with owner : vandelj

toolshed version soon

------------

This tools suite is developed by the Team_ "Molecular analysis of gene regulation in cardiometabolic diseases" at U1011_-Nuclear receptors, cardiovascular diseases and diabetes. U1011_ is a french mixed research unit of the University of Lille, the Inserm, the Universitary Hospital of Lille (CHU) and the Institut Pasteur de Lille.

.. _Team: https://u1011.pasteur-lille.fr/lunite/theme-4-analyse-moleculaire-de-la-regulation-des-genes-dans-le-syndrome-cardiometabolique/

.. _U1011: http://u1011.pasteur-lille.fr/accueil/

Tools compatibilities
===================

To insure full compatibility between GIANT tools, please be sure to use compatible tool versions as described in the table below. Usually, minor changes generate increment of the version third number (e.g v0.3.1 to v0.3.2). Major changes are indicated through increment of the version second number (e.g v0.2.2 to v0.3.0).
More generally, to avoid compatibility issues, try to get the most up to date version of each tool from the Galaxy tool shed.

+------------------------+----------------------+--------------------+----------------------+--------------------+-------------------------------+-------------------+
| giant_factor_generator | giant_plot_functions | giant_aptsummarize | giant_limma_analysis | giant_volcano_plot | giant_hierarchical_clustering | giant_gsea_format |
+========================+======================+====================+======================+====================+===============================+===================+
| v0.1.1                 | v0.1                 | v0.1               | v0.1                 | NA                 | v0.1                          | v0.1              |
+------------------------+----------------------+--------------------+----------------------+--------------------+-------------------------------+-------------------+
| v0.1.1                 | v0.1.1               | v0.1               | v0.2.0               | v0.1.0             | v0.1                          | v0.1              |
+------------------------+----------------------+--------------------+----------------------+--------------------+-------------------------------+-------------------+
| v0.1.1                 | v0.1.1               | v0.1               | v0.3.0               | v0.2.0             | v0.2.0                        | v0.2.0            |
+------------------------+----------------------+--------------------+----------------------+--------------------+-------------------------------+-------------------+
| v0.1.1                 | v0.1.2               | v0.1.1             | v0.3.6               | v0.2.3             | v0.4.0                        | v0.2.0            |
+------------------------+----------------------+--------------------+----------------------+--------------------+-------------------------------+-------------------+

