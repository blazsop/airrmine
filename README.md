AIRRMINE - Adaptive Immune Receptor Repertoire data processing system
=====================================================================

Developed by
 * Peter Blazso      - lead bioinformatician, code maintenance
 * Krisztian Csomos  - concepts and testing
 * Boglarka Ujhazi   - data management


0. Directory structure and definitions
--------------------------------------

Ready to use scripts (without compilation) are found in the following
directories:

* airrnat
  AIRR dataset analytics and network assessment wrap-around and custom
  algorithms

* immchaintracer
  Immune Receptor clonal lineage tracker and assembler toolchain

* immchaintracer_analysis
  Immune Receptor clonal lineage analysis and visualization toolset

* samples
  Demonstration datasets for testing purposes.

The same order of analytic steps are recommended since the subsequent
steps are largely dependent on the output of the previous ones.


1. System requirements (software has been tested on)
----------------------------------------------------

Operating system:
These scripts are written, used and tested in UNIX/Linux environments
(specifically on Ubuntu Linux 20).
Command line interface is sufficient to use every piece of code.

IMGT V-Quest was used for Adaptive Immune Receptor sequence alignment.
Website: http://www.imgt.org/IMGT_vquest

Command line toolset requirements:
* bash    5.1
* make    4.3 (GNU make)
* gawk    5.1 (GNU awk)
* R       4.0.4  ("Lost Library Book")
* python  2.7.18
* python3-pip   (for easy installaion of Immcantation tools)
* phylip  3.697 (PHYLogeny Inference Package)

Immcantation software packages
Website: https://immcantation.readthedocs.io
* changeo      1.1.0

R library dependencies:
* alakazam     1.1.0
* colorspace   2.0.1
* dplyr        1.0.7
* igraph       1.2.6
* openxlsx     4.2.4
* parallel     4.0.4
* readr        1.4.0
* shazam       1.0.2
* tidyr        1.1.3
* treemap      2.4.2

2. Installation guide 
---------------------
Your computer needs to be online during the whole installation process.
It takes typically 20-30 minutes on a modern multiprocessing computer.

2.1. Operating system and software packages
Please see your operating system's manual and install the command line
toolset or make sure they are available. It may need administratory 
privileges on the system. 
For example on Ubuntu Linux use apt (e.g.: apt install bash).

2.2. Immcantation framework
You can install "changeo" with the following command:
 pip install changeo

2.3. R libraries
Run "R" and at the R command prompt install the required libraries:
 install.packages(c("alakazam","colorspace","dplyr","igraph","openxlsx","parallel","readr","shazam","tidyr","treemap"))

4. Demonstration material and testing installation
--------------------------------------------------
Sample data are processed in an hour on a contemporary, multiprocessing desktop
computer. 

4.1. Create the directory enviroments
* copy "airrnat", "immchaintracer" and "immchaintracer_analysis" directories
with their contents to an empty target directory (e.g. airr)
* copy the contents of "samples" to the same target directory (e.g. airr)

4.2. Go into "airrnat" directory and use "make" to run the different steps of
the AIRRNAT pipeline or if you just simply run "make" it will run the whole
processing toolchain on the raw and pre-processed datasets.
Hint: if you have computing power and want to speed up the whole process, just
use "-j" flag to utilize more CPU cores (e.g. make -j4 will start 4 concurrent
processing threads).

If everything worked correctly, you will be able to find the following
additional directories containing results:

* "coll"    - collapsed and filtered datasets 

* "clone"   - clonally assigned AIRR datasets
  file suffix "_allclone.tsv" : clonal affiliations added
  file suffix "_withgl.tsv"   : putative germline sequence added
  file suffix "_agsel.tsv"    : BASELINe antigene selection data added
  file suffix "_rsmut.xlsx"   : R/S mutation statistics
  file suffix "_cstats.tsv"   : general compartment indices and statistics

* "dest"    - additional data on compartments
  file suffix "_genes.xlsx"   : V(D)J gene usage
  file suffix "_tmap.tiff"    : treemap of random 2500 sequences
  file suffix "_tmap_cl.tiff" : treemap of the most abundant 2500 clones

* "graph"   - AIRR lineage connectivity data and network diagrams
  file suffix "_AIRRND.svg"   : raw AIRR network diagram
  file suffix "_AIRRND_mutations.svg" : mutation "heatmap"
  file suffix "_AIRRND_itypes.svg" : color-coded isotpyes
  file suffix "_AIRRND_agsel.svg"  : heatmap of BASELINe sigma


4.3. You can run "make" in "immchaintracer" directory, as well.

Expected data directories with contents are:

* "filtered" - mini "repertoire" subsets are made based on similarity criteria
  One file corresponds to a single Adaptive Immune Receptor sequence and contains 
  all the other similar sequences that can be captured around it.

* "clone"    - clonally assigned "mini" repertoires
* "lin"      - extracted clonal lineages

4.4. If you copy "lin" directory (with contents) from step 4.3. to
"immchaintracer_analysis" and run "make" in it you can generate various plots
and statistics of single lineage trees.

5. Re-run analytics on real data
--------------------------------

5.1. Download the supplementary raw data published by our group the same way as
demo data found in "samples" directory.
5.2. Simply replace demo data with real ones in "airrnat" directory.
5.3. Start "make" in "airrnat" directory. 
5.4. When "airrnat" pipeline finishes you shall copy/move collapsed datasets
into "prep" directory in "immchaintracer".
5.5. You also need to download and put the corresponding single clone files
here, as well.
5.6. Run "make" in "immchaintracer".
5.7. Copy "lin" directory just like in 4.4. and run analytics.
