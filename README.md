<h2>AIRRMINE - Adaptive Immune Receptor Repertoire data processing system</h2>

<h3>Directory structure and definitions</h3>

Ready to use scripts are found in the following directories:

* <b>airrnat</b> :
  AIRR dataset analytics and network assessment wrap-around and custom
  algorithms

* <b>immchaintracer</b> :
  Immune Receptor clonal lineage tracker and assembler toolchain

* <b>immchaintracer_analysis</b> :
  Immune Receptor clonal lineage analysis and visualization toolset

* <b>samples</b> :
  Demonstration datasets for testing purposes.

The same order of analytic steps are recommended since the subsequent
steps are largely dependent on the output of the previous ones.


<h3>1. System requirements (software has been tested on)</h3>

<h4>1.1. Operating system</h4>
These scripts are written, used and tested in UNIX/Linux environments
(specifically on Ubuntu Linux 21.04.). Command line interface is
sufficient to use every piece of code.

<h4>1.2. AIRR sequence alignments</h4>
External IMGT V-Quest (http://www.imgt.org/IMGT_vquest) service was and
can be used for Adaptive Immune Receptor sequence alignment.

<h4>1.3. Command line toolset</h4>

* bash    5.1
* make    4.3 (GNU make)
* gawk    5.1 (GNU awk)
* R       4.0.4  ("Lost Library Book")
* python  2.7.18
* python3-pip   (for easy installaion of Immcantation tools)
* phylip  3.697 (PHYLogeny Inference Package)

<h4>1.4. Immcantation software packages</h4>
Immcantation framework components and R libraries (https://immcantation.readthedocs.io)
are extensively used in the code.

* changeo      1.1.0

<h4>1.5. R library dependencies</h4>

* alakazam     1.1.0
* colorspace   2.0.1
* dplyr        1.0.7
* igraph       1.2.6
* openxlsx     4.2.4
* readr        1.4.0
* shazam       1.0.2
* tidyr        1.1.3
* treemap      2.4.2

<h3>2. Installation guide</h3>

Your computer needs to be online during the whole installation process.
It takes typically 20-30 minutes on a modern multiprocessing computer.

<h4>2.1. Operating system and command line software packages</h4>
Please see your operating system's manual and install the command line
toolset or make sure they are available. It may need administratory 
privileges on the system. 
For example on Ubuntu Linux use apt (e.g.: <code>sudo apt install bash</code>).

<h4>2.2. Immcantation framework</h4>
You can install "changeo" with the following command:
<code>pip install changeo</code>

<h4>2.3. R libraries</h4>
Run <code>R</code> and at the R command prompt install the required libraries:
<code>install.packages(c("alakazam","colorspace","dplyr","igraph","openxlsx","readr","shazam","tidyr","treemap"))</code>

<h3>3. Demonstration material and testing installation</h3>
Sample data are processed in an hour on a contemporary, multiprocessing desktop
computer. 

<h4>3.1. Create the directory enviroments</h4>

* Unzip "samples.zip".

* Copy "airrnat", "immchaintracer" and "immchaintracer_analysis" directories
with their contents to an empty target directory (e.g. airr).

* Copy the <u>contents</u> of "samples" to the same target directory (e.g. airr).

<h4>3.2. Run AIRRNAT pipeline</h4>
Go into "airrnat" directory and use `make` command to run the different steps of
the AIRRNAT pipeline or if you just simply run `make` without command line arguments
it will run the whole processing toolchain on the raw and pre-processed datasets.

<i>Hint: if you have computing power and want to speed up the whole process, just
use "-j" flag to utilize more CPU cores (e.g. <code>make -j4</code> will start 4 concurrent
processing threads).</i>

If everything worked correctly, you will be able to find the following
additional directories containing results:

* <b>coll</b> : collapsed and filtered datasets 

* <b>clone</b> : clonally assigned AIRR datasets
  - file suffix "_allclone.tsv" : clonal affiliations added
  - file suffix "_withgl.tsv"   : putative germline sequence added
  - file suffix "_agsel.tsv"    : BASELINe antigene selection data added
  - file suffix "_rsmut.xlsx"   : R/S mutation statistics
  - file suffix "_cstats.tsv"   : general compartment indices and statistics

* <b>dest</b> : additional data on compartments
  - file suffix "_genes.xlsx"   : V(D)J gene usage
  - file suffix "_tmap.tiff"    : treemap of random 2500 sequences
  - file suffix "_tmap_cl.tiff" : treemap of the most abundant 2500 clones

* <b>graph</b> : AIRR lineage connectivity data and network diagrams
  - file suffix "_AIRRND.svg"   : raw AIRR network diagram
  - file suffix "_AIRRND_mutations.svg" : mutation "heatmap"
  - file suffix "_AIRRND_itypes.svg" : color-coded isotpyes
  - file suffix "_AIRRND_agsel.svg"  : heatmap of BASELINe sigma

<h4>3.3. Run ImmChainTracer pipeline</h4>
You can run <code>make</code> in "immchaintracer" directory, as well.

These directories with processed data are created:

* <b>filtered</b> : mini, pre-filtered "repertoire" subsets are made based on
  similarity criteria. One file corresponds to a single (cloned) Adaptive
  Immune Receptor sequence and contains all the other similar sequences that
  can be found with it in a bulk AIRR sequence pool.

* <b>clone</b> : clonally assigned "mini" repertoires

* <b>lin</b> : extracted clonal lineages

<h4>3.4. Run analysis on captured lineages</h4>

If you copy "lin" directory (with contents) from step 3.3. to 
"immchaintracer_analysis" and run <code>make</code> in it you can generate
various plots and statistics of single lineage trees.

<h3>4. Re-run analytics on real data</h3>

* Download the supplementary raw data published by our group the same way as
demo data found in "samples" directory.
* Simply replace demo data with real ones in "airrnat" directory.
* Start <code>make</code> in "airrnat" directory. 
* When "airrnat" pipeline finishes you shall copy/move collapsed datasets
into "prep" directory in "immchaintracer".
* You also need to download and put the corresponding single clone files
here, as well.
* Run <code>make</code> in "immchaintracer".
* Copy "lin" directory just like in 3.4. and run analytics.

<h3>Developed by</h3>

 * Peter Blazso : lead bioinformatician, code maintenance
 * Krisztian Csomos : concepts and testing
 * Boglarka Ujhazi : data management
 * Jolan E. Walter : principal investigator
