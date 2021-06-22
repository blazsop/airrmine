# ------------------------------------------------------------------------------
# ImmChainTracer   v1.00
# Germline sequence tracker pipeline in adaptive immune receptor repertoires
#
# by Peter Blazso
# 2021-06-21
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# GLOBAL SETTINGS
# ------------------------------------------------------------------------------
SHELL     := /bin/bash   # DO NOT EDIT! Using "bash" is mandatory!

# unique name/label of the processing pipeline
LABEL     :=

# immune receptor chain type
CHAIN     := IgH

# order of repositories for the extraction of germline-related sequences
REPORDER  := MN AN MB NE

# source/target directories
SCRIPTDIR := scripts
IMGTDIR   := imgt
IREPDIR   := irep
TSVDIR    := tsv
PREPDIR   := prep
DESTDIR   := dest
FILTDIR   := filtered
CLONEDIR  := clone
LINDIR    := lin

# filename suffices
GL_SUFF   := _germlines_$(CHAIN).tsv
GL_P_SUFF := _germlines_$(CHAIN)_prepared.tsv
IREP_SUFF := _$(CHAIN)_irep.csv

# source/target files
#IMGT      = $(wildcard $(IMGTDIR)/*.txz)     # individual IMGT output files
#TSV       = $(wildcard $(TSVDIR)/*.tsv) 
#TSVnoGL   = $(filter-out %$(GL_SUFF), $(TSV))
PREP      = $(wildcard $(PREPDIR)/*.tsv)
PREPnoGL  = $(filter-out %$(GL_P_SUFF), $(PREP))
REPS      = $(subst _collapsed,,$(notdir $(basename $(PREPnoGL)))) # raw repertoire names
#GLS       = $(wildcard $(TSVDIR)/*$(GL_SUFF))
#GLSBASE   = $(subst $(GL_SUFF),,$(notdir $(GLS)))
GLFPREP   = $(shell gawk -f $(SCRIPTDIR)/glrellist.gawk $(PREPnoGL))
GLFS      = $(basename $(GLFPREP))            # raw germline-filtered names
LINS      = $(shell gawk -f $(SCRIPTDIR)/linlist.gawk $(PREPDIR)/*$(GL_P_SUFF))
LINEAGES  = $(addprefix $(LINDIR)/, $(addsuffix .tsv, $(LINS)))

#PREPGLS   = $(addprefix $(PREPDIR)/, $(addsuffix $(GL_P_SUFF), $(GLSBASE)))
DISTTH    = $(addprefix $(DESTDIR)/, $(addsuffix _dist_thrs.data, $(REPS)))
FILTREPS  = $(addprefix $(FILTDIR)/, $(addsuffix _filtered.list, $(REPS)))
CLONES    = $(addprefix $(CLONEDIR)/, $(addsuffix _clone.tsv, $(GLFS)))

# ------------------------------------------------------------------------------
.PHONY: all clean

all: filtreps clones lineages addcopy 
	@echo done.

# ------------------------------------------------------------------------------
# 0. CONVERT (EXTRACT) IMGT annotated DB files into TSV databases (optional)
# ------------------------------------------------------------------------------
#$(TSVDIR)/%_db-pass.tsv: $(IMGTDIR)/%.txz
#	@echo Converting IMGT database $< into TSV file...
#	@mkdir -p $(TSVDIR)
#	@MakeDb.py imgt -i $< --outdir $(TSVDIR)

# ------------------------------------------------------------------------------
# ENTRY POINT
# 1. CALCULATE Hamming distances, distance thresholds and generate histograms
# ------------------------------------------------------------------------------
#$(DISTTH) : $(DESTDIR)/%_dist_thrs.data: $(TSVDIR)/%.tsv
#	@echo Calculating Hamming distances and determine distance threshold on $< ...
#	@mkdir -p $(DESTDIR)
#	@Rscript --vanilla $(SCRIPTDIR)/findthreshold.R $< $@

# ------------------------------------------------------------------------------
# 2. PREPARE GERMLINE databases for subsequent analysis
# ------------------------------------------------------------------------------
#$(PREPGLS) : $(PREPDIR)/%_prepared.tsv : $(TSVDIR)/%.tsv 
#	@echo Preparing germline database from $* ...
#	@mkdir -p $(PREPDIR)
#	@Rscript --vanilla $(SCRIPTDIR)/prepgls.R $< $@

# ------------------------------------------------------------------------------
# 3. FILTER AND DISSECT repertoires and COMBINE them with germline sequences
# ------------------------------------------------------------------------------
$(FILTREPS) : $(FILTDIR)/%_filtered.list : $(PREPDIR)/%_collapsed.tsv 
	@echo Pre-filtering and dissecting repertoire $* ...
	@mkdir -p $(FILTDIR)
	@gawk -v outdir=$(FILTDIR) -f $(SCRIPTDIR)/multifilter.gawk $< > $@

# ------------------------------------------------------------------------------
# 4. DEFINE CLONES in filtered, germline combined repertoires
# ------------------------------------------------------------------------------
$(CLONES) : $(CLONEDIR)/%_clone.tsv : | filtreps
	@echo Defining clones in germline-related, filtered repertoire $* ...
	@mkdir -p $(CLONEDIR)
	@$(eval GLSUFFIX=$(lastword $(subst _, ,$*)))
	@$(eval DISTFILE=$(subst collapsed_$(GLSUFFIX),dist_thrs.data,$(DESTDIR)/$*))
	@$(eval DIST=$(shell cat $(DISTFILE)))
	@DefineClones.py --failed --nproc 1 --mode gene --act set --norm len --model ham --dist $(DIST) -d $(FILTDIR)/$*.tsv --outdir $(CLONEDIR)

# ------------------------------------------------------------------------------
# 5. COLLECT germline-related sequences from clonally defined sequence groups
# ------------------------------------------------------------------------------
$(LINEAGES):  | clones
	@mkdir -p $(LINDIR)
	@$(eval GL=$(lastword $(subst _, ,$(subst .tsv,,$@))))
	@$(eval SUBJ=$(firstword $(subst _, ,$(lastword $(subst /, ,$@)))))
	@$(eval SRCFILES=$(shell ls $(CLONEDIR)/$(SUBJ)_*_collapsed_$(GL)_clone-pass.tsv))
	@echo "Collecting clone-related sequences of $(GL) in subject $(SUBJ) ..."
	@gawk -f $(SCRIPTDIR)/extractgl.gawk -v seqid=$(GL) $(firstword $(SRCFILES)) > $@
	@for REP in $(REPORDER) ; \
	do \
	  for FILE in $(SRCFILES) ; \
	    do \
	      FNAME=`basename $$FILE;`; \
	      if [[ "$$FNAME" =~ "_$${REP}_" ]] ; then \
	        echo " -> processing: $$FNAME ..."; \
	        gawk -f $(SCRIPTDIR)/extractclone.gawk $$FILE >> $@; \
	      fi; \
	    done; \
	done

#tsvdb: $(TSVDB)
#prepgls: $(PREPGLS)
distth: $(DISTTH)
filtreps: $(FILTREPS)
clones: $(CLONES)
lineages: $(LINEAGES)

clean:
	@rm -rf $(DESTDIR) $(CLONEDIR) $(LINDIR)
	@echo clean.

debug:
#	@echo TSV: $(TSV)
	@echo GLS: $(GLS)
#	@echo GLSBASE: $(GLSBASE)
	@echo PREPGLS: $(PREPGLS)
	@echo DISTTH: $(DISTTH)
	@echo REPS: $(REPS)
	@echo PREP: $(PREP)
	@echo PREPnoGL: $(PREPnoGL)
	@echo GLFPREP: $(GLFPREP)
	@echo FILTREPS: $(FILTREPS)
	@echo CLONES: $(CLONES)
	@echo LINEAGES: $(LINEAGES)
