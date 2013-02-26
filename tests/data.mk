# HMP Sample anterior nares
ifdef SRS018585
KMIN=27
KMAX=31
STEPSIZE=2
OUT=test-out/masm-SRS018585
FASTQ1=data/SRS018585/SRS018585.denovo_duplicates_marked.trimmed.1.fastq
FASTQ2=data/SRS018585/SRS018585.denovo_duplicates_marked.trimmed.2.fastq
FASTQBASE=SRS018585
DO_QTRIM=no # any value other than yes, seqs have already been trimmed

# Download the test data
$(FASTQ1):
	mkdir -p data
	cd data && \
	wget http://downloads.hmpdacc.org/data/Illumina/anterior_nares/SRS018585.tar.bz2 && \
	tar -xjf SRS018585.tar.bz2
$(FASTQ2): $(FASTQ1)
	if [ -f $@ ]; then \
		touch $@; \
	else \
		exit 1; \
	fi
endif
