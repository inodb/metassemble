SRS018585=letsgo 
include data.mk

include /gulo/glob/inod/github/metassemble/scripts/Makefile

define diff_files
@echo TEST $1
# Compare each output in test-out with the corresponding output in ref-out
for file in $^; do \
	diff -q $$file `echo $$file | sed 's/test-out/ref-out/'`; \
done
touch $@
@echo TEST $1 - SUCCEEDED
endef

define check_exist
@echo TEST $1
# Check if each file exists. This test is used if output can vary, when the
# program is multi-threaded for instance
for file in $^; do \
	if [ ! -f $$file ]; then \
		echo Output $$file does not exist; \
		exit 1; \
	fi; \
done
touch $@
@ echo TEST $1 - SUCCEEDED
endef


$(OUT).qtrim.done: $(FASTQ_TRIM_IL)
	$(call check_exist,QTRIM)

$(OUT).velvet.done: $(VELVETG_OUT_NOSCAF) $(VELVETG_OUT_SCAF)
	$(call check_exist,VELVET)

$(OUT).metavelvet.done: $(METAVELVETG_OUT_NOSCAF) $(METAVELVETG_OUT_SCAF)
	$(call check_exist,METAVELVET)

$(OUT).ray.done: $(RAY_CONTIGS_OUT) $(RAY_SCAFFOLDS_OUT)
	$(call check_exist,RAY)

$(OUT).minimus2.done: $(MINIMUS2_OUT_VELVET_NOSCAF)/$(MERGE_FILENAME)
	$(call check_exist,MINIMUS2)

$(OUT).newbler.done: $(NEWBLER_OUT_VELVET_NOSCAF)/$(MERGE_FILENAME)
	$(call check_exist,NEWBLER)

$(OUT).bambus2.done: $(foreach contigs, $(filter %$(KMAX)/$(CONTIG_FILENAME),$(VELVETG_OUT_NOSCAF)), $(shell dirname $(contigs))/bambus2/bambus2.scaffold.linear.fasta)
	$(call check_exist,BAMBUS2)

testqtrim: $(OUT).qtrim.done
testvelvet: $(OUT).velvet.done
testmetavelvet: $(OUT).metavelvet.done
testminimus2: $(OUT).minimus2.done
testnewbler: $(OUT).newbler.done
testray: $(OUT).ray.done

TESTDONEFILES=$(OUT).velvet.done $(OUT).metavelvet.done $(OUT).ray.done $(OUT).minimus2.done $(OUT).newbler.done $(OUT).bambus2.done
# only add qtrim if qtrim is performed
ifeq ($(DO_QTRIM),yes)
$(TESTDONEFILES)=$(OUT).qtrim.done $(TESTDONEFILES)
endif

test: $(TESTDONEFILES)
clean:
	rm -f $(TESTDONEFILES)
.PHONY: testqtrim testvelvet testmetavelvet testminimus2 testnewbler testray test clean
