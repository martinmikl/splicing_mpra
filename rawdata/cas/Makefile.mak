SEQUENCE_LENGTH = 140

##########################################################################################
### 			  Mapping reads from 300PE sequencing runs		       ###
##########################################################################################


concatenate_reads:
	cat *R1_001.fastq \
	| paste - - - - \
	| awk '{OFS="\t"}{print $$1 "_" $$2, $$3, $$4, $$5}' \
	> all_reads1.tab; \
	cat *R2_001.fastq \
	| paste - - - - \
	| awk '{OFS="\t"}{print $$1 "_" $$2, $$3, $$4, $$5}' \
	> all_reads2.tab; \
	
### Trim read length in order to increase accuracy of mapping by barcodes only ###

trim_reads:
	cat all_reads1.tab \
	| awk '{OFS="\t"}{if (length($$2)>=$(SEQUENCE_LENGTH)) print $$1, substr($$2,1,$(SEQUENCE_LENGTH)), $$3, substr($$4,1,$(SEQUENCE_LENGTH)); else print $$1, $$2, $$3, $$4}' \
	> trimmed_reads1.tab; \
	cat trimmed_reads1.tab \
	| tr '_' ' ' \
	| tr '\t' '\n' \
	> trimmed_reads1.fastq; \
	cat all_reads2.tab \
	| awk '{OFS="\t"}{if (length($$2)>=$(SEQUENCE_LENGTH)) print $$1, substr($$2,1,$(SEQUENCE_LENGTH)), $$3, substr($$4,1,$(SEQUENCE_LENGTH)); else print $$1, $$2, $$3, $$4}' \
	> trimmed_reads2.tab; \
	cat trimmed_reads2.tab \
	| tr '_' ' ' \
	| tr '\t' '\n' \
	> trimmed_reads2.fastq; \

split_reads:
	cat trimmed_reads1.fastq \
	| split -l 200000 -d -a 3 - Split1-; \
	cat trimmed_reads2.fastq \
	| split -l 200000 -d -a 3 - Split2-; \

concatenate_trim_split:
	make concatenate_reads
	make trim_reads
	make split_reads
	
run_python_mapping:
	$(eval num_files := $(shell ls -l | grep -o "Split1-.*" | wc -l))
	$(eval num_list := $(shell numgen.pl -s 0 -z 3 -n $(num_files)))
	echo -n "" \
	$(foreach num, $(num_list), \
		python2.7 ../../code/functions/mapcasrna.py $(num)
	)\


