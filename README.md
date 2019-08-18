All the code is written for Python2.7, with the following packages (versions): pandas 0.20.3, numpy 1.13.1, seaborn 0.6, scipy 0.17, and sklearn 0.18.2.

Contents of the code folder:
DESIGN_* - Code used for the generation of the oligonucleotide library. DESIGN_selection_of_contexts.py starts from RNAseq data from K562 cells (produced by Thomas Gingeras' lab, available from the ENCODE website, accession ENCFF000HFA) and creates all the input files for the specific splicing type library design files.
Mapping RNA - Code used for mapping the RNAseq reads to the library variants in a splice-isoform specific way. The mapXXX.py files map the reads (available from NCBI GEO (accession GSE132064)) to the respective library variant through the unique DNA barcode, the SequenceXXX.py files map the reads to the splice isoforms.
Mapping Protein - /rawdata/five_protein/CoverageFIVE.tab and /rawdata/ir_protein/CoverageIR.tab contain the raw data of DNA reads mapped to library variants across bins, as described in the Methods section of the manuscript. The SequenceXXX.py files contain the subsequent processing and filtering of the raw data (as described in the Methods section), leading to the protein-based splicing value used in the manuscript.
ANALYSIS_COMBINED.py - This file creates dataframes and supplementary tables for each splicing type containing all the information used subsequently for analysis/plotting the figures of the manuscript and training the model.
Figure[1/2/3/4/5/6] - Produces most figures (with the exception of the ones related to the predictor) starting from Supplementary tables (S5-S8).
ML_prepare - split train and test data for all splicing types (contains randomized dropping of variants with identical sequence. If this is re-run, the model needs to be trained again on the new training sets in order to avoid overlap with the test set). 
MODEL_training - Code used for training a computational model predicting splicing ratios based on different sets of features. Features can be generated and a model can be trained on any set of input sequences and corresponding splicing value (same input is required as for the example provided below to test the predictor: a csv file with identifies, DNA sequence and splice site coordinates). 
MODEL_testing - Testing the model on the test set of the library and on data generated in other MPRAs
MODEL_interpretation - Interpretation of the models trained on the library variants
functions - collection of functions used, needs to be added to the PYTHONPATH

Data availability:
All the raw data (sequencing reads) are available in the NCBI GEO (accession GSE132064)

General instructions:
The "/code/functions" folder needs to be added to the PYTHONPATH
All code has to be run from the "code" folder
Not all the intermediate dataframes have been provided in order to not unnecessarily inflate the size of this repository, but it should be possible to generate everything from the data provided or available in public repositories.

Specific instructions for testing the predictor:
A python script is provided to test the predictor with custom input sequences. The feature generation step is computationally intensive, therefore it is recommended only test a number of sequences (or parallelize/be patient). 
The input needs to be a csv file with four columns:
column 1: any identifier
column 2: DNA sequence containing either a (potentially retained) intron, a cassette exon and the surrounding intronic sequences, tandem 5' or tandem 3' splice sites
column 3 and 4: the coordinates (0-based) of the splice sites in the DNA sequence provided in column 2, i.e. intron start (column 3) and intron end (column 4) for retained introns, exon start (column 3) and exon end (column 4) for cassette exons, donor 1 (column 3) and donor 2 (column 4) for tandem 5' splice sites and acceptor 1 (column 3) and acceptor 2 (column 4) for tandem 3' splice sites.
csv files for testing are provided in the folders "/code/data/[ir/cas/five/three]" and contain three variants from the respective test sets (generated in the file "MODEL_testing.py"). Any other input files generated based on our variants or other DNA sequences can be used as input. The only restrictions on the length of the DNA sequence is that at least 25 nt upstream of the first splice site and 25 nt downstream of the second splice sites must be included (and everything in between). However, it can be expected that the model performs best with exon/intron lengths that are close to the ones of the training set (intron lengths between 120 and 40 bp, exon lengths between 30 and 90 bp, difference between splice sites between 2 and 80 bp).
Example output files are given in the "code" folder ("newpreds[ir/cas/five/three]") and contain predictions of splicing log-ratios between the two options (spliced out vs retained for introns, included vs skipped for cassette exons, and 2nd splice site vs 1st splice site for tandem 5' and 3' splice sites) based on all feature sets. 

Obligatory input: 
	[ir/cas/five/three] depending on the splicing type of the provided sequences
	path to the csv file containing the 
Optional:
	output file name (default: predictions[.csv])

Example code (to be run from the "code" folder, with "/code/functions" in the PYTHONPATH):
python2.7 run_predictor.py ir ./data/ir/ireval_first3.csv newpredsir
python2.7 run_predictor.py cas ./data/cas/caseval_first3.csv newpredscas
python2.7 run_predictor.py five ./data/five/fiveeval_first3.csv newpredsfive
python2.7 run_predictor.py three ./data/three/threeeval_first3.csv newpredsthree