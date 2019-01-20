# RepeatResolver
Data and Code for the paper "Resolving repeat families with long reads"

This repository contains the scripts and tools used to compute the results presented in “Resolving repeat families with long reads”. 

In this paper I describe a pipeline for the resolution of interspersed repeats in long read genome assemblies. In this pipeline we build a multiple sequence alignment out of the repetitive sections of the reads, refine this multiple sequence alignment, and conduct a statistical analysis to extract copy differences to then use these copy differences in a clustering algorithm to pull the different repeat copy groups apart. We then use several clusterings of sections of the multiple sequence alignment to resolve a large repetitive region by computing the most likely connections of the unique sequences on one side of the repeat to the unique sequences on the other side. 
This pipeline is then assessed with simulated data and Drosophila melanogaster transposons.

The script used to create the simulated data sets is included. 

The transposon data sets can be found at: https://github.com/PhilippBongartz/DrosophilaHistoneComplex/tree/master/Transposons
Under this link the files MidTransposonMMA_x.tar.bz2 contain MSAs for all transposon data sets, while the files TransposonCopies_x contain the respective ground truth information. 

In the following we present the scripts and tools (both in **bold**) in the order in which they have to be run to recreate the benchmarking of the pipeline via simulated data.



**DataSimulator.py**

DataSimulator has several flags to adjust the parameters of the simulated data set:
-c <40>      determines the sequencing coverage.
-n <100>     can be used to choose the number of copies in the repeat family.
-d <1>       the percentage minimal differences between copies.
-l <30000>   the repeat length, which will be bolstered by 10kbp unique flanking sequences.
-t <Tree>    Tree, Distributed, EquiDistant - the three choices of repeat structure types.

python DataSimulator.py -c 40 -n 100 -d 1 -l 30000 -t Tree
would be the command identical to the default parameter values and creates one of the data sets analysed in the paper.

The output consists of four data sets, containing the reads, the readplacement information, the read copy information and the repeat template. These files are named by incorporating repeat type, percentage of mininmal differences and the repeat length into a data set name, for example the default would be Tree_1perc_30000kb + .fasta, _ReadPlacements, _ReadCopynumbers and _Template.fasta. 

DataSetName_Template.fasta
contains the underlying repeat template. From this template all repeat copies are sampled, by introducing single nucleotide polymorphisms.

DataSetName_ReadCopynumbers
Contains the number of the copy from which each read has been sampled.

DataSetName_ReadPlacements
Contains the starting place in the copy from which the read has been sampled.

DataSetName.fasta
Contains the reads sampled from the copies.




TOOLS:

**ReadCutter.c**

Is a tool that takes the reads and the template as input and cuts the reads into instances of the repeat and instances of flanking sequences. It outputs DataSetName_Seq.fasta, which contains the sequences into which the reads have been cut. And DataSetName_ReadSeqInfo which contains the information into which sequences each read has been cut.

gcc ReadCutter.c -o ReadCutter -lm
./ReadCutter DataSetName_Template.fasta DataSetName.fasta




**InitialAligner.c**

The InitialAligner creates a multiple sequence alignment by aligning the DataSetName_Seq.fasta to the DataSetName_Template.fasta. It outputs a multiple sequence alignment DataSetName_MSA and the information which sequences where aligned and which didn’t fit the template, DataSetName_ReadSeqInfo. It can be run in parallel, which is recommended. 

gcc InitialAligner.c -o InitialAligner -lpthread
./InitialAligner DataSetName_Template.fasta DataSetName_Seq.fasta -p <number of available cores>




**PW_ReAligner.c**

PW_ReAligner refines the initial multiple sequence alignment by optimising the sum-of-pairs-score. It outputs by default MSAreal, which can be changed with -o DataSetName_MSAreal. This is the computational bottleneck, as recursive realigning is difficult to parallelize. Running it over a weekend is recommended. Full convergence may not be absolutely necessary, the current result is written out after each realignment round and can be used even while the MSA is further refined. 

gcc -mcmodel=medium PW_ReAligner.c -o PW_ReAligner
./PW_ReAligner DataSetName_MSA -o DataSetName_MSAreal





**Window.py**

Calculates coverage for different sites of the multiple sequence alignment and provides a subdivision of the MSA into equally spaced sections with a certain coverage. The flag -c <0.90> determines the fraction of the average coverage of the MSA that should still be available at each point in a section. The flag -p <6> determines into how many sections the MSA will be divided. The output is a list of MSA sites that can be used to subdivide the MSA. For the simulated 30kbp repeats we have used 6 parts a la 5kbp. For the transposon MSA 1 part would allow to recreate the single-step resolution. 
python Window.py -c <0.90> -p 6





**MaxCorrelation.c**

MaxCorrelation calculates the statistical significance of the deviations from the majority base in each column of a given multiple sequence alignment. This allows us to detect differences between repeat copies despite the high error rate of the reads. The input is a refined MSA “DataSetName_MSAreal”, the output is a list “MaxCorrsOf_DataSetName_MSAreal” of the maximal statistical significance of all intersections with all base groups at all sites for each base group at each site. The flag -p <1> allows parallel computation. The flag -c <30> is used to determine coverage, which is used to discard base groups that are two small to belong to even a single copy. MaxCorrelation uses the gnu scientific library. 

gcc -I/usr/local/include -c MaxCorrelation.c
gcc -L/usr/local/lib MaxCorrelation.o -lgsl -lgslcblas -lm -o MaxCorrelation -lpthread
./MaxCorrelation DataSetName_MSAreal -c <30> -p <number of available cores>





**RepeatResolver.c**

The RepeatResolver creates a clustering of a given multiple sequence alignment using the statistically significant differences within the specified window. Input is the refined MSA “DataSetName_MSAreal”. It also reads in the list of significant differences “MaxCorrsOf_DataSetName_MSAreal” that fits the MSA path. The flag -p <1> allows parallel computation. The flag -c <30> is used to determine coverage, it can be chosen below the actual coverage, which leads to a clustering with smaller clusters. The flag -f <start end> specifies a window for resolution. Only the sequences that cover all the columns between start and end are clustered and only the differences in the columns between start and end are used for the clustering.
Start and end are should be consecutive pairs of MSA sites taken from the output of Window.py. 
So if the output of Window.py is “x y z”, we run RepeatResolver once with -f x y and once with -f y z. 
RepeatResolver outputs three clusterings. 
KmeansSubdivisionOf_start_end_DataSetName_MSAreal contains the final clustering that used a kmeans variant to further subdivide the earlier clustering. The earlier clusterings are also written out as DropoffSubdivisionOf_start_end_DataSetName_MSAreal and RelDropSubdivisionOf_start_end_DataSetName_MSAreal, were the former is the clustering using only the subdivision by refined groups and the later is the clustering using the recursive subdivision. 

gcc -I/usr/local/include -c RepeatResolver.c
gcc -L/usr/local/lib RepeatResolver.o -lgsl -lgslcblas -lm -o RepeatResolver -lpthread
./RepeatResolver DataSetName_MSAreal -c <30> -f start end -p <number of available cores>





**SimDataAssessment.py**
Calculates the results for a simulated data set. SimDataAssessment.py assumes that all files belonging to the data set are in the working directory, while the folder containing the calculated clusterings is given as a command line parameter. It detects and reads in all files and calculates the number of resolved copies for each clustering as single-step resolution and for the overall repeat by using all clusterings. 

python SimDataAssessment.py /Clusterings




**TransposonAssessment.py**
Calculates the results for a transposon data set. The Kmeans clustering path is given as a command line parameter and the other files are loaded depending on that. TransposonAssessment.py expects the MSA, the MaxCorrs-file, the three Clusterings, as well as the ground truth files “TransposonCopies_x” to be in the working directory. 

TransposonAssessment.py KmeansSubdivisionOf_start_end_MidTransposonMMA_x_real



