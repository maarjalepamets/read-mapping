# Index Based Genome Mapper

Genome mappers can be used for mapping short next-generation sequencing reads to reference genomes, as well as locating any sub-sequence of interest on a genome. Mapping with index based mapper is a two-step process. First, all <i>n</i>-mers in full genome must be indexed. Secondly, this index is used to locate queries/reads listed in a separate file. This way indices could be used for several times since creating an index can be very time-consuming for larger genomes.

## Build instructions

Required package: gcc 4.7 

In the root source directory:
<pre>
$ make all
</pre>

## Usage instructions

### Example files

example.tgz contains the following files that can be used for testing the application:
* pseudomonas_full_genome.fna - input for Indexer in FastA format
* pseudomonas_10.index - output of Indexer, input for Mapper
* pseudomonas.names - output of Indexer, input for Mapper
* queries - input for Mapper

### Indexer

Indexer takes a genome file(s) as a compulsory parameter, output name and wordlength are optional parameters. An example of the commandline looks as follows:
<pre>
$ ./indexer -i example/pseudomonas_full_genome.fna -o pseudomonas -n 10
</pre>
This creates the two files that are referred to as output files of the Indexer in the examples directory.

Additional help:
<pre>
$ ./indexer --help
</pre>

### Mapper

Mapper takes both Indexer's output files and a query file as an input. Number of mismatches and the length of the step are optional parameters. An example of the commandline looks as follows:
<pre>
$ ./mapper -i example/pseudomonas_10.index -g example/pseudomonas.names -q example/queries -mm 2
</pre>

This creates the output that looks as follows:
<pre>
0	gi|15595198|ref|NC_002516.1|	1765971	2	F
1	gi|15595198|ref|NC_002516.1|	1764846	0	F
2	gi|15595198|ref|NC_002516.1|	3975100	0	F
3	gi|15595198|ref|NC_002516.1|	319766	1	F
4	gi|15595198|ref|NC_002516.1|	320675	0	F
</pre>
First column indicates the number of the query, second is the name of the chromosome, third is the location on that chromosome where the query mapped, fourth is the number of mismatches/errors and the last column indicates the strand (forward (F), reverse (R)).

Additional help:
<pre>
$ ./mapper --help
</pre>


## Results

results.tgz contains the results of the validation process. We generated 25,000 random reads from the human reference genome build 37 with up to 2 mismatches and mapped them back to the genome using our Mapper. The reads are given in a file "reads" and the results in "human37.results".