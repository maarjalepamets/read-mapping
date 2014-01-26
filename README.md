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
Directory "example/" contains the following files that can be used for testing the application:
* 
* 
* 
* 

### Indexer

### Mapper

