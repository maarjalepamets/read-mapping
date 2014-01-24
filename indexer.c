/*
 * Genome Indexer
 * (using radix sort)
 *
 * Text Algorithms, 2013/2014
 *
 * Author: Maarja Lepamets
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

#include "utils.h"

int debug = 0;

/* reading in the genome in FastA format */
/* if name is not NULL it will be assigned a copy of sequence name */
unsigned readfastafile(const char *filename, wordtable *table, unsigned loc, char **name);

/* fills the table with words and their locations */
/* if name is not NULL it will be assigned a copy of sequence name */
unsigned fillwordtable(const char *data, wordtable *table, struct stat st, unsigned loc, char **name);

/* mask has as many 1's as the wordlength is */
unsigned createmask(int wordlength);

/* wrapper for in-place radix sort */
void sortwords(wordtable *table);

/* */
unsigned countunique(wordtable *table);

/* */
void findstartpositions(wordtable *table);

/* */
void sortlocations(wordtable *table);

/* for combining two sorted lists of words (with locations) */
void combineindices(wordtable *table, wordtable *other);

/* */
void writetoindex(wordtable *table, const char *outputname, sequences seqinfo);

int main (int argc, const char *argv[])
{
	int wordlen;
	int i, inputbeg = -1, inputend = -1;
	unsigned location;
	const char *outputname;
	sequences seqinfo;
	wordtable tables[2];
	wordtable *table = &tables[0];
	wordtable *temptable = &tables[1];
	char ofsname[256];
	FILE *ofs;

	memset(table, 0, sizeof(wordtable));
	memset(temptable, 0, sizeof(wordtable));

	/* default value */
	wordlen = 16;
	outputname = "output";

	/* parsing commandline arguments */
	for (i = 1; i < argc; ++i) {
		if (!strcmp(argv[i], "-i")) {
			
			if (!argv[i + 1]) break;
			/* get the locations of the input files */
			inputbeg = ++i;
			while (i < argc - 1 && argv[i + 1][0] != '-') {
				++i;
			}
			inputend = i;
		} else if (!strcmp(argv[i], "-d")) {
			debug += 1;
		} else if (!strcmp(argv[i], "-o")) {
			outputname = argv[i + 1];
			++i;
		} else if (!strcmp(argv[i], "-n")) {
			char *e;
			wordlen = strtol (argv[i + 1], &e, 10);
			if (*e != 0) {
				fprintf(stderr, "Invalid input: %s!\n", argv[i + 1]);
				exit(1);
			}
			++i;
		} else {
			fprintf(stderr, "Unknown argument: %s\n", argv[i]);
			exit(1);
		}
	}

	/* checking the parameter */
	if (wordlen > 16) {
		fprintf(stderr, "Seed size too large!\n");
		exit(1);
	}
	if (inputbeg == -1) {
		fprintf(stderr, "No input files given!\n");
		exit(1);
	}

	seqinfo.seq_name = (char **) malloc((inputend - inputbeg + 1) * sizeof(const char *));
	seqinfo.start_pos = (unsigned *) malloc((inputend - inputbeg + 1) * sizeof(unsigned));
	seqinfo.n = inputend - inputbeg + 1;

	/* work-flow */
	location = 0;
	sprintf(ofsname, "%s.names", outputname);
	ofs = fopen (ofsname, "w");
	for (i = inputbeg; i <= inputend; ++i) { /* iteration over the input file parameters */
		char *seqname;
		unsigned currentloc;
		char *p;

		if (debug > 0) fprintf (stderr, "Reading: %s...\n", argv[i]);

		temptable->wordlength = wordlen;
		if (i > inputbeg) {
			location += 10000;
		}
		currentloc = location;
		seqinfo.seq_name[i - inputbeg] = (char *) argv[i];
		seqinfo.start_pos[i - inputbeg] = location;
		location = readfastafile(argv[i], temptable, location, &seqname);
		if (temptable->nwords == 0) continue;
		sortwords(temptable);
		findstartpositions(temptable);
		sortlocations(temptable);
		if (i > inputbeg) {
			combineindices(table, temptable);
		} else {
			*table = *temptable;
		}
		memset(temptable, 0, sizeof(wordtable));
		for (p = seqname; *p; ++p) {
			if (*p <= ' ') {
				*p = 0;
				break;
			}
		}
		fprintf (ofs, "%s %s %u\n", seqname, argv[i], currentloc);
	}
	fclose (ofs);

	if (debug) {
		unsigned int i, j;
		for (i = 0; i < table->nwords - 1; ++i) {
			printf("järjestus %s, start %u\n", word2string(table->words[i], table->wordlength), table->starts[i]);
			for (j = table->starts[i]; j < table->starts[i + 1]; ++j) {
				printf("asukohad: %u ", table->locations[j]);
			}
			printf("\n");
		}
		printf("järjestus %s, start %u\n", word2string(table->words[i], table->wordlength), table->starts[i]);
		for (j = table->starts[i]; j < table->nloc; ++j) {
			printf("asukohad: %u ", table->locations[j]);
		}
		printf("\n");
	}

	writetoindex(table, outputname, seqinfo);
	fprintf(stdout, "Done!\n");
	return 0;
}

unsigned readfastafile(const char *filename, wordtable *table, unsigned locprev, char **name)
{
	struct stat st;				/* file statistics */
	int status, handle;
	const char *data;
	unsigned loc;

	/* memory-mapping a file */
	status = stat(filename, &st);
	if (status < 0) {
		fprintf (stderr, "Cannot get the statistics of file %s!\n", filename);
		exit (1);
	}
	handle = open(filename, O_RDONLY);
	if (handle < 0) {
		fprintf (stderr, "Cannot open file %s!\n", filename);
		exit (1);
	}
	data = (const char *) mmap(NULL, st.st_size, PROT_READ, MAP_SHARED, handle, 0);
	if (data == (const char *) -1) {
		fprintf (stderr, "Cannot memory-map file %s!\n", filename);
		exit (1);
	}
	close(handle);

	if (st.st_size > table->nword_slots) {
		table->nword_slots = st.st_size;
		table->nloc_slots = st.st_size;
		table->words = (unsigned *) realloc(table->words, table->nword_slots * sizeof(unsigned));
		table->locations = (unsigned *) realloc(table->locations, table->nloc_slots * sizeof(unsigned));
	}
	loc = fillwordtable(data, table, st, locprev, name);

	munmap((void *) data, st.st_size);
	return loc;
}

unsigned fillwordtable(const char *data, wordtable *table, struct stat st, unsigned locprev, char **name)
{
	unsigned word, location;
	unsigned mask;
	int m, wordlength;
	int isheader = 0;
	off_t i;
	/* Initiqalized to suppress warning */
	off_t namestart = 0;

	wordlength = table->wordlength;
	mask = createmask(wordlength);
	word = 0;
	location = locprev;
	m = 0;

	/* find the unsigned integer corresponding to a word with a given length */
	for (i = 0; i < st.st_size; ++i) {
		if (data[i] == '>') {
			isheader = 1;
			namestart = i + 1;
			if (i != 0) {
				fprintf(stderr, "Only one FastA sequence per file is allowed!\n");
				exit(1);
			}
		}
		if (isheader) {
			if (data[i] == '\n') {
				isheader = 0;
				if (name) {
					*name = (char *) malloc (i - namestart);
					memcpy (*name, data + namestart, i - namestart);
					(*name)[i - namestart] = 0;
				}
			}
		} else {
			if (data[i] < 'A') {
				continue;
			} else if (strchr(alphabet, data[i]) == NULL) {
				word = 0;
				m = 0;
				location += 1;
				continue;
			}
			word <<= 2;
			word |= getnuclvalue(data[i]);
			m += 1;
			if (m > wordlength) {
				word &= mask;
				m = wordlength;
			}
			if (m == wordlength) {
				table->words[table->nwords] = word;
				table->locations[table->nloc] = location + 1 - table->wordlength;
				table->nwords += 1;
				table->nloc += 1;
				/* location += 1; */
			}
			location += 1;
		}

	}
	/*printf("alguses nloc %u\n", table->nloc);*/
	return location;
}

unsigned createmask(int wordlength)
{
	int i;
	unsigned mask = 0;

	for (i = 0; i < 2 * wordlength; ++i) {
		mask = (mask << 1) | 1;
	}
	return mask;
}

void sortwords(wordtable *table)
{
	int firstshift = 0, i;
	if (table->nwords == 0) return;

	/* calculate the number of shifted positions for making radix sort faster (no need to sort digits that are all zeros)*/
	for (i = 8; i < table->wordlength * 2; i+= 8) firstshift += 8;

	hybridInPlaceRadixSort256(table->words, table->words + table->nwords, table->locations, firstshift);

	return;
}

unsigned countunique(wordtable *table)
{
	unsigned i, count;
	count = 0;
	for (i = 0; i < table->nwords; ++i) {
		if (i == 0 || table->words[i] != table->words[i - 1])
			++count;
	}
	return count;
}

void findstartpositions(wordtable *table)
{
	unsigned ri, wi, count, nunique;

	if (table->nwords == 0) return;
	nunique = countunique(table);

	if (nunique > table->nstart_slots) {
		table->nstart_slots = nunique;
		table->starts = (unsigned *) realloc(table->starts, table->nstart_slots * sizeof(unsigned));
	}

	wi = 1;
	count = 1;
	table->starts[0] = 0;
	for (ri = 1; ri < table->nwords; ++ri) {
		if (table->words[ri] == table->words[ri - 1]) {
			count += 1;
		} else {
			table->words[wi - 1] = table->words[ri - 1];
			table->starts[wi] = count;
			count += 1;
			++wi;
		}
	}
	table->words[wi - 1] = table->words[ri - 1];
	table->nwords = nunique;
	table->nstarts = nunique;
	return;
}

void sortlocations(wordtable *table)
{
	unsigned i;
	for (i = 0; i < table->nwords - 1; ++i) {
		hybridInPlaceRadixSort256(table->locations + table->starts[i], table->locations + table->starts[i + 1], NULL, 24);
	}
	hybridInPlaceRadixSort256(table->locations + table->starts[i], table->locations + table->nloc, NULL, 24);
}

void combineindices(wordtable *table, wordtable *other)
{
	long long i, j, k, l;
	int n1, n2;
	unsigned count, location_pos, count1 = 0, count2 = 0;

	if (table->nloc + other->nloc > table->nloc_slots) {
		table->nloc_slots = table->nloc + other->nloc;
		table->locations = (unsigned *) realloc(table->locations, table->nloc_slots * sizeof(unsigned));
	}
	location_pos = table->nloc + other->nloc - 1;

	j = 0;
	count = 0;
	for (i = 0; i < other->nwords; ) {
		if (j < table->nwords && table->words[j] == other->words[i]) {
			++i;
			++j;
		} else if (j < table->nwords && table->words[j] < other->words[i]) {
			++j;
		} else {
			++i;
			++count;
		}
	}

	if (table->nwords + count > table->nword_slots) {
		table->nword_slots = table->nwords + count;
		table->words = (unsigned *) realloc(table->words, table->nword_slots * sizeof(unsigned));
	}
	if (table->nstarts + count > table->nstart_slots) {
		table->nstart_slots = table->nwords + count;
		table->starts = (unsigned *) realloc(table->starts, table->nstart_slots * sizeof(unsigned));
	}

	i = other->nwords - 1;
	j = table->nwords - 1;
	n1 = table->nloc - table->starts[j];
	n2 = other->nloc - other->starts[i];
	for (k = table->nwords + count - 1; k >= 0; --k) {
		if (i >= 0 && j >= 0 && table->words[j] == other->words[i]) {
			for (l = n2 - 1; l >= 0; --l) {
				table->locations[location_pos] = other->locations[other->starts[i] + l];
				location_pos--;
				count1++;
			}
			for (l = n1 - 1; l >= 0; --l) {
				table->locations[location_pos] = table->locations[table->starts[j] + l];
				location_pos--;
				count1++;
			}
			n1 = (j > 0) ? table->starts[j] - table->starts[j - 1] : 0;
			n2 = (i > 0) ? other->starts[i] - other->starts[i - 1] : 0;
			table->words[k] = table->words[j];
			table->starts[k] = location_pos + 1;
			--j; --i;
		} else if ((i < 0) || (j >= 0 && table->words[j] > other->words[i])) {
			for (l = n1 - 1; l >= 0; --l) {
				table->locations[location_pos] = table->locations[table->starts[j] + l];
				location_pos--;
				count1++;
			}
			n1 = (j > 0) ? table->starts[j] - table->starts[j - 1] : 0;
			table->words[k] = table->words[j];
			table->starts[k] = location_pos + 1;
			--j;

		} else {
			for (l = n2 - 1; l >= 0; --l) {
				table->locations[location_pos] = other->locations[other->starts[i] + l];
				location_pos--;
				count2++;
			}
			n2 = (i > 0) ? other->starts[i] - other->starts[i - 1] : 0;
			table->words[k] = other->words[i];
			table->starts[k] = location_pos + 1;
			--i;
		}
	}
	table->nwords += count;
	table->nstarts += count;
	table->nloc += other->nloc;
}

void writetoindex(wordtable *table, const char *outputname, sequences seqinfo)
{
	unsigned long long i;
	char fname[256];
	FILE *f;
	info h;
	if (table->nwords == 0) return;

	h.wordsize = table->wordlength;
	h.nwords = table->nwords;
	h.nlocations = table->nloc;
	h.seqs = seqinfo;

	if (debug) {
		long long i;
		for (i = 0; i < seqinfo.n; ++i) {
			printf("%s %u\n", seqinfo.seq_name[i], seqinfo.start_pos[i]);
		}
	}

	printf("%s %u\n", h.seqs.seq_name[0], h.seqs.start_pos[0]);

	sprintf(fname, "%s_%d.index", outputname, table->wordlength);
	f = fopen(fname, "w");
	fwrite(&h, sizeof(info), 1, f);

	for (i = 0; i < table->nwords; ++i) {
		fwrite(&table->words[i], sizeof(table->words[i]), 1, f);
	}
	for (i = 0; i < table->nstarts; ++i) {
		fwrite(&table->starts[i], sizeof(table->starts[i]), 1, f);
	}
	for (i = 0; i < table->nloc; ++i) {
		fwrite(&table->locations[i], sizeof(table->locations[i]), 1, f);
	}
	return;
}
