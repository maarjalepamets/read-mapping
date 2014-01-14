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

#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

#include "utils.h"

const char *alphabet = "ACGTUacgtu";

/* reading in the genome in FastA format (also Multi-FastA) */
unsigned readfastafile(const char *filename, wordtable *table, unsigned loc);

/* fills the table with words and their locations */
unsigned fillwordtable(const char *data, wordtable *table, struct stat st, unsigned loc);

/* mask has as many 1's as the wordlength is */
unsigned createmask(int wordlength);

/* wrapper for in-place radix sort */
void sortwords(wordtable *table);

/* */
unsigned countunique(wordtable *table);

/* */
void findstartpositions(wordtable *table);

/* */
void combineindices(wordtable *table, wordtable *other);

int main (int argc, const char *argv[])
{
	int i, inputbeg, inputend, location = 0;
	const char *outputname;
	wordtable tables[2];
	wordtable *table = &tables[0];
	wordtable *temptable = &tables[1];
	memset(table, 0, sizeof(wordtable));
	memset(temptable, 0, sizeof(wordtable));

	/* default value */
	int n = 16;

	/* parsing commandline arguments */
	for (i = 1; i < argc; ++i) {
		if (!strcmp(argv[i], "-i")) {
			/* get the locations of the input files */
			inputbeg = ++i;
			while (i < argc - 1 && argv[i + 1][0] != '-') {
				++i;
			}
			inputend = i;
		} else if (!strcmp(argv[i], "-o")) {
			outputname = argv[i + 1];
			++i;
		} else if (!strcmp(argv[i], "-n")) {
			n = atoi(argv[i + 1]);
			++i;
		} else {
			fprintf(stderr, "Unknown argument: %s\n", argv[i]);
			exit(1);
		}
	}

	/* checking the parameter */
	if (n > 16) {
		fprintf(stderr, "Seed size too large!\n");
		exit(1);
	}

	/* work-flow */
	for (i = inputbeg; i <= inputend; ++i) { /* iteration over the input file parameters */
		temptable->wordlength = n;
		if (i > inputbeg) {
			location += 10000;
		}
		location = readfastafile(argv[i], temptable, location);
		if (temptable->nwords == 0) continue;
		sortwords(temptable);
		findstartpositions(temptable);
		if (i > inputbeg) {
			combineindices(table, temptable);
		} else {
			*table = *temptable;
		}
		memset(temptable, 0, sizeof(wordtable));
	}

	/*writetoindex(table, outputname);*/

	return 0;
}

unsigned readfastafile(const char *filename, wordtable *table, unsigned locprev)
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
	loc = fillwordtable(data, table, st, locprev);

	munmap((void *) data, st.st_size);
	return loc;
}

unsigned fillwordtable(const char *data, wordtable *table, struct stat st, unsigned locprev)
{
	unsigned word, location;
	unsigned mask;
	int m, wordlength;
	int isheader = 0;
	off_t i;

	wordlength = table->wordlength;
	mask = createmask(wordlength);
	word = 0;
	location = locprev;
	m = 0;

	/* find the unsigned integer corresponding to a word with a given length */
	for (i = 0; i < st.st_size; ++i) {
		if (data[i] == '>') {
			isheader = 1;
			m = 0;
			word = 0;
			location += 10000;
		}
		if (isheader) {
			if (data[i] == '\n')
				isheader = 0;
		} else {
			if (strchr(alphabet, data[i]) == NULL && data[i] != '\n') {
				word = 0;
				m = 0;
				location += 10000;
				continue;
			} else if (data[i] == '\n') {
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
				table->locations[table->nloc] = location;
				table->nwords += 1;
				table->nloc += 1;
				location += 1;
			}
		}

	}
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
	int firstshift = 0;
	if (table->nwords == 0) return;

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
	wi = 0;
	count = 1;

	if (table->nwords == 0) return;

	nunique = countunique(table);

	if (nunique > table->nstart_slots) {
		table->nstart_slots = nunique;
		table->starts = (unsigned *) realloc(table->starts, table->nstart_slots * sizeof(unsigned));
	}

	for (ri = 1; ri < table->nwords; ++ri) {
		if (table->words[ri] == table->words[ri - 1]) {
			++count;
		} else {
			table->words[wi] = table->words[ri - 1];
			table->starts[wi] = count;
			count = 1;
			++wi;
		}
	}
	table->words[wi] = table->words[ri - 1];
	table->starts[wi] = count;
	table->nwords = nunique;
	return;
}

/* for combining two sorted lists of words (with locations) */
void combineindices(wordtable *table, wordtable *other)
{
	long long i, j, k, count_new;

	count_new = 0;
	j = 0;

	if (table->nloc + other->nloc > table->nloc_slots) {
		table->nloc_slots = table->nloc + other->nloc;
		table->locations = (unsigned *) realloc(table->locations, table->nloc_slots * sizeof(unsigned));
	}

	for (i = 0; (unsigned long long) i < other->nwords; ) {
		if (table->words[j] == other->words[i]) {
			++i; ++j;
		} else if (table->words[j] < other->words[i]) {
			++j;
		} else {
			++i; ++count_new;
		}
	}

	if (count_new > 0) {

		if (table->nwords + count_new > table->nword_slots) {
			table->nword_slots = table->nwords + count_new;
			table->words = (unsigned *) realloc(table->words, table->nword_slots * sizeof(unsigned));
		}
		if (table->nstarts + count_new > table->nstart_slots) {
			table->nstart_slots = table->nwords + count_new;
			table->starts = (unsigned *) realloc(table->starts, table->nstart_slots * sizeof(unsigned));
		}

		i = other->nwords - 1;
		j = table->nwords - 1;
		for (k = table->nwords + count - 1; k >= 0; --k) {
			if ((i < 0) || (j >= 0 && table->words[j] > other->words[i])) {
				table->words[k] = table->words[j];
				table->frequencies[k] = table->frequencies[j];
				--j;
			} else {
				table->words[k] = other->words[i];
				table->frequencies[k] = other->frequencies[i];
				--i;
			}
		}
		table->nwords += count;
	}
}
