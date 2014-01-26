/*
 * Genome Indexer
 * (using radix sort)
 *
 * Authors: Maarja Lepamets, Fanny-Dhelia Pajuste
 */

#ifndef INDEXCREATER_H_
#define INDEXCREATER_H_

#ifndef __UTILS_CPP__
extern const char *alphabet;
#endif

/* Maximum number of mismatched regions */
#define MAX_REGIONS 4

typedef struct _wordtable {
	int wordlength;
	unsigned nword_slots;
	unsigned nwords;
	unsigned nstart_slots;
	unsigned nstarts;
	unsigned nloc_slots;
	unsigned nloc;
	unsigned *words;
	unsigned *starts;
	unsigned *locations;
} wordtable;

typedef struct _sequence {
	int n;
	char **seq_name;
	unsigned *start_pos;
} sequences;

typedef struct _info {
	int wordsize;
	unsigned nwords;
	unsigned nlocations;
	sequences seqs;
} info;

/* Mismatched region of candidate */
typedef struct _region {
	unsigned qstart;
	unsigned qend;
	unsigned loc;
} region;

typedef struct _candidate {
	unsigned loc;
	unsigned mmis;
	unsigned length;
	unsigned nregions;
	region reg[MAX_REGIONS];
} candidate;

typedef struct _queryblock {
	const char *query;
	candidate *candidates;
	unsigned ncandidates;
} queryblock;

typedef struct _chromosome {
	char *name;
	char *filename;
	unsigned start;
	unsigned length;
	char *sequence;
} Chromosome;

/*
 * functions are defined and commented in utils.c
 */
int getnuclvalue(char nucl);
unsigned int getreversecomplementstr (char *dst, const char *seq, unsigned len);

void hybridInPlaceRadixSort256(unsigned *begin, unsigned *end, unsigned *beg_location, unsigned shift);
char* word2string(unsigned w, int wordlength);

unsigned get_seeds (const char *query, unsigned *words, unsigned nwords, unsigned wordlen, unsigned m, unsigned *seeds);
unsigned search_word (unsigned word, unsigned *words, unsigned nwords);

unsigned find_candidates (queryblock *qb, unsigned *words, unsigned nwords, unsigned *starts, unsigned *locations, unsigned nlocations,
		unsigned wordlen, unsigned m, unsigned *seeds, unsigned max_candidates, unsigned mmis);


#endif /* INDEXCREATER_H_ */
