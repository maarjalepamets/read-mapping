/*
 * Genome Indexer
 * (using radix sort)
 *
 * Text Algorithms, 2013/2014
 *
 * Author: Maarja Lepamets
 */

#ifndef INDEXCREATER_H_
#define INDEXCREATER_H_

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


int getnuclvalue(char nucl);
void hybridInPlaceRadixSort256(unsigned *begin, unsigned *end, unsigned *beg_location, unsigned shift);




#endif /* INDEXCREATER_H_ */
