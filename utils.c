/*
 * Genome Indexer
 * (using radix sort)
 *
 * Authors: Maarja Lepamets, Fanny-Dhelia Pajuste
 */

#define __UTILS_CPP__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const char *alphabet = "ACGTUacgtu";

/* get the bit value of the nucleotide (two bits) */
int getnuclvalue(char nucl)
{
	static int bit1 = 1 << 2;
	static int bit2 = 3 << 1;

	if (nucl & bit1) {
		return ((nucl >> 4) | 2) & 3;
	}
	return (nucl & bit2) >> 1;
}

/* Get reverse complement of sequence given as string
 * Returns 1 if all nucleotides are OK, 0 if there are unrecognized symbols
 * Does not 0-terminate destination string
 */
unsigned int getreversecomplementstr (char *dst, const char *seq, unsigned len)
{
	unsigned i, valid;
	static char *rev = NULL;
	if (!rev) {
		rev = (char *) malloc (256);
		for (i = 0; i < 256; i++) rev[i] = 'N';
		rev['a'] = rev['A'] = 'T';
		rev['c'] = rev['C'] = 'G';
		rev['g'] = rev['G'] = 'C';
		rev['t'] = rev['T'] = 'A';
	}
	valid = 1;
	for (i = 0; i < len; i++) {
		dst[i] = rev[(unsigned char) seq[len - 1 - i]];
		if (dst[i] == 'N') valid = 0;
	}
	return valid;
}

char* word2string(unsigned w, int wordlength)
{
	char *sequence = (char *)malloc(wordlength + 1);
	int i, temp;

	for (i = 0; i < wordlength; ++i) {
		temp = w & 3;
		sequence[wordlength - i - 1] = alphabet[temp];
		w >>= 2;
	}
	sequence[wordlength] = 0;
	return sequence;
}

/*
 * for debugging
 */
void word2bits(unsigned a)
{
	unsigned mask = (unsigned ) 1 << 31;
	while (mask != 0) {
		if ((mask & a) != 0)
			putchar('1');
		else putchar('0');
		mask = mask >> 1;
	}
	printf("\n");
}



 /* this implementation is based on:
	http://www.drdobbs.com/architecture-and-design/algorithm-improvement-through-performanc/221600153?pgno=1

	In-place MSD hybrid radix sort (with insertion sort) */

/* used for small buckets */
void insertionSort(unsigned *begin, unsigned *end, unsigned *beg_location)
{
	unsigned *p, *q;
	unsigned temp, temp_loc;
	for (p = begin + 1; p != end; ++p) {
		for (q = p; q != begin && *q < *(q - 1); --q) {
			temp = *q;
			*q = *(q - 1);
			*(q - 1) = temp;
			if (beg_location) {
				temp_loc = beg_location[q - begin];
				beg_location[q - begin] = beg_location[q - begin - 1];
				beg_location[q - begin - 1] = temp_loc;
			}
		}
	}
}

void hybridInPlaceRadixSort256(unsigned *begin, unsigned *end, unsigned *beg_location, unsigned shift)
{
	unsigned *p;
	unsigned temp, temp_location, position, digit;
	size_t bins[256];	/* for counts */
	size_t positions[256]; /* for starting positions */
	size_t binsize[256]; /* for counting the filled positions */
	int i;

	if (end - begin <= 32) {
		insertionSort(begin, end, beg_location);
		return;
	}

	memset(bins, 0, sizeof(bins));
	memset(binsize, 0, sizeof(binsize));

	/* calculating counts for every bin */
	for (p = begin; p != end; ++p) {
		digit = (*p >> shift) & 255;
		bins[digit]++;
	}

	/* finding starting positions for each bin */
	positions[0] = 0;
	for (i = 1; i < 256; ++i) {
		positions[i] = positions[i - 1] + bins[i - 1];
	}

	/* swapping words and locations */
	for (i = 0; i < end - begin; )  {
		p = begin + i;
		digit = (*p >> shift) & 255;

		if (bins[digit] <= binsize[digit]) {
			++i;
			continue;
		}
		position = positions[digit] + binsize[digit];
		binsize[digit]++;
		if (p == begin + position) {
			++i;
			continue;
		}
		temp = *p;
		*p = begin[position];
		begin[position] = temp;

		if (beg_location) {
			temp_location = beg_location[i];
			beg_location[i] = beg_location[position];
			beg_location[position] = temp_location;
		}

	}

	/* recursive step */
	if (shift > 0) {
		for (i = 0; i < 256; ++i) {
			if (bins[i] > 1) {
				if (beg_location) {
					hybridInPlaceRadixSort256(begin + positions[i],
							begin + positions[i] + bins[i], beg_location + positions[i], shift - 8);
				} else {
					hybridInPlaceRadixSort256(begin + positions[i],
							begin + positions[i] + bins[i], NULL, shift - 8);
				}
			}
		}
	}
	return;
}
