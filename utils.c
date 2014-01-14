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

/* get the value of the reverse complement of the given word */
/* TODO: pole päris õige */
unsigned long long getreversecomplement(unsigned long long w, unsigned length)
{
	int i, mask, v;
	unsigned long long revcompl = 0L;
	w = ~w;
	mask = 3;
	for (i = 0; i < length; ++i) {
		v = w & mask;
		revcompl <<= 2;
		revcompl |= v;
		w >>= 2;
	}
	return revcompl;
}



 /* this implementation is based on:
	http://www.drdobbs.com/architecture-and-design/algorithm-improvement-through-performanc/221600153?pgno=1

	In-place MSD hybrid radix sort (with insertion sort) */

/* used for small buckets */
void insertionSort(unsigned *begin, unsigned *end)
{
	unsigned *p, *q;
	unsigned temp;
	for (p = begin + 1; p != end; ++p) {
		for (q = p; q != begin && *q < *(q - 1); --q) {
			temp = *q;
			*q = *(q - 1);
			*(q - 1) = temp;
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
		insertionSort(begin, end);
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
	for (i = 0; i < end - begin; ++i)  {
		p = begin + i;
		digit = (*p >> shift) & 255;
		if (bins[digit] <= binsize[digit]) {
			++p;
			continue;
		}
		position = positions[digit] + binsize[digit];
		binsize[digit]++;
		if (p == begin + position) {
			++p;
			continue;
		}
		temp = *p;
		temp_location = beg_location[i];
		*p = begin[position];
		beg_location[i] = beg_location[position];
		begin[position] = temp;
		beg_location[position] = temp_location;

	}

	/* recursive step */
	if (shift > 0) {
		for (i = 0; i < 256; ++i) {
			if (bins[i] > 1) {
				hybridInPlaceRadixSort256(begin + positions[i],
						begin + positions[i] + bins[i], beg_location + positions[i], shift - 8);
			}
		}
	}
	return;
}
