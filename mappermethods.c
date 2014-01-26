/*
 * Genome Mapper
 *
 * Authors: Maarja Lepamets, Fanny-Dhelia Pajuste
 */

#include <stdio.h>
#include <string.h>
#include <malloc.h>

#include "utils.h"

#define debug 0

typedef unsigned u32;

/*
 * Per seed location function
 *
 * n          - index of current seed
 * locations  - list of all locations in genome
 * pos        - per seed array of indices into locations
 * m          - step between seeds
 *
 * returns    - the match location of given seed (from locations array)
 */

u32 loc (u32 n, u32 *locations, u32 *pos, u32 m) {
	return locations[pos[n]] - n * m;
}

/* Find all candidate locations of given query
 *
 * query      - query sequence (read)
 * words      - pointer to sorted array of all words in genome
 * nwords     - number of words (and starts)
 * starts     - starting indices of given word
 * locations  - list of all locations in genome
 * nlocations - number of distinct locations
 * wordlen    - word length
 * m          - step between seeds
 * seeds      - array where seeds will be written (has to be big enough to fit all)
 * candidates - array where candidate locations will be written
 * max_candidates - the size of candidate array
 * cutoff     - minimum ratio of confirming seeds to consider given position as candidate
 *
 * returns    - number of candidate locations
 */

u32 find_candidates (queryblock *qb, u32 *words, u32 nwords, u32 *starts, u32 *locations, u32 nlocations, u32 wordlen, u32 m, u32 *seeds, u32 max_candidates, u32 mmis) {
	/* Per seed arrays */
	/* pos is index into locations array we are currently processing */
	static u32 *pos = NULL;
	/* end is the end index (one past last) of givevn seed locations */
	static u32 *end = NULL;
	static u32 possize = 0;
	u32 nseeds, i, minloc, ncandidates;
	int minn;
	int cutoff;

	nseeds = get_seeds (qb->query, words, nwords, wordlen, m, seeds);
	cutoff = (wordlen % m == 0) ? nseeds - (wordlen / m) * mmis : nseeds - (wordlen / m + 1) * mmis;
	if (cutoff <= 0) cutoff = 1;
	if (debug) fprintf(stderr, "Siide: %u, cutoff: %d\n", nseeds, cutoff);

	if (nseeds == 0) {
		if (debug) fprintf (stderr, "Query %s gave 0 seeds\n", qb->query);
		return 0;
	}

	if (debug > 1) {
		fprintf (stderr, "Query %s gave %d seeds\n", qb->query, nseeds);
		if (debug > 2) {
			for (i = 0; i < nseeds; i++) {
				if (seeds[i] == nwords) continue;
				u32 j;
				fprintf (stderr, "Seed %d index %u word %u sequence ", i, seeds[i], words[seeds[i]]);
				for (j = 0; j < wordlen; j++) {
					static const char *n = "ACGT";
					fprintf (stderr, "%c", n[(words[seeds[i]] >> (wordlen - 1 - j)) & 3]);
				}
				fprintf (stderr, "\n");
			}
		}
	}

	/* Ensure we have enough room for pos and count arrays */
	if (nseeds > possize) {
		possize = nseeds;
		pos = (u32 *) realloc (pos, possize * sizeof (u32));
		end = (u32 *) realloc (end, possize * sizeof (u32));
	}

	/* Initialize per seed arrays */
	for (i = 0; i < nseeds; i++) {
		if (seeds[i] == nwords) continue;
		pos[i] = starts[seeds[i]];
		if (seeds[i] < (nwords - 1)) {
			end[i] = starts[seeds[i] + 1];
		} else {
			end[i] = nlocations;
		}
	}

	/* Main iteration */
	ncandidates = 0;
	while (ncandidates < max_candidates) {
		int nfound;
		/* Find minimum n */
		minn = -1;
		minloc = 0xffffffff;
		for (i = 0; i < nseeds; i++) {
			if (seeds[i] == nwords) continue;
			if (pos[i] < end[i]) {
				u32 l = loc (i, locations, pos, m);
				if ((minloc == 0xffffffff) || (l < minloc)) {
					minn = i;
					minloc = l;
				}
			}
		}
		if (minn < 0) {
			/* No more locations for any seed */
			break;
		}
		if (debug > 1) fprintf (stderr, "Found match at %u\n", minloc);
		/* Get all seeds that confirm minloc */
		/* We increase pos values here as well because we do not need them until next iteration */
		candidate cand;
		cand.loc = minloc;
		cand.mmis = 0;
		cand.length = (unsigned) strlen(qb->query);
		cand.nregions = 1;
		cand.reg[0].loc = minloc;
		cand.reg[0].qstart = 0;
		cand.reg[0].qend = cand.length;
		nfound = 0;
		if (debug > 1) fprintf (stderr, "Seeds: ");
		for (i = 0; i < nseeds; i++) {
			if (seeds[i] == nwords) continue;
			if (pos[i] < end[i]) {
				u32 l = loc (i, locations, pos, m);
				int delta = (long long) l - (long long) minloc;
				if ((delta >= -(int) mmis) && (delta <= (int) mmis)) {
					int sloc;
					/* This seed confirms given location */
					nfound += 1;
					/* Update query region list */
					sloc = m * i;
					if ((sloc > cand.reg[cand.nregions - 1].qstart) && ((sloc + wordlen) < cand.reg[cand.nregions - 1].qend)) {
						/* Split region */
						cand.reg[cand.nregions - 1].qend = sloc;
						if (cand.nregions < MAX_REGIONS) {
							cand.nregions += 1;
							cand.reg[cand.nregions - 1].loc = minloc + sloc + wordlen;
							cand.reg[cand.nregions - 1].qstart = sloc + wordlen;
							cand.reg[cand.nregions - 1].qend = cand.length;
						}
					} else if (sloc > cand.reg[cand.nregions - 1].qstart) {
						/* Clip region end */
						cand.reg[cand.nregions - 1].qend = sloc;
					} else if ((sloc + wordlen) < cand.reg[cand.nregions - 1].qend) {
						/* Clip region start */
						cand.reg[cand.nregions - 1].qstart = sloc + wordlen;
						cand.reg[cand.nregions - 1].loc = minloc + sloc + wordlen;
					}
					/* Advance pos value */
					pos[i] += 1;
					if (debug > 1) fprintf (stderr, "%u ", i);
				}
			}
		}
		if (debug > 0) fprintf (stderr, "\n");
		if ((cand.reg[cand.nregions - 1].qstart >= cand.length) || (cand.reg[cand.nregions - 1].qend <= 0)) {
			/* Region became empty */
			cand.nregions -= 1;
		}
		if (nfound >= cutoff) {
			qb->candidates[ncandidates++] = cand;
			if (debug > 1) {
				fprintf (stderr, "Found candidate location %u\n", minloc);
			}
		}
	}

	return ncandidates;
}

/*
 * Get list of seed indices (in words array)
 *
 * query      - search query
 * words      - pointer to sorted array of all words in genome
 * nwords     - number of words (and starts)
 * wordlen    - word length
 * m          - step between seeds
 * seeds      - array where seeds will be written (has to be big enough to fit all)
 *
 * returns    - number of seeds
 */

u32 get_seeds (const char *query, u32 *words, u32 nwords, u32 wordlen, u32 m, u32 *seeds) {
	static int *nucl = NULL;
	u32 qlen, pos, idx;

	/* Initialize nucleotide lookup table */
	if (nucl == NULL) {
		int i;
		nucl = (int *) malloc (256 * sizeof (int));
		for (i = 0; i < 256; i++) nucl[i] = -1;
		nucl['a'] = nucl['A'] = 0;
		nucl['c'] = nucl['C'] = 1;
		nucl['g'] = nucl['G'] = 2;
		nucl['t'] = nucl['T'] = 3;
	}

	qlen = strlen (query);
	pos = 0;
	idx = 0;
	while (pos < (qlen - wordlen)) {
		u32 word = 0;
		u32 i, index;
		for (i = 0; i < wordlen; i++) {
			if (nucl[(unsigned char) query[pos + i]] < 0) break;
			word <<= 2;
			word |= nucl[(unsigned char) query[pos + i]];
		}
		/* If we did not complete full iteration there was invalid nucleotide */
		if (i == wordlen) {
		  /* Find index of given word */
		  index = search_word (word, words, nwords);
		} else {
		  index = nwords;
		}
		seeds[idx++] = index;
		pos += m;
	}

	return idx;
}

/*
 * Binary search
 *
 * word       - current word
 * words      - pointer to sorted array of all words in genome
 * nwords     - number of words (and starts)
 *
 * returns    - the index of current word or nwords if not found
 */

u32 search_word (u32 word, u32 *words, u32 nwords) {
	u32 s, e;
	/* Do binary search */
	s = 0;
	e = nwords - 1;
	while ((e - s) > 1) {
		int m = (s + e) / 2;
		if (words[m] == word) {
			return m;
		} else if (words[m] < word) {
			s = m;
		} else {
			e = m;
		}
	}
	return nwords;
}
