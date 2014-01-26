/*
 * Genome Mapper
 *
 * Authors: Maarja Lepamets, Fanny-Dhelia Pajuste
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

#include "utils.h"

#define MAX_READ_LENGTH 1000
#define MAX_CANDIDATES 10000

int debug = 0;

static int editDistance(const char *query, const char *sequence, unsigned *qstart, char *s, char *q);
void printindex(const char *data);
const char* filemmap(const char *filename, struct stat *st);
static void mapperwrapper(const char *queryfile, const char *index, info *h, int mmis, int d, Chromosome *chr, unsigned nchr);
/* Return edit distance */
static unsigned adjustmapping (unsigned queryidx, candidate *cand, const char *query, unsigned qlen, Chromosome *chr, unsigned int nchr, unsigned int nmm, unsigned *loc, char *s, char *q, unsigned reverse);

int main (int argc, const char *argv[])
{
	int i;
	int mmis = 0;
	int step = 5;
	const char *indexfile = NULL, *queryfile = NULL, *namefile = NULL;
	struct stat stindex;
	const char *ind;
	info *h;
	Chromosome chr[256];
	int nchr = 0;
	const char *chrmap;

	for (i = 1; i < argc; ++i) {
		if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "--debug")) {
			debug += 1;
		} else if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "--input")) {
			if (!argv[i + 1] || argv[i + 1][0] == '-') {
				fprintf(stderr, "Error: No index file specified!\n");
				exit(1);
			}
			indexfile = argv[i + 1];
			++i;
		} else if (!strcmp(argv[i], "-q") || !strcmp(argv[i], "--query")) {
			if (!argv[i + 1] || argv[i + 1][0] == '-') {
				fprintf(stderr, "Error: No query file specified!\n");
				exit(1);
			}
			queryfile = argv[i + 1];
			++i;
		} else if (!strcmp(argv[i], "-g") || !strcmp(argv[i], "--genome")) {
			if (!argv[i + 1] || argv[i + 1][0] == '-') {
				fprintf(stderr, "Error: No genome (.names file) specified!\n");
				exit(1);
			}
			namefile = argv[i + 1];
			++i;
		} else if (!strcmp(argv[i], "-mm") || !strcmp(argv[i], "--mismatches")) {
			if (!argv[i + 1]) {
				fprintf(stderr, "Warning: No number of mismatches specified! Using the default value: %d.\n", mmis);
				break;
			}
			char *e;
			mmis = strtol (argv[i + 1], &e, 10);
			if (*e != 0) {
				fprintf(stderr, "Invalid input: %s! Must be an integer.\n", argv[i + 1]);
				exit(1);
			}
			++i;
		} else if (!strcmp(argv[i], "-step")) {
			if (!argv[i + 1]) {
				fprintf(stderr, "Warning: No step length specified! Using the default value: %d.\n", step);
				break;
			}
			char *e;
			step = strtol (argv[i + 1], &e, 10);
			if (*e != 0) {
				fprintf(stderr, "Invalid input: %s! Must be an integer.\n", argv[i + 1]);
				exit(1);
			}
			++i;
		} else {
			fprintf(stderr, "Unknown argument: %s\n", argv[i]);
			exit(1);
		}
	}

	/* checking parameters */
	if (!indexfile || !queryfile || !namefile) {
		fprintf(stderr, "Error: Some of the input files are missing!\n");
		exit(1);
	}
	if (mmis < 0 || mmis > 10) {
		fprintf(stderr, "Error: Number of mismatches must be between 0 and 10!\n");
		exit(1);
	}
	if (step < 1 || step > 10) {
		fprintf(stderr, "Error: Step length must be between 1 and 10!\n");
		exit(1);
	}

	/* Parse chromosome description file */
	chrmap = filemmap (namefile, &stindex);
	if (chrmap) {
		unsigned s = 0;
		if (debug > 0) fprintf (stderr, "Chromosome locations:\n");
		while (s < stindex.st_size) {
			unsigned ns, ne, fs, fe, ps, pe;
			ns = ne = s;
			while (chrmap[ne] > ' ') ne += 1;
			fs = ne;
			while ((fs < stindex.st_size) && (chrmap[fs] <= ' ')) fs += 1;
			fe = fs;
			while (chrmap[fe] > ' ') fe += 1;
			ps = fe;
			while ((ps < stindex.st_size) && (chrmap[ps] <= ' ')) ps += 1;
			pe = ps;
			while (chrmap[pe] > ' ') pe += 1;
			if ((ne > ns) && (fe > fs) && (pe > ps)) {
				chr[nchr].name = (char *) malloc (ne - ns + 1);
				memcpy (chr[nchr].name, chrmap + ns, ne - ns + 1);
				chr[nchr].name[ne - ns] = 0;
				chr[nchr].filename = (char *) malloc (fe - fs + 1);
				memcpy (chr[nchr].filename, chrmap + fs, fe - fs + 1);
				chr[nchr].filename[fe - fs] = 0;
				chr[nchr].start = strtol(chrmap + ps, NULL, 10);
				chr[nchr].length = 0;
				chr[nchr].sequence = NULL;
				if (debug > 0) fprintf (stderr, "Chromosome %s at %s position %u\n", chr[nchr].name, chr[nchr].filename, chr[nchr].start);
				nchr += 1;
			}
			s = pe;
			while ((s < stindex.st_size) && (chrmap[s] <= ' ')) s += 1;
		}
	}

	ind = filemmap(indexfile, &stindex);
	h = (info *) ind;

	mapperwrapper(queryfile, ind, h, mmis, step, chr, nchr);

	return 0;
}

static void mapperwrapper(const char *queryfile, const char *index, info *h, int mmis, int d, Chromosome *chr, unsigned nchr)
{
	unsigned int j;
	unsigned *words, *starts, *locations;
	unsigned nwords, nloc;
	unsigned ncandidates;
	int wordlength;
	unsigned queryidx;

	char readfw[MAX_READ_LENGTH];
	FILE *q;

	unsigned *seeds = NULL;
	unsigned nseed_slots = 0;
	queryblock qb;
	qb.candidates = (candidate *) malloc (MAX_CANDIDATES * sizeof(candidate));

	q = fopen(queryfile, "r");
	if (q == NULL) {
		fprintf(stderr, "mapperwrapper: Cannot open file %s.\n", queryfile);
		exit (1);
	}

	wordlength = h->wordsize;
	nwords = h->nwords;
	nloc = h->nlocations;

	words = (unsigned *)(index + sizeof(info));
	starts = (unsigned *)(index + sizeof(info) + nwords * sizeof(unsigned));
	locations = (unsigned *)(index + sizeof(info) + nwords * sizeof(unsigned) + nwords * sizeof(unsigned));

	queryidx = 0;
	while (fscanf(q, "%s\n", readfw) != EOF) {
		unsigned int len = (unsigned int)strlen(readfw);
		if (nseed_slots < len) {
			nseed_slots = len;
			seeds = (unsigned *) realloc (seeds, nseed_slots * sizeof(unsigned));
		}
		if (len > 0) {
			char r[256];
			unsigned nmatched;
			if (debug > 1) fprintf (stderr, "Query: %s\n", readfw);

			qb.query = readfw;
			ncandidates = find_candidates (&qb, words, nwords, starts, locations, nloc, wordlength, d, seeds, MAX_CANDIDATES, mmis);
			if (debug > 0) {
				fprintf (stderr, "Found %u candidates:\n", ncandidates);
				if (debug > 1) {
					for (j = 0; j < ncandidates; j++) {
						unsigned k;
						fprintf (stderr, "Candidate %u location %u length %u nregions %u\n", j, qb.candidates[j].loc, qb.candidates[j].length, qb.candidates[j].nregions);
						if (debug > 2) {
							for (k = 0; k < qb.candidates[j].nregions; ++k) {
								fprintf (stderr, "    Region %u qstart %u qend %u loc %u\n", k, qb.candidates[j].reg[k].qstart, qb.candidates[j].reg[k].qend, qb.candidates[j].reg[k].loc);
							}
						}
					}
				}
			}
			nmatched = 0;
			for (j = 0; j < ncandidates; j++) {
				char s[256], q[256];
				unsigned qstart = 0, editdist;
				editdist = adjustmapping (queryidx, &qb.candidates[j], qb.query, len, chr, nchr, mmis, &qstart, s, q, 0);
				if (editdist <= mmis) nmatched += 1;
			}
			/* Reverse complement */
			getreversecomplementstr (r, readfw, len);
			r[len] = 0;
			memset(readfw, 0, sizeof(readfw));
			memset(qb.candidates, 0, MAX_CANDIDATES * sizeof(candidate));
			qb.query = r;
			if (debug > 1) fprintf (stderr, "Reverse Query: %s\n", qb.query);
			ncandidates = find_candidates (&qb, words, nwords, starts, locations, nloc, wordlength, d, seeds, MAX_CANDIDATES, mmis);
			if (debug > 2) fprintf(stderr, "kandidaatide arv: %u, neist esimene: %d, mismatche %d\n", ncandidates, qb.candidates[0].loc, qb.candidates[0].mmis);
			if (debug > 1) {
				fprintf (stderr, "Candidates:\n");
				for (j = 0; j < ncandidates; j++) {
					unsigned k;
					fprintf (stderr, "Loc %u len %u nregions %u\n", qb.candidates[j].loc, qb.candidates[j].length, qb.candidates[j].nregions);
					for (k = 0; k < qb.candidates[j].nregions; ++k) {
						fprintf (stderr, "  Region %u qstart %u qend %u loc %u\n", k, qb.candidates[j].reg[k].qstart, qb.candidates[j].reg[k].qend, qb.candidates[j].reg[k].loc);
					}
				}
			}
			for (j = 0; j < ncandidates; j++) {
				char s[256], q[256];
				unsigned qstart = 0, editdist;
				editdist = adjustmapping (queryidx, &qb.candidates[j], qb.query, len, chr, nchr, mmis, &qstart, s, q, 1);
				if (editdist <= mmis) nmatched += 1;
			}
			if (!nmatched) fprintf (stderr, "%d\t-\n", queryidx);
			memset(readfw, 0, sizeof(readfw));
			memset(qb.candidates, 0, MAX_CANDIDATES * sizeof(candidate));
		}
		queryidx += 1;
	}

}

static unsigned adjustmapping (unsigned queryidx, candidate *cand, const char *query, unsigned qlen, Chromosome *chr, unsigned int nchr, unsigned int nmm, unsigned *qstart, char *s, char *q, unsigned reverse)
{
	unsigned i, editdist;
	unsigned sloc, slen;
	char seq[256];
	for (i = 1; i < nchr; i++) {
		if (cand->loc < chr[i].start) break;
	}
	i -= 1;
	/* Load chromosome if not already loaded */
	if (chr[i].sequence == NULL) {
		struct stat st;
		unsigned d, s, header;
		const char *seq = filemmap (chr[i].filename, &st);
		if (seq == NULL) {
			fprintf (stderr, "Cannot mmap %s\n", chr[i].filename);
			exit (1);
		}
		chr[i].sequence = (char *) malloc (st.st_size);
		d = 0;
		header = 0;
		for (s = 0; s < st.st_size; s++) {
			if (seq[s] == '>') header = 1;
			if (header) {
				if (seq[s] == '\n') header = 0;
			} else {
				if (seq[s] >= 'A') chr[i].sequence[d++] = seq[s];
			}
		}
		chr[i].length = d;
		munmap ((void *) seq, st.st_size);
	}

	sloc = cand->loc - nmm;
	slen = qlen + 2 * nmm;
	memcpy (seq, chr[i].sequence + sloc - chr[i].start, slen);
	seq[slen] = 0;
	/* editdist = editDistanceMiddle (q, s); */
	*qstart = nmm;
	editdist = editDistance (query, seq, qstart, s, q);
	if (debug > 0) {
		fprintf (stderr, "Location %u Distance %u Query start %u\n", cand->loc, editdist, *qstart);
		fprintf (stderr, "Query: %s\n", q);
		fprintf (stderr, "Seq:   %s\n", s);
	}
	if (editdist <= nmm) {
		fprintf (stderr, "%u\t%s\t%u\t%d\t%s\n", queryidx, chr[i].name, cand->loc -nmm + *qstart - chr[i].start, editdist, (reverse) ? "R" : "F");
	}
	return editdist;
}

/* fixme: Update qstart */

static int editDistance (const char *query, const char *seq, unsigned *qstart, char *s, char *q)
{
	static int *d = NULL;
	static int dsize = 0;
	int qlen, slen, qi, si, dist, lasts;
	qlen = strlen (query);
	slen = strlen (seq);
	if (debug > 2) fprintf (stderr, "editDistance: Query %s (len = %d), sequence %s (len = %d)\n", query, qlen, seq, slen);
	if ((qlen + 1) * (slen + 1) > dsize) {
		dsize = (qlen + 1) * (slen + 1);
		d = (int *) realloc (d, dsize * sizeof (int));
	}
	/* Fill first column */
	for (si = 0; si <= slen; si++) {
		/*d[si * (qlen + 1) + 0] = (si < (int) *qstart) ? *qstart - si : si - *qstart;*/
		d[si * (qlen + 1) + 0] = (si < (int) 2 * *qstart) ? 0 : si - 2 * *qstart;
	}
	/* Fill first row */
	for (qi = 1; qi <= qlen; qi++) {
		d[qi] = 99;
	}
	for (qi = 1; qi <= qlen; qi++) {
		for (si = 1; si <= slen; si++) {
			unsigned dl, dtl, dt;
			dl = d[si * (qlen + 1) + qi - 1] + 1;
			dtl = d[(si - 1) * (qlen + 1) + qi - 1];
			if ((qi < qlen) && (si < slen) && (query[qi] != seq[si])) dtl += 1;
			dt = d[(si - 1) * (qlen + 1) + qi] + 1;
			d[si * (qlen + 1) + qi] = ((dl < dtl) && (dl < dt)) ? dl : ((dt < dl) && (dt < dtl)) ? dt : dtl;
		}
	}
	if (debug > 3) {
		/* Print table */
		for (si = 0; si <= slen; si++) {
			for (qi = 0; qi <= qlen; qi++) {
				fprintf (stderr, "%2u ", d[si * (qlen + 1) + qi]);
			}
			fprintf (stderr, "\n");
		}
	}
	dist = 99;
	lasts = 0;
	for (si = 0; si <= slen; si++) {
		if (d[si * (qlen + 1) + qlen] < dist) {
			dist = d[si * (qlen + 1) + qlen];
			lasts = si;
		}
	}
	if (s && q) {
		unsigned sp, qp, i;
		/* Print alignment */
		sp = 0;
		qp = 0;
		qi = qlen - 1;
		si = slen - 1;
		while (si >= lasts) {
			s[sp++] = seq[si--];
			q[qp++] = '-';
		}
		while ((qi >= 0) || (si >= 0)) {
			unsigned dl, dtl, dt;
			dl = (qi > 0) ? d[si * (qlen + 1) + qi - 1] : 99;
			dtl = ((qi > 0) && (si > 0)) ? d[(si - 1) * (qlen + 1) + qi - 1] : 99;
			dt = (si > 0) ? d[(si - 1) * (qlen + 1) + qi] : 99;
			if (qi < 0) {
				/* End of query */
				s[sp++] = seq[si--];
				q[qp++] = '-';
			} else if (si < 0) {
				/* End of sequence */
				s[sp++] = '-';
				q[qp++] = query[qi--];
			} else if (qi == 0) {
				s[sp++] = seq[si--];
				q[qp++] = query[qi--];
			} else if ((dtl <= dl) && (dtl <= dt)) {
				/* Match or ungapped mismatch */
				s[sp++] = seq[si--];
				q[qp++] = query[qi--];
			} else if (si == 0) {
				s[sp++] = seq[si--];
				q[qp++] = '-';
			} else if ((dl <= dtl) && (dl <= dt)) {
				/* Gap in sequence */
				s[sp++] = '-';
				q[qp++] = query[qi--];
			} else {
				/* Gap in query */
				s[sp++] = seq[si--];
				q[qp++] = '-';
			}
		}
		for (i = 0; i < qp / 2; i++) {
			char t = q[i];
			q[i] = q[qp - 1 - i];
			q[qp - 1 - i] = t;
		}
		for (i = 0; i < sp / 2; i++) {
			char t = s[i];
			s[i] = s[sp - 1 - i];
			s[sp - 1 - i] = t;
		}
		s[sp] = 0;
		q[qp] = 0;
	}
	return dist;
}

void printindex(const char *data)
{
	unsigned *words, *starts, *locations;
	unsigned wordlength, nwords, nlocations;
	info *h;
	unsigned i, j;

	h = (info *)data;
	wordlength = h->wordsize;
	nwords = h->nwords;
	nlocations = h->nlocations;

	words = (unsigned *)(data + sizeof(info));
	starts = (unsigned *)(data + sizeof(info) + nwords * sizeof(unsigned));
	locations = (unsigned *)(data + sizeof(info) + nwords * sizeof(unsigned) + nwords * sizeof(unsigned));

	for (i = 0; i < nwords - 1; ++i) {
		fprintf(stdout, "%s\t%u\n", word2string(words[i], wordlength), starts[i]);
		for (j = starts[i]; j < starts[i + 1]; ++j) {
			fprintf(stdout, "%u ", locations[j]);
		}
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "%s\t%u\n", word2string(words[i], wordlength), starts[i]);
	for (j = starts[i]; j < nlocations; ++j) {
		fprintf(stdout, "%u ", locations[j]);
	}
	fprintf(stdout, "\n");

}

const char* filemmap(const char *filename, struct stat *st)
{
	int status, handle;
	const char *data;

	/* memory-mapping a file */
	status = stat(filename, st);
	if (status < 0) {
		fprintf (stderr, "Cannot get the statistics of file %s!\n", filename);
		exit (1);
	}
	handle = open(filename, O_RDONLY);
	if (handle < 0) {
		fprintf (stderr, "Cannot open file %s!\n", filename);
		exit (1);
	}
	data = (const char *) mmap(NULL, st->st_size, PROT_READ, MAP_SHARED, handle, 0);
	if (data == (const char *) -1) {
		fprintf (stderr, "Cannot memory-map file %s!\n", filename);
		exit (1);
	}
	close(handle);

	return data;
}


