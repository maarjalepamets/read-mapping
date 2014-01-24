/*
 * matcher.c
 *
 *  Created on: 14 Jan 2014
 *      Author: Fanny-Dhelia
 */

#include <stdio.h>
#include <string.h>

#define debug 0

unsigned min(unsigned a, unsigned b, unsigned c) {
	if (a <= b && a <= c) {
		return a;
	}
	if (b <= a && b <= c) {
		return b;
	} else {
		return c;
	}
}

void swap(unsigned* arr1, unsigned* arr2, unsigned len){
	unsigned temp;
	unsigned k;
	for (k = 0; k < len + 1; k++) {
		temp = arr1[k];
		arr1[k] = arr2[k];
		arr2[k] = temp;
	}
}

void reverse(char* s, unsigned len){
	unsigned i;
	char temp;
	for (i = 0; i < len/2; i++) {
		temp = s[i];
		s[i] = s[len-1-i];
		s[len-1-i] = temp;
	}
}


int editDistanceMiddle(char* sg, char* sr) {
	unsigned arr1[strlen(sr) + 1];
	unsigned arr2[strlen(sr) + 1];
	unsigned i;
	unsigned j;
	unsigned k;
	for (i = 0; i < strlen(sg); i++) {
		arr1[i] = i;
	}
	arr2[0] = 1;
	for (j = 0; j < strlen(sg); j++) {
		for (i = 0; i < strlen(sr); i++) {
			if (sg[j] == sr[i]) {
				arr2[i + 1] = arr1[i];
			} else {

				arr2[i + 1] = min(arr1[i], arr1[i + 1], arr2[i]) + 1;
			}
		}
		swap(arr1, arr2, strlen(sr));
		arr2[0] = j+2;
		if (debug > 0) {
			for (k = 0; k < strlen(sr) + 1; k++) {
				printf("%d ", arr1[k]);
			}
			printf("\n");
		}
	}
	return arr1[strlen(sr)];
}

int editDistanceEnd(char* sg, char* sr, unsigned *loc) {
	unsigned arr1[strlen(sr) + 1];
	unsigned arr2[strlen(sr) + 1];
	unsigned i;
	unsigned j;
	unsigned k;
	unsigned best = strlen(sg)+1;
	unsigned last;
	for (i = 0; i < strlen(sg); i++) {
		arr1[i] = i;
	}
	arr2[0] = 1;
	for (j = 0; j < strlen(sg); j++) {
		for (i = 0; i < strlen(sr); i++) {
			if (sg[j] == sr[i]) {
				arr2[i + 1] = arr1[i];
			} else {

				arr2[i + 1] = min(arr1[i], arr1[i + 1], arr2[i]) + 1;
			}
		}
		swap(arr1, arr2, strlen(sr));
		arr2[0] = j+2;
		for (k = 0; k < strlen(sr) + 1; k++) {
			printf("%d ", arr1[k]);
		}
		printf("\n");
		if (j >= strlen(sr)-1){
			if (arr1[strlen(sr)] < best){
				best = arr1[strlen(sr)];
				last = j;
			}
		}
	}
	*loc = last+1-strlen(sr);
	return best;
}

int editDistanceBeginning(char* sg, char* sr, unsigned *loc){
	char s1[strlen(sg)];
	char s2[strlen(sr)];
	unsigned k;
	for (k = 0; k < strlen(sg) + 1; k++) {
		s1[k] = sg[k];
	}
	for (k = 0; k < strlen(sr) + 1; k++) {
		s2[k] = sr[k];
	}
	reverse(s1, strlen(sg));
	reverse(s2, strlen(sr));
	return editDistanceEnd(s1, s2, loc);
}

#if 0
int main() {
	unsigned loc;
	unsigned ed = editDistanceBeginning("ATGGA", "ATGA", &loc);
	printf("Edit Distance: %d, loc: %d\n", ed, loc);
	ed = editDistanceMiddle("ATGGA", "TGGC");
	printf("Edit Distance: %d \n", ed);
	ed = editDistanceEnd("ATGCGA", "GCGA", &loc);
	printf("Edit Distance: %d, loc: %d\n", ed, loc);
	return 0;
}
#endif
