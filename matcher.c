/*
 * matcher.c
 *
 *  Created on: 14 Jan 2014
 *      Author: Fanny-Dhelia
 */

#include <stdio.h>
#include <string.h>

int min(int a, int b, int c) {
	if (a <= b && a <= c) {
		return a;
	}
	if (b <= a && b <= c) {
		return b;
	} else {
		return c;
	}
}


int editDistance1(char* sg, char* sr) {
	unsigned mismatches = strlen(sg)-strlen(sr);
	unsigned arr1[strlen(sr) + 1];
	unsigned arr2[strlen(sr) + 1];
	unsigned temp;
	unsigned i;
	unsigned j;
	unsigned k;
	unsigned first = 0;
	for (i = 0; i < strlen(sg); i++) {
		arr1[i] = i;
	}
	for (j = 0; j < strlen(sg); j++) {
		first++;
		if (mismatches > 0) {
			first = 0;
			mismatches--;
		}
		arr2[0] = first;
		for (i = 0; i < strlen(sr); i++) {
			if (sg[j] == sr[i]) {
				arr2[i + 1] = arr1[i];
			} else {

				arr2[i + 1] = min(arr1[i], arr1[i + 1], arr2[i]) + 1;
			}
		}
		for (k = 0; k < strlen(sr) + 1; k++) {
			temp = arr1[k];
			arr1[k] = arr2[k];
			arr2[k] = temp;
		}
		arr2[0] = first;
		for (k = 0; k < strlen(sr) + 1; k++) {
			printf("%d ", arr1[k]);
		}
		printf("\n");
	}
	return arr1[strlen(sr)];
}

int editDistance2(char* sg, char* sr) {
	unsigned arr1[strlen(sr) + 1];
	unsigned arr2[strlen(sr) + 1];
	unsigned temp;
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
		for (k = 0; k < strlen(sr) + 1; k++) {
			temp = arr1[k];
			arr1[k] = arr2[k];
			arr2[k] = temp;
		}
		arr2[0] = j+2;
		for (k = 0; k < strlen(sr) + 1; k++) {
			printf("%d ", arr1[k]);
		}
		printf("\n");
	}
	return arr1[strlen(sr)];
}

int editDistance3(char* sg, char* sr) {
	unsigned mismatches = strlen(sg)-strlen(sr);
	unsigned arr1[strlen(sr) + 1];
	unsigned arr2[strlen(sr) + 1];
	unsigned temp;
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
		for (k = 0; k < strlen(sr) + 1; k++) {
			temp = arr1[k];
			arr1[k] = arr2[k];
			arr2[k] = temp;
		}
		arr2[0] = j+2;
		for (k = 0; k < strlen(sr) + 1; k++) {
			printf("%d ", arr1[k]);
		}
		printf("\n");
		if (j >= strlen(sr)){
			if (arr1[strlen(sr)] < best){
				best = arr1[strlen(sr)];
				last = j;
			}
		}
	}
	return best;
}

int main() {
	unsigned ed = editDistance1("ATGGA", "ATGA");
	printf("Edit Distance: %d \n", ed);
	unsigned ed2 = editDistance2("ATGGA", "TGGC");
	printf("Edit Distance: %d \n", ed2);
	unsigned ed3 = editDistance3("ATGCGA", "TGCG");
	printf("Edit Distance: %d \n", ed3);
	return 0;
}
