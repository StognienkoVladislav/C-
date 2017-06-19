#include <stdio.h> 
#include <omp.h> 
#include <iostream>

using namespace std;

#define SIZE 7

int multiplyVectors() {
	int vec1[SIZE];
	int vec2[SIZE];
	int result = 0;

#pragma omp parallel for
	for (int i = 0; i < SIZE; i++) {
		vec1[i] = rand()%(i+1) + 3;
		vec2[i] = rand()%(i+1) + 2;
	}

	for (int i = 0; i < SIZE; i++) {
		cout << vec1[i] << " ";
	}
	cout << "\n";
	for (int i = 0; i < SIZE; i++) {
		cout << vec2[i] << " ";
	}

#pragma omp parallel for
	for (int i = 0; i < SIZE; i++) {
		result += vec1[i] * vec2[i];
	}

	return result;
}


int main(int argc, char **argv)
{
	cout << "\nFIRST TASK: " << multiplyVectors() << endl;

	system("pause");
	return 0;
}
