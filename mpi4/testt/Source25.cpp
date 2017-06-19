#include <stdio.h>
#include "mpi.h"
#include <cmath>
#include <stdlib.h>

#define N 4
int main(int argc, char* argv[]) {
	int rank, size, global_result;
	int vect1[N];
	int vect2[N];

	for (int i = 0; i < N; i++) {
		vect1[i] = i;
		vect2[i] = i;
	}
	vect1[0] = -1;
	vect2[0] = 56;
	/*----------------------------------------------*/
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	printf("rank : %i - ", rank);

	if (rank == 0) {
		int start = 0;
		int end;

		if (size == 1) {
			end = N;
		}
		else {
			end = N / size;
		}

		int *localVect1 = new int[end - start];
		int *localVect2 = new int[end - start];

		for (int i = 0; i < end; i++) {
			localVect1[i] = vect1[i];
			localVect2[i] = vect2[i];
		}
		//////////////////
		int result = 0;
		for (int i = 0; i < end; i++) {
			result += localVect1[i] * localVect2[i];
		}

		printf("Result %i \n", result);
		//////////////////
		for (int i = 1; i < size; i++) {
			start = (N / size) * i;
			if (i == size - 1) {
				end = N;
			}
			else {
				end = (N / size) * (i + 1);
			}

			int *buf = new int[2 * (end - start)];

			int iter = 0;
			for (int j = start; j < end; j++) {
				buf[iter] = vect1[j];
				iter++;
			}

			for (int j = start; j < end; j++) {
				buf[iter] = vect2[j];
				iter++;
			}
			//printf("%i ---- %i", iter, (2 * (end - start)));
			MPI_Send(buf, 2 * (end - start), MPI_INT, i, 10, MPI_COMM_WORLD);
		}
		MPI_Reduce(&result, &global_result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		printf("Final result : %i \n", global_result);
	}
	if (rank != 0) {
		int bufElems;

		MPI_Probe(MPI_ANY_SOURCE, 10, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_INT, &bufElems);

		int *buffer = new int[bufElems];
		MPI_Recv(buffer, bufElems, MPI_INT, 0, 10, MPI_COMM_WORLD, &status);
		/////////////////////////
		int local_result = 0;
		for (int i = 0; i < bufElems / 2; i++) {
			local_result += buffer[i] * buffer[i + bufElems / 2];
		}
		printf("Locale result : %i \n", local_result);
		/////////////////////////
		MPI_Reduce(&local_result, &global_result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		//MPI_Recv(buf, 256, MPI_INT, 0, 10, MPI_COMM_WORLD, &status);
	}

	MPI_Finalize();
	return 0;	
}

