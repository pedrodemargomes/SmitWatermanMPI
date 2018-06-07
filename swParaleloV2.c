#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define MATCH 3
#define MISS -3
#define PENALTY -2

int numSeq1,numSeq2;
char *seq1,*seq2;
int *M;

inline void score(int bi, int bj) {
	if(seq1[bi-1] == seq2[bj-1]) {
		M[bi*(numSeq2+1)+bj] = M[(bi-1)*(numSeq2+1)+bj-1] + MATCH;
	}
	else {
		M[bi*(numSeq2+1)+bj] = M[(bi-1)*(numSeq2+1) + bj-1] + MISS;
	}
	if(M[bi*(numSeq2+1)+bj] < M[(bi-1)*(numSeq2+1) + bj] + PENALTY ) {
		M[bi*(numSeq2+1)+bj] = M[(bi-1)*(numSeq2+1) + bj] + PENALTY;
	}
	if(M[bi*(numSeq2+1)+bj] < M[bi*(numSeq2+1)+ bj-1] + PENALTY ) {
		M[bi*(numSeq2+1)+bj] = M[bi*(numSeq2+1)+ bj-1] + PENALTY;
	}
	if(M[bi*(numSeq2+1)+bj] < 0)
		M[bi*(numSeq2+1)+bj] = 0;
}

int main(int argc, char *argv[]) {
	int i,j,ii,jj;
	int tamBlocos;
	int numBlocos;

	int my_rank,n_procs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

	if(my_rank == 0) {
		scanf("%d %d",&numSeq1,&numSeq2);
	}

	MPI_Bcast(&numSeq1, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&numSeq2, 1, MPI_INT, 0, MPI_COMM_WORLD);

	tamBlocos = numSeq2/n_procs;
	numBlocos = n_procs;

	seq1 = malloc(numSeq1*sizeof(char));
	seq2 = malloc(numSeq2*sizeof(char));

	if(my_rank == 0) {
		scanf("%s",seq1);
		scanf("%s",seq2);
	}

	//printf("%d\n",my_rank );

	MPI_Bcast(seq1, numSeq1, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(seq2, numSeq2, MPI_CHAR, 0, MPI_COMM_WORLD);

	//while(1) {}

	//printf("%d\n",my_rank );
	M = calloc( (numSeq1+1) * (numSeq2+1), sizeof(int));
	for(i =0;i< (numSeq1+1) * (numSeq2+1);i++)
		M[i] = 5;
	
	while(1) {}

	MPI_Request req;

	for(j=0;j<numBlocos;j++) {

		// Se nao eh processador que calcula a primeira linha de blocos
		if(my_rank != 0) {
			// Recebe a a linha de cima do processo my_rank-1
			MPI_Recv(&(M[ (my_rank*tamBlocos)*(numSeq1+1) + j*tamBlocos +1 ]), tamBlocos, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		// ++ Calcula bloco ++
		for(ii=my_rank*tamBlocos+1;ii<(my_rank+1)*tamBlocos+1;ii++) {
			for(jj=j*tamBlocos+1;jj<(j+1)*tamBlocos+1;jj++) {				
				score(ii,jj);
			}
		}
		// +++++++++++++++++++
		
		// Se nao eh processador que calcula a ultima linha de blocos
		if(my_rank != n_procs-1) {
			// Manda a linha de baixo do processo para my_rank+1
			MPI_Isend(&(M[ ((my_rank+1)*tamBlocos)*(numSeq1+1) + j*tamBlocos +1 ]), tamBlocos, MPI_INT, my_rank+1, 0, MPI_COMM_WORLD, &req );
		}

	}

	if(my_rank == 0) {
		for(i = 1;i<n_procs;i++) {
			MPI_Recv(&(M[ i*(tamBlocos)*(numSeq1+1) + (numSeq1+1) ]), (numSeq1+1)*tamBlocos, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
		}
	} else {
		MPI_Isend(&(M[ my_rank*(tamBlocos)*(numSeq1+1) + (numSeq1+1) ]), (numSeq1+1)*tamBlocos, MPI_INT, 0, 0, MPI_COMM_WORLD, &req );
	}


	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	#ifdef DEBUG
	if(my_rank == 0) {
		printf("rank = %d\n",my_rank);
		for(i=0;i<numSeq1+1;i++) {
			for(j=0;j<numSeq2+1;j++) {
				printf("%d ",M[i*(numSeq1+1)+j]);
			}
			printf("\n");
		}
	}
	#endif

	//while(1) {}

	return 0;
}



