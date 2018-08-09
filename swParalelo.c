#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <time.h>

#define MATCH 3
#define MISS -3
#define PENALTY -2

#define MAX 20000
#define MAX_MATRIX MAX+1

struct Pacote {
	int ki;
	int kj;
};

struct Maior { 
	int val; 
	int i;
	int j;
};    


int numSeq1,numSeq2;
char *seq1,*seq2;
int *M; //[MAX_MATRIX][MAX_MATRIX];

int tamBloco;
struct Maior maior;
struct Pacote pkg;

int s_block;
int pi,pj; // Primeiros i e j da diagonal

int numDiagonais;
int numElementos;

int testa;

inline int max(int a, int b) {
	if(a > b)
		return a;
	return b;
}

inline int min(int a, int b) {
	if(a > b)
		return b;
	return a;
}

inline int numElementosDiagonal(int i, int numSeq1, int numSeq2) {
	if( i < numSeq1+1 && i < numSeq2+1 )
		return i;
	else if(i < max(numSeq1+1,numSeq2+1) )
		return min(numSeq1,numSeq2);
	else  
		return 2*min(numSeq1+1,numSeq2+1) - i + abs(numSeq1+1 -(numSeq2+1) ) - 2;
}

inline void calcPrimElemDiagonal(int i,int *pi,int *pj,int numSeq1) {
	if (i < numSeq1+1) {
		*pi = i;
 		*pj = 1;
	} else {
		*pi = (numSeq1+1) - 1;
		*pj = i -(numSeq1+1) + 2;
	}
}

inline void score(int bi, int bj){
	if(seq1[bi-1] == seq2[bj-1]) {
		M[bi*(numSeq2+1)+bj] = M[(bi-1)*(numSeq2+1)+bj-1] + MATCH;
	}
	else {
		M[bi*(numSeq2+1) + bj] = M[(bi-1)*(numSeq2+1) + bj-1] + MISS;
	}
	if(M[bi*(numSeq2+1)+ bj] < M[(bi-1)*(numSeq2+1) + bj] + PENALTY ) {
		M[bi*(numSeq2+1)+ bj] = M[(bi-1)*(numSeq2+1) + bj] + PENALTY;
	}
	if(M[bi*(numSeq2+1)+ bj] < M[bi*(numSeq2+1)+ bj-1] + PENALTY ) {
		M[bi*(numSeq2+1)+ bj] = M[bi*(numSeq2+1)+ bj-1] + PENALTY;
	}
	if(M[bi*(numSeq2+1)+ bj] < 0)
		M[bi*(numSeq2+1)+ bj] = 0;
}

int main(int argc, char *argv[]) {
	int i,j,ki,kj;
	int my_rank,n_procs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

	tamBloco= 4;//omp_get_num_procs(); numero de blocos por linha

	if(my_rank == 0) {

		scanf("%d %d",&numSeq1,&numSeq2);
		
		seq1 = malloc(numSeq1*sizeof(char));
		seq2 = malloc(numSeq2*sizeof(char));

		scanf("%s",seq1);
		scanf("%s",seq2);

		M = calloc( (numSeq1+1) * (numSeq2+1), sizeof(int));

		s_block = ceil((double)numSeq2/tamBloco);

		
		numDiagonais = numSeq1/s_block+numSeq2/s_block-1; // Numero de diagonais a percorrer
	}

	MPI_Bcast(&s_block, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&numSeq1, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&numSeq2, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Datatype col_matrix;
	MPI_Type_vector(numSeq1+1, 1, numSeq2+1, MPI_INT, &col_matrix);
	MPI_Type_commit(&col_matrix);

	MPI_Status status;

	if(my_rank != 0) {
		M = calloc( (numSeq1+1) * (numSeq2+1) * 2, sizeof(int));
		seq1 = malloc(numSeq1*sizeof(char));
		seq2 = malloc(numSeq2*sizeof(char));
	}

	MPI_Bcast(seq1, numSeq1, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(seq2, numSeq2, MPI_CHAR, 0, MPI_COMM_WORLD);

	MPI_Datatype mysubarray;
	MPI_Type_vector(s_block, s_block, numSeq2+1, MPI_INT, &mysubarray );
	//MPI_Type_create_subarray(2, array_size, array_subsize, array_start, MPI_ORDER_C, MPI_INT, &mysubarray);
	MPI_Type_commit(&mysubarray);



	if(my_rank == 0) {
		int slaves[100];
		printf("numDiagonais = %d\n", numDiagonais);
		for(i =1; i <= numDiagonais; i++) {
			numElementos = numElementosDiagonal(i,numSeq1/s_block,numSeq2/s_block);
			calcPrimElemDiagonal(i,&pi,&pj,numSeq1/s_block);

			printf("numElementos = %d\n", numElementos);
			for(j = 1;j <= numElementos; j++) {
				pkg.ki = pi - j +1;
				pkg.kj = pj + j -1;

				MPI_Recv(NULL, 0, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);

				slaves[j] = status.MPI_SOURCE;

				printf("ki = %d kj = %d\n", pkg.ki,pkg.kj);
				MPI_Send(&pkg, 2, MPI_INT, status.MPI_SOURCE, 2, MPI_COMM_WORLD);
				MPI_Send(&(M[((pkg.ki-1)*s_block)*(numSeq2+1) + (pkg.kj-1)*s_block+1]), s_block, MPI_INT, status.MPI_SOURCE, 3, MPI_COMM_WORLD); // Manda linha de cima do bloco
				MPI_Send(&(M[((pkg.ki-1)*s_block+1)*(numSeq2+1) + (pkg.kj-1)*s_block]), 1, col_matrix, status.MPI_SOURCE, 4, MPI_COMM_WORLD); // Manda coluna da esquerda do bloco
			}
			for(j = 1;j <= numElementos; j++) {
				pkg.ki = pi - j +1;
				pkg.kj = pj + j -1;

				MPI_Recv(&(M[((pkg.ki-1)*s_block+1)*(numSeq2+1) + (pkg.kj-1)*s_block+1]), 1, mysubarray, slaves[j], 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}

			// TESTE
			int ii,jj;
			for (ii = 0; ii < numSeq1+1; ++ii) {
				for (jj = 0; jj < numSeq2+1; ++jj) {
					printf("%d ", M[ii*(numSeq2+1) + jj]);
				}
				printf("\n");
			}
			// TESTE <<

		}
		pkg.ki = -1;
		for(i= 1;i<n_procs; i++)
			MPI_Send(&pkg, 2, MPI_INT, i, 0, MPI_COMM_WORLD);
	} else {
		int bi,bj; // i e j do bloco
		for(;;) {
			MPI_Ssend(NULL, 0, MPI_INT, 0, 1, MPI_COMM_WORLD);

			MPI_Recv(&pkg,2,MPI_INT,0,2,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if(pkg.ki == -1)
				break;

			MPI_Recv(&(M[((pkg.ki-1)*s_block)*(numSeq2+1) + (pkg.kj-1)*s_block+1]), s_block, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Recebe linha de cima do bloco
			MPI_Recv(&(M[((pkg.ki-1)*s_block+1)*(numSeq2+1) + (pkg.kj-1)*s_block]), 1, col_matrix, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Recebe coluna da esquerda do bloco
			// ++++++++ BLOCO +++++++++
			for(bi=(pkg.ki-1)*s_block+1; bi < (pkg.ki-1)*s_block+1+s_block; bi++) {
				for(bj=(pkg.kj-1)*s_block+1;bj < (pkg.kj-1)*s_block+1+s_block; bj++) {
					score(bi,bj);
				}
			}
			// +++++++++++++++++++++++++
			MPI_Send(&(M[((pkg.ki-1)*s_block+1)*(numSeq2+1) + (pkg.kj-1)*s_block+1]), 1, mysubarray, 0, 5, MPI_COMM_WORLD);
		}
	}

	MPI_Finalize();

	#ifdef DEBUG
	if(my_rank == 0) {
		printf("numSeq1 = %d\nnumSeq2 = %d\n my_rank = %d\n",numSeq1,numSeq2,my_rank);
		for(i=0;i<numSeq1+1;i++) {
			for(j=0;j<numSeq2+1;j++) {
				printf("%d ",M[i*(numSeq2+1)+j]);
			}
			printf("\n");
		}
	}
	#endif

	return 0;
}



