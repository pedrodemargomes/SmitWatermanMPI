#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define true (1==1)
#define false (0==0)

#define MATCH 3
#define MISS -3
#define PENALTY -2

struct Direcao {
	int i;
	int j;
};

struct Maior {
	int i;
	int j;
	int val;
	int rank;
};

int numSeq1,numSeq2;
char *seq1,*seq2;
char *respSeq1,*respSeq2;
int *M;
struct Direcao *direc;
struct Maior maior;

inline void scoreMaster(int bi, int bj) {
	if(seq1[bi-1] == seq2[bj-1]) {
		M[bi*(numSeq2+1)+bj] = M[(bi-1)*(numSeq2+1)+bj-1] + MATCH;
	} else {
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

inline void scoreSlave(int tamBlocos, int my_rank,int bi, int bj) {
	if(seq1[ (bi-1) + tamBlocos*my_rank ] == seq2[ (bj-1) ]) {
		M[bi*(numSeq2+1)+bj] = M[(bi-1)*(numSeq2+1)+bj-1] + MATCH;
		direc[bi*(numSeq2+1)+bj].i = bi-1;
		direc[bi*(numSeq2+1)+bj].j = bj-1;
	} else {
		M[bi*(numSeq2+1)+bj] = M[(bi-1)*(numSeq2+1) + bj-1] + MISS;
		direc[bi*(numSeq2+1)+bj].i = bi-1;
		direc[bi*(numSeq2+1)+bj].j = bj-1;
	}
	if(M[bi*(numSeq2+1)+bj] < M[(bi-1)*(numSeq2+1) + bj] + PENALTY ) {
		M[bi*(numSeq2+1)+bj] = M[(bi-1)*(numSeq2+1) + bj] + PENALTY;
		direc[bi*(numSeq2+1)+bj].i = bi-1;
		direc[bi*(numSeq2+1)+bj].j = bj;
	}
	if(M[bi*(numSeq2+1)+bj] < M[bi*(numSeq2+1)+ bj-1] + PENALTY ) {
		M[bi*(numSeq2+1)+bj] = M[bi*(numSeq2+1)+ bj-1] + PENALTY;
		direc[bi*(numSeq2+1)+bj].i = bi;
		direc[bi*(numSeq2+1)+bj].j = bj-1;
	}
	if(M[bi*(numSeq2+1)+bj] < 0) {
		M[bi*(numSeq2+1)+bj] = 0;
	}

	if( M[bi*(numSeq2+1)+bj] > maior.val ) {
		maior.val = M[bi*(numSeq2+1)+bj];
		maior.i = bi;
		maior.j = bj;
	}

}


void complexcompare(struct Maior * in, struct Maior * inout, int * len, MPI_Datatype * dptr) {
	int i;
  	struct Maior c;
  	for (i = 0; i < *len; i++) {
  	  if(in->val > inout->val) {
  	    	/* compares the module of the complex in *in and the complex in *inout */
  	    	*inout = *in;
  	 	}
  		in++;
		inout++;
	}
}

int main(int argc, char *argv[]) {
	MPI_Op myOp;
	MPI_Datatype Maior;

	int i,j,ii,jj;
	int tamBlocos;
	int numBlocos;

	double tiG,tfG;
	double tiGet, tfGet;
	double tiTotal,tfTotal;

	int my_rank,n_procs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

	MPI_Type_contiguous(4, MPI_INT, &Maior);
	MPI_Type_commit( &Maior );

	MPI_Op_create((MPI_User_function *) complexcompare, true, &myOp);

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

	MPI_Bcast(seq1, numSeq1, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(seq2, numSeq2, MPI_CHAR, 0, MPI_COMM_WORLD);

	M = calloc( (tamBlocos+1)*(numSeq1+1), sizeof(int));
	direc = malloc( (tamBlocos+1)*(numSeq1+1)* sizeof(struct Direcao));
	
	//while(1) {}

	MPI_Request reqSend;

	maior.val = 0; maior.i = 0; maior.j = 0; maior.rank = my_rank;

	for(j=0;j<numBlocos;j++) {

		// Se nao eh processador que calcula a primeira linha de blocos
		if(my_rank != 0) {
			// Recebe a a linha de cima do processo my_rank-1
			MPI_Recv(&(M[ j*tamBlocos +1 ]), tamBlocos, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}


		// printf("++++++++++++++ rank = %d\n",my_rank);
		// for(ii=0;ii<tamBlocos+1;ii++) {
		// 	for(jj=0;jj<numSeq2+1;jj++) {
		// 		printf("%d ",M[ii*(numSeq1+1)+jj]);
		// 	}
		// 	printf("\n");
		// }

		// ++ Calcula bloco ++
		for(ii=1;ii<tamBlocos+1;ii++) {
			for(jj=j*tamBlocos+1;jj<(j+1)*tamBlocos+1;jj++) {				
				//printf("%d %d\n",ii,jj);
				//printf("%d\n",ii*(numSeq2+1)+jj);
				scoreSlave(tamBlocos,my_rank,ii,jj);
			}
		}
		// +++++++++++++++++++
		
		// Se nao eh processador que calcula a ultima linha de blocos
		if(my_rank != n_procs-1) {
			// Manda a linha de baixo do processo para my_rank+1
			MPI_Isend(&(M[ (tamBlocos)*(numSeq1+1) + j*tamBlocos +1 ]), tamBlocos, MPI_INT, my_rank+1, 0, MPI_COMM_WORLD, &reqSend );
		}

	}

	// if(my_rank == 0) {
	//  	printf("++++++++++++++ rank = %d\n",my_rank);
	// 	for(ii=0;ii<tamBlocos+1;ii++) {
	// 		for(jj=0;jj<numSeq2+1;jj++) {
	// 			printf("%d ",M[ii*(numSeq1+1)+jj]);
	// 		}
	// 		printf("\n");
	// 	}
	// }


	struct Maior maiorGlobal;
	MPI_Allreduce(&maior, &maiorGlobal, 1, Maior, myOp, MPI_COMM_WORLD);

	//printf("maior global =  %d maior i = %d maior j = %d maior rank = %d\n", maiorGlobal.val,maiorGlobal.i,maiorGlobal.j,maiorGlobal.rank);
	//while(1){}

	int count = 0;
	respSeq1 = malloc( (numSeq1*2+1) * sizeof(char));
	respSeq2 = malloc( (numSeq2*2+1) * sizeof(char));

	if( maiorGlobal.rank == my_rank ) {

		i = maiorGlobal.i; j = maiorGlobal.j;
		while(i > 0 && M[i*(numSeq2+1)+j] != 0 ) {
			if( direc[i*(numSeq2+1)+j].i == i-1 && direc[i*(numSeq2+1)+j].j == j-1 ) {
				respSeq1[count] = seq1[ (i-1) + tamBlocos*my_rank ];
				respSeq2[count] = seq2[j-1];
				i = i-1;j = j-1; // Diagonal
				// printf("Diagonal\n");
			} else if( direc[i*(numSeq2+1)+j].i == i-1 && direc[i*(numSeq2+1)+j].j == j ) {
				respSeq1[count] = seq1[ (i-1) + tamBlocos*my_rank ];
				respSeq2[count] = '-';
				i = i-1; // Cima
				// printf("Cima\n");
			} else {
				respSeq1[count] = '-';
				respSeq2[count] = seq2[j-1];
				j = j-1; // Esquerda
				// printf("Esquerda\n");
			}
			count++;
		}
		
		//printf("Manda para rank = %d\n", my_rank-1);
		
		if(M[i*(numSeq2+1)+j] == 0) {
			// Terminou 
			#ifdef DEBUG
			respSeq1[count] = '\0';
			respSeq2[count] = '\0';
			printf("%s\n",respSeq1 );
			printf("%s\n",respSeq2 );
			#endif
			if(my_rank != 0) {
				j = -1;
				MPI_Send(&j, 1, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD ); // FIM
			}
		} else {
			if(my_rank != 0) {
				MPI_Send(&j, 1, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD );
				MPI_Send(&count, 1, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD );
				MPI_Send(respSeq1, count, MPI_CHAR, my_rank-1, 0, MPI_COMM_WORLD );
				MPI_Send(respSeq2, count, MPI_CHAR, my_rank-1, 0, MPI_COMM_WORLD );
			}
		}
	} else if( my_rank < maiorGlobal.rank ) {
		MPI_Recv(&j,1,MPI_INT,my_rank+1,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if( j == -1) {
			// j = -1;
			if(my_rank != 0)
				MPI_Send( &j , 1, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD ); // FIM

			MPI_Finalize();
			return 0;

		} else {
			MPI_Recv(&count, 1, MPI_INT, my_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
			MPI_Recv(respSeq1, count, MPI_CHAR, my_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
			MPI_Recv(respSeq2, count, MPI_CHAR, my_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
		}

		i = tamBlocos;
		while(i > 0 && M[i*(numSeq2+1)+j] != 0 && j != -1) {
			if( direc[i*(numSeq2+1)+j].i == i-1 && direc[i*(numSeq2+1)+j].j == j-1 ) {
				respSeq1[count] = seq1[ (i-1) + tamBlocos*my_rank ];
				respSeq2[count] = seq2[j-1];
				i = i-1;j = j-1; // Diagonal
				// printf("Diagonal\n");
			} else if( direc[i*(numSeq2+1)+j].i == i-1 && direc[i*(numSeq2+1)+j].j == j ) {
				respSeq1[count] = seq1[ (i-1) + tamBlocos*my_rank ];
				respSeq2[count] = '-';
				i = i-1; // Cima
				// printf("Cima\n");
			} else {
				respSeq1[count] = '-';
				respSeq2[count] = seq2[j-1];
				j = j-1; // Esquerda
				// printf("Esquerda\n");
			}
			count++;
		}

		// respSeq1[count] = '\0';
		// respSeq2[count] = '\0';
		// printf("respSeq1 = %s\n",respSeq1 );
		// printf("respSeq2 = %s\n",respSeq2 );

		
		//printf("Manda para rank = %d\n", my_rank-1);
		if(M[i*(numSeq2+1)+j] == 0) {
			// Terminou 
			#ifdef DEBUG
			respSeq1[count] = '\0';
			respSeq2[count] = '\0';
			printf("%s\n",respSeq1 );
			printf("%s\n",respSeq2 );
			// while(1){}
			#endif DEBUG
			if(my_rank != 0) {
				j = -1;
				MPI_Send( &j , 1, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD ); // FIM
			}
		} else {
			if(my_rank != 0) {
				MPI_Send(&j, 1, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD );
				MPI_Send(&count, 1, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD );
				MPI_Send(respSeq1, count, MPI_CHAR, my_rank-1, 0, MPI_COMM_WORLD );
				MPI_Send(respSeq2, count, MPI_CHAR, my_rank-1, 0, MPI_COMM_WORLD );
			}
		}

	} else {
		//printf("FIM rank = %d\n", my_rank);
	}


	// #ifdef DEBUG
	// if(my_rank == 1) {
	// 	for(i=1;i<tamBlocos+1;i++) {
	// 		for(j=0;j<numSeq2+1;j++) {
	// 			printf("%d ",M[i*(numSeq1+1)+j]);
	// 		}
	// 		printf("\n");
	// 	}
	// }
	// #endif

	MPI_Finalize();

	return 0;
}



