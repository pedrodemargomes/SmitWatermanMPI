#include <stdio.h>
#include <stdlib.h>
#include <time.h>
int main (){
	srand(time(NULL));
	int x, y, r;
	scanf("%d %d", &x, &y);
	for (int i=0;i <x; i++){
		r=rand()%4;
		if (r==0){printf("T");}
		if (r==1){printf("A");}
		if (r==2){printf("C");}
		if (r==3){printf("G");}
	}
	printf("\n");	
	for (int i=0; i<y; i++){
		r=rand()%4;
		if (r==0){printf("T");}
		if (r==1){printf("A");}
		if (r==2){printf("C");}
		if (r==3){printf("G");}
	}
	printf("\n");	
	return 0;
}