#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include <mympi.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define BLOCK_LOW(id, p, n) ((id)*(n)/(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id) +1, p, n) -1)
#define BLOCK_SIZE(id, p, n) (BLOCK_HIGH(id, p, n) - BLOCK_LOW(id, p, n) + 1)
#define BLOCK_OWNER(j, p, n) (((p)*((j)+1)-1)/(n))

int main(int argc, char* argv[]){
	int count;		// local prime count
	double elapsed_time;	// parallel execution time
	int first;		// index of first multiple
	int global_count;	// global prime count
	int high_value;		// highest value on this proc
	int i;
	int id;			// process id number
	int index;		// index of current prime 
	int low_value;		// lowest value on this proc
	char* marked;		// portion of 2, ... 'n'
	int n;			// sievimg from 2, ... 'n'
	int p;			// number of processes
	int proc0_size;		// size of proc 0's subarray
	int prime;		// current prime
	int size;		// elements in 'marked'
	int *pinakas;
	MPI_Status st;
	
	
	MPI_Init(&argc, &argv[]);
	
	//start the timer
	MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time = -MPI_Wtime();
	
	//get id ant o
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	//check for command line arguments
	if(argc != 2){
		if(!id)
			printf( "Usage: %s <m>\n", argv[0]);
		MPI_Finalize();
		exit(1);
	}
	n = atoi(argv[1]);
	
	if(id == 0){
		if((pinakas = malloc(sizeof(int) * p)) == NULL){
			printf("Cannot allocate enough memory\n");
			MPI_Finalize();
			exit(1);
		}
	}
	
	//get low and high values
	low_value = 2+BLOCK_LOW(id, p, n-1);
	high_value = 2+BLOCK_HIGH(id, p, n-1);
	size = BLOCK_SIZE(id, p, n-1);
	
	//check all primes used for sieving are contained in block 0
	proc0_size=(n-1)/p;
	if((2+proc0_size)<(int)sqrt((double)n)){
		if(!id)
			printf("Too many processes\n");
		MPI_Finalize();
		exit(1);
	}
	
	//Ok. allocate process's space
	marked = (char*) malloc(size);
	if(marked == NULL){
		printf("Cannot allocate enough memory\n");
		MPI_Finalize();
		exit(1);
	}
	
	//Initialize
	for(i=0; i<size; i++)
		marked[i] = 0;
	if(!id)
		index = 0;
	prime = 2;
	
	//start sieving
	do{
		if(prime*prime > low_value)
			first = prime*prime - low_value;
		else{
			if(!(low_value % prime))
				first = 0;
			else
				first = prime - (low_value % prime);
		}
		
		for(i=first; i<size; i+= prime)
			marked[i] = 1;
		if(!id){
			while(marked[++index]);
			prime = index+2;
		}
		MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}while (prime*prime <= n);
	
	//collect the results
	count = 0;
	for(i=0; i<size; i++){
		if(!marked[i])
			count++;
	}

	MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	if(id == 0){
		pinakas[0] = count;
		for(i=1; i<p; i++){
			MPI_Recv(&(pinakas[i]), 1, MPI_INT, i, 1, MPI_COMM_WORLD, &st);
		}
	}
	else{
		MPI_Send(&count, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
	}
	
	//stop the timer
	elapsed_time += MPI_Wtime();
	
	//print the results
	if(!id){
		printf("process\t|primes\t|pososto 1\t|pososto 2\n");
		printf("--------+-------+---------------+----------\n");
		for(i = 0; i< p; i++){
			printf("%d\t|%d\t|%f\t|%f\n", i, pinakas[i], ((float)pinakas[i])/((float)size), ((float)pinakas[i])/((float) global_count));
		}
		printf("%d primes are less than or equal to %d.\n", global_count, n);
		printf("Total_elapsed_time: %10.6f.\n", elapsed_time);
	}
	MPI_Finalize();
	return 0;
}