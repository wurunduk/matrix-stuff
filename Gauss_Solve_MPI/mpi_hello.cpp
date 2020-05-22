/*
 * compile:
 *   mpicxx mpi_hello.cpp
 * run:
 *   mpirun -n 4 ./a.out
 */
#include "mpi.h"
#include <stdio.h>
#include <string.h>

#define BUF_LEN 256

int main (int argc, char *argv[])
{
	int my_rank;
	int p;
	int source;
	int dest;
	int tag = 0;
	char message[BUF_LEN];
	MPI_Status status;

	//Initializes the MPI execution environment 
	MPI_Init (&argc, &argv);

	//Determines the rank of the calling process in the communicator. 
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

	//Returns the size of the group associated with a communicator. 
	MPI_Comm_size (MPI_COMM_WORLD, &p);

	sprintf (message, "Hello from process %d!", my_rank);

	if (my_rank != 0)
	{
		dest = 0;
		//Performs a standard-mode blocking send. 
		MPI_Send (message, strlen (message) + 1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
	}
	else
	{
		printf ("%s\n", message);
		for (source = 1; source < p; source++)
		{
			MPI_Recv (message, BUF_LEN, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
			printf ("%s\n", message);
		}
	}

	//Terminates MPI execution environment. 
	MPI_Finalize ();
	return 0;
}

