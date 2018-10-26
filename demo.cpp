#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]){
  int rank, size;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int number;
  if (rank == 0){
    number = -1;
    MPI_Send(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    MPI_Send(&number, 1, MPI_INT, 2, 0, MPI_COMM_WORLD);

    printf("Sent number: %d\n", number);
  } else {
    MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("Process: %d Received number: %d\n", rank, number);
  }
  MPI_Finalize();
  return 0;
}
