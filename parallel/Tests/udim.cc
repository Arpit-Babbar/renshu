#include <iostream>
#include <cmath>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

#include <cstdlib> // Random number

#include "../../include/array3d.h"
#include "../../include/array3d.cc"
// TODO - Use a proper make file

using namespace std;

#define p_dim 3

int main(int argc, char** argv)
{
  // Objective - make dimensional independent code. Not gonna work because
  // of loops :/
  int myid, numprocs, ierr; // rank, size renamed for problem
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  int N = 32;
  int spat_dim[p_dim]; // N x N x N grid
  int proc_dim[p_dim]; // np1 X np2 X np2
  int pbc_check[p_dim]; // 0/1 to indicate open/periodic bc in particular dimension

  // The book does this by taking input from user, which is why they do
  // it only in id 0. We should make an input file. Input from user will be
  // too slow.
  if (myid == 0)
  {
    for (int i = 0; i < p_dim; i++)
    {
      pbc_check[i] = 0; // Open bc in all directions
      spat_dim[i]  = N; // Assuming equal size in all directions
      proc_dim[i]  = 0; // Zero means that MPI_Dims_create specifies these
    }
  }
  // Broadcast respective arrays from root = 0 to all other processes
  ierr = MPI_Bcast(spat_dim,       // Starting address of buffer
                   p_dim,          // Number of entries in buffer
                   MPI_INTEGER,    // MPI_Datatype
                   0,              // Root
                   MPI_COMM_WORLD  // MPI_Comm
                   );
  ierr = MPI_Bcast(proc_dim, p_dim, MPI_INTEGER, 0, MPI_COMM_WORLD);
  ierr = MPI_Bcast(pbc_check, p_dim, MPI_INTEGER, 0, MPI_COMM_WORLD);

  double xmin = 0.0, xmax = 1.0;
  double ymin = 0.0, ymax = 1.0;
  double zmin = 0.0, zmax = 1.0;
  (void)xmax,(void)ymax,(void)zmax; // To supress warnings
  double h = 1.0/N; // assume dx = dy = dx

  ierr = MPI_Dims_create(numprocs,   // Total processes
                         p_dim,      // Physical Dimension
                         proc_dim    // array storing (np1,np2,np3)
                         );
  if (myid==0)
    for (int i = 0; i < p_dim; i++)
      printf("proc_dim[%d] = %d\n", i, proc_dim[i]);

  int let_reorder = 1; // Allow reordering to handle all numprocs, sending
                       // surplus to MPI_COMM_NULL

  MPI_Comm GRID_COMM_WORLD;

  ierr = MPI_Cart_create(MPI_COMM_WORLD,   // Old communicator
                         p_dim,            // number of dimensions
                         proc_dim,         // no. of procs in each dimension
                         pbc_check,        // periodicity indicator
                         let_reorder,      // Allow reordering(1) or not(0)
                         &GRID_COMM_WORLD  // output new communicator
                         );

  int myid_grid, num_p_grid; // Rank and size are now different in
                             // GRID_COMM_WORLD and must be respecified

  MPI_Comm_rank(GRID_COMM_WORLD, &myid_grid);
  MPI_Comm_size(GRID_COMM_WORLD, &num_p_grid);

  // The number of points in each dimension allocated to current process
  // and coordinates of current process

  int local_dim[p_dim]/* = {0}*/, mycoord[p_dim]/* = {0}*/;

  // We obtain mycoord from MPI_Cart_coords()
  ierr = MPI_Cart_coords(GRID_COMM_WORLD,   //communicator
                         myid_grid,         // my rank
                         p_dim,             // Dimension of grid
                         mycoord            // output array of cart coords
                         );

  // local_dim = physical grid points allocated to process
  for (int i = 0; i < p_dim; i++)
  {
    local_dim[i]  = (int) floor(spat_dim[i] / proc_dim[i]);
                  // Total points in dim     / # processes in dimension
    // We just have to put them in some process. Obviously, this number will
    // be less number of processes in this dimension. So, we iteratively start
    // to put these points from 0th rank in this dimension to as many as there are.
    if (mycoord[i] < spat_dim[i] % proc_dim[i])
                     // total unused points
    {
      local_dim[i] += 1; // Adding unused point
    }
    // proc_dim>unused_pts, so #(procs for which if holds)=#( unused pts)
  }
  int Ni = local_dim[2], Nj = local_dim[1], Nk = local_dim[0];
  printf("For rank %d with coordinates (%d,%d,%d), Nk, Nj, Ni = %d, %d, %d\n",
          myid, mycoord[0], mycoord[1], mycoord[2], Nk, Nj, Ni);
  Array3D phi[2];
  phi[0].resize(Ni+2,Nj+2,Nk+2); phi[1].resize(Ni+2,Nj+2,Nk+2);

  int totmsgsize[p_dim], MaxBufLen=0;

  // j-k plane
  totmsgsize[2] = local_dim[0] * local_dim[1]; // Nk*Nj
  MaxBufLen     = fmax(totmsgsize[2], MaxBufLen);
  // i-k plane
  totmsgsize[1] = local_dim[0] * local_dim[2]; // Nk*Ni
  MaxBufLen     = fmax(totmsgsize[1], MaxBufLen);
  // i-j plane
  totmsgsize[0] = local_dim[1] * local_dim[2]; // Nj*Ni
  MaxBufLen     = fmax(totmsgsize[0], MaxBufLen);

  double fieldSend[MaxBufLen] = {0.0}; // usage: phi -> fieldSend -> MPI_Send
  double fieldRecv[MaxBufLen] = {0.0}; // usage: MPI_Recv -> fieldRecv -> phi

  // left, right physical limits, i.e., values to be updated
  int udim[2][p_dim] = {0};
  int disp = -1;
  // TO BE DEBUGGED BY UNIT TESTS
  for (int dir = 0; dir < p_dim; dir++)
  {
    int source, dest;
    ierr = MPI_Cart_shift(GRID_COMM_WORLD,
                          dir,             // physical direction
                          disp,            // displacement
                          &source,         // output negative dir proc
                          &dest            // output positive dir proc
                          );
    if (dest != MPI_PROC_NULL) // non-boundary, neighbour on 'left'
    {
      udim[0][dir] = 1;
      assert(mycoord[dir]>0);
      printf("Along %d, for rank %d with coordinates (%d,%d,%d), neighbor on left\n",
             dir, myid_grid, mycoord[0],mycoord[1],mycoord[1]);
    }
    else                       // boundary, no neighbour of 'left'
    {
      udim[0][dir] = 2;
      assert(mycoord[dir]==0);
      printf("Along %d, for rank %d with coordinates (%d,%d,%d), no neighbor on left\n",
             dir, myid_grid, mycoord[0],mycoord[1],mycoord[1]);
    }
    if (source != MPI_PROC_NULL) // non-boundary, neighbour on 'right'
    {
      udim[1][dir] = local_dim[dir];
      assert(mycoord[dir]<proc_dim[dir]-1);
      printf("Along %d, for rank %d with coordinates (%d,%d,%d), neighbor on right\n",
             dir, myid_grid, mycoord[0],mycoord[1],mycoord[1]);
    }
    else                         // boundary, no neighbour on 'right'
    {
      udim[1][dir] = local_dim[dir] - 1;
      assert(mycoord[dir]==proc_dim[dir]-1);
      printf("Along %d, for rank %d with coordinates (%d,%d,%d), no neighbor on right\n",
             dir, myid_grid, mycoord[0],mycoord[1],mycoord[1]);
    }
  } // udim[][0] -> Nk, udim[][1] -> Nj, udim[][2] -> Ni
  ierr = MPI_Finalize();
}

// What about splitting as 1 x 1 x np?
