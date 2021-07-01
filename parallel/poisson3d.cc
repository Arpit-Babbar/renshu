#include <iostream>
#include <cmath>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

#include "../include/array3d.h"
#include "../include/array3d.cc"
// TODO - Use a proper make file

using namespace std;

// Define physical dimension globally
#define p_dim 3
// Solve 3D Poisson's equation with source term f(x) = 1 using Jacobi's method

// Visualizing the steady state equation as time dependent equation where we
// are evolving from t0 to t1, we shall store solution_old in solution(:,t0)
// and solution_new in (:,t1) while running the iterates of Jacobi's method.

// Please see the tutorials on mpi_cart before this

double f(double x, double y, double z) // Source term
{
  double value = 1.0;
  return value;
}

// Copy what I am sending to my neighbour to MaxBufLen
void CopySendBuf(Array3D phi[],                // Solution old and new
                 int t0,
                 double disp, double dir, double fieldSend[], int MaxBufLen);

void CopyRecvBuf(Array3D phi[],                   // Solution old and new
                 int t0,
                 double disp, double dir, double fieldRecv[], int MaxBufLen);

void Jacobi_sweep(int udim[][p_dim], // local_dim(pts in dimension)
                  Array3D phi[], int t0, int t1,
                  double xmin, double ymin, double zmin,
                  double h,
                  double *maxdelta);
// Are Ni, Nj, Nk wrong ranges for some reason? Sir's code has udim.

int main(int argc, char** argv)
{
  // Objective - make dimensional independent code. Not gonna work because
  // of loops :/
  int myid, numprocs, ierr; // rank, size renamed for problem
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  int pbc_check[p_dim]; // 0/1 to indicate open/periodic bc in particular dimension
  int spat_dim[p_dim], proc_dim[p_dim];

  // The book does this by taking input from user, which is why they do
  // it only in id 0. We should make an input file. Input from user will be
  // too slow.
  int N = 60; // WHY DOES IT ONLY WORK FOR N LESS THAN 80???
  if (myid == 0)
  {
    for (int i = 0; i < p_dim; i++)
    {
      pbc_check[i] = 0; // Open bc in all directions
      spat_dim[i]  = N; // Assuming equal size in all directions
      proc_dim[i]  = 0; // Zero bc means that MPI_Dims_create specifies these
    }
  }
  double xmin = 0.0, xmax = 1.0;
  double ymin = 0.0, ymax = 1.0;
  double zmin = 0.0, zmax = 1.0;
  (void)xmax,(void)ymax,(void)zmax; // To supress warnings
  double h = 1.0/N; // assume dx = dy = dx

  // Broadcast respective arrays from root = 0 to all other processors
  ierr = MPI_Bcast(spat_dim,       // Starting address of buffer
                   p_dim,          // Number of entries in buffer
                   MPI_INTEGER,    // MPI_Datatype
                   0,              // Root
                   MPI_COMM_WORLD  // MPI_Comm
                   );
  ierr = MPI_Bcast(proc_dim, p_dim, MPI_INTEGER, 0, MPI_COMM_WORLD);
  ierr = MPI_Bcast(pbc_check, p_dim, MPI_INTEGER, 0, MPI_COMM_WORLD);

  ierr = MPI_Dims_create(numprocs,
                         p_dim,
                         proc_dim);
  if (myid==0)
    for (int i = 0; i < p_dim; i++)
      printf("proc_dim[%d] = %d\n", i, proc_dim[i]);

  int let_reorder = 1; // Allow reordering to handle all numprocs, sending
                       // surplus to MPI_COMM_NULL

  MPI_Comm GRID_COMM_WORLD;

  ierr = MPI_Cart_create(MPI_COMM_WORLD,
                         3,
                         proc_dim,
                         pbc_check,
                         let_reorder,
                         &GRID_COMM_WORLD
                         );

  int myid_grid, num_p_grid; // Rank and size are now different in
                             // GRID_COMM_WORLD and must be respecified

  MPI_Comm_rank(GRID_COMM_WORLD, &myid_grid);
  MPI_Comm_size(GRID_COMM_WORLD, &num_p_grid);

  if (GRID_COMM_WORLD == MPI_COMM_NULL)
  {
    printf("All processes are null, choose a different numprocs.\n");
    // Is this correct? It's not doing it's job, code runs through even for
    // np = 7. TODO - Fix with assertions.
  }

  // The number of points in each dimension allocated to current process
  // and coordinates of current process

  int local_dim[p_dim], mycoord[p_dim];

  // We obtain mycoord from MPI_Cart_coords()
  ierr = MPI_Cart_coords(GRID_COMM_WORLD,   //communicator
                         myid_grid,         // my rank
                         p_dim,             // Dimension of grid
                         mycoord            // Array to store coordinates in
                         );

  for (int i = 0; i < p_dim; i++)
  {
    local_dim[i]  = (int) floor (spat_dim[i] / proc_dim[i]);
                  // Total points in dim    / points allocated to proc in dim
    // spat_dim[i] % proc_dim[i] gives the total unused points
    // We just have to put them in some process. Obviously, this number will
    // be less number of processes in this dimension. So, we iteratively start
    // to put these points from 0th rank in this dimension to as many as there are.
    if (mycoord[i] < spat_dim[i] % proc_dim[i])
    {
      local_dim[i] += 1; // Adding unused points
    }
  }

  int Ni = local_dim[2], Nj = local_dim[1], Nk = local_dim[0];
  printf("For rank %d, Ni, Nj, Nk = %d, %d, %d\n", myid, Ni, Nj, Nk);
  // Don't know why the ordering has been changed!!
  Array3D phi[2];
  phi[0].resize(Ni+1,Nj+1,Nk+1); phi[1].resize(Ni+1,Nj+1,Nk+1);
  int totmsgsize[p_dim], MaxBufLen=0;

  // Might have developed a bug copying from Fortran's indexing from 1
  // j-k plane
  totmsgsize[2] = local_dim[0] * local_dim[1];
  MaxBufLen     = fmax(totmsgsize[2], MaxBufLen);
  // i-k plane
  totmsgsize[1] = local_dim[0] * local_dim[2];
  MaxBufLen     = fmax(totmsgsize[1], MaxBufLen);
  // i-j plane
  totmsgsize[0] = local_dim[1] * local_dim[2];
  MaxBufLen     = fmax(totmsgsize[0], MaxBufLen);

  double fieldSend[MaxBufLen]; // Place where we'd store and send
  double fieldRecv[MaxBufLen]; // Place where we'd recieve and move to phi

  // At the boundary partitions of the grid, there may be no neighbours.
  // At those boundaries, we'd just leave the Dirichlet BC as it is.
  // If we don't do that, since MPI doesn't put any halo layers on those sides,
  // we'd get a segmentation fault.

  // So, we figure out whether this is the case for our process
  // and create its left, right limits udim[0,dir] and udim[1,dir]
  // for all 3 directions dir
  int udim[2][p_dim];
  int disp = -1;
  // TO BE DEBUGGED BY UNIT TESTS
  for (int dir = 0; dir < p_dim; dir++)
  {
    int source, dest;
    ierr = MPI_Cart_shift(GRID_COMM_WORLD, dir, disp, &source, &dest);
    if (dest != MPI_PROC_NULL) // non-boundary, neighbour on 'left'
      udim[0][dir] = 1;
    else                       // boundary, no neighbour of 'left'
      udim[0][dir] = 2;
    if (source != MPI_PROC_NULL) // non-boundary, neighbour on 'right'
      udim[1][dir] = local_dim[dir];
    else                         // boundary, no neighbour on 'right'
      udim[1][dir] = local_dim[dir] - 1;
  } // udim[][0] -> Nk, udim[][1] -> Nj, udim[][2y] -> Ni


  int t0=0, t1=1;
  double maxdelta=0.0; // Difference b/w two jacobi iterates, measures convergence

  int tag = 0; // Unused?
  int itermax = 10000; // Doing 10 iterations for no reason, should use maxdelta
  double eps = 1e-10;
  int source, dest; // Neighbouring processors with which we'd trade.
                    // For example, if my cartesian rank is (1,1),
                    // my neighbours in 0 direction are (2,1) and (0,1).
                    // These will be obtained from MePI_Cart_shift
  int iter = 0;
  while (iter < itermax)
  {
    maxdelta = 0.0;
    for (int disp = -1; disp <= 1; disp = disp + 2) // disp = -1,1 to cover both directions
    {
      for (int dir = 0; dir < 3; dir++) // Looping over all directions.
      {
        MPI_Request req;
        MPI_Status status;
        ierr = MPI_Cart_shift(GRID_COMM_WORLD,
                              dir,
                              disp,
                              &source,
                              &dest); // Get neighbours' ranks in source and dest
        // If the boundaries are open, certain neighbours of the extreme
        // processors are to be ignored.
        if (source != MPI_PROC_NULL) // if source exists
        {

          MPI_Irecv(fieldRecv,         // Recieve buff
                    totmsgsize[dir], // Upper bound of recieve size
                    MPI_DOUBLE_PRECISION,
                    source,
                    tag,
                    GRID_COMM_WORLD,
                    &req
                    );
        }
        if (dest != MPI_PROC_NULL) // if destination exists
        {

        CopySendBuf(phi, t0,
                    disp, dir,
                    fieldSend, MaxBufLen); // Move send buffer to fieldSend
        ierr = MPI_Send(fieldSend, totmsgsize[dir], MPI_DOUBLE_PRECISION,
                 dest, tag, GRID_COMM_WORLD);
        }
        if (source != MPI_PROC_NULL)
        {
          ierr = MPI_Wait(&req, &status);
          CopyRecvBuf(phi, t0,
                      disp, dir,
                      fieldRecv, MaxBufLen);
        }
        Jacobi_sweep(udim, phi, t0, t1,
                     xmin, ymin, zmin, h, &maxdelta);
        ierr = MPI_Allreduce(MPI_IN_PLACE,
                             &maxdelta,
                             1,           // buffer size
                             MPI_DOUBLE,
                             MPI_MAX,
                             GRID_COMM_WORLD
                            );
        int tmp = t0; t0 = t1; t1 = tmp; // Swap t0 and t1
      }
    }
    iter += 1;
    if (myid==0)
      printf("iter = %d, eps = %.16f, maxdelta = %.16f\n", iter, eps, maxdelta);
    if (maxdelta < eps)
      break;
  }
  ierr = MPI_Finalize();
  printf("ierr = %d \n",ierr);
  return 0;
}

// Copy what I am sending to my neighbour to MaxBufLen
void CopySendBuf(Array3D phi[], int t0,
                 double disp, double dir, double fieldSend[], int MaxBufLen)
{
  int Ni = phi[t0].sizex()-1;
  int Nj = phi[t0].sizey()-1;
  int Nk = phi[t0].sizez()-1;
  int i1, i2, j1, j2, k1, k2; // left, right limits in appropriate directions
  int start = 0;
  assert ((dir >= 0 && dir <3 ) && "CopySendBuf: Incorrect dir");
  assert((disp == 1 || disp == -1) && "CopySendBuf: Incorrect disp");

  // Specifying the limits of moving phi to fieldSend[]
  if (dir == 0) // i-j plane
  {
    i1 = start + 1, i2 = Ni - 1; // non-ghost indices
    j1 = start + 1, j2 = Nj - 1; // non-ghost indices
    if (disp == -1)
      k1 = 1,    k2 = 1;         // non-ghost indices
    else
      k1 = Nk-1, k2 = Nk-1;      // non-ghost indices
  }
  else if (dir == 1) // i-k plane
  {
    i1 = start + 1, i2 = Ni - 1; // non-ghost indices
    k1 = start + 1, k2 = Nk - 1; // non-ghost indices
    if (disp == -1)
      j1 = 1,    j2 = 1;         // non-ghost indices
    else
      j1 = Nj-1, j2 = Nj - 1;    // non-ghost indices
  }
  else if (dir == 2) // j-k plane
  {
    j1 = start + 1, j2 = Nj - 1; // non-ghost indices
    k1 = start + 1, k2 = Nk - 1; // non-ghost indices
    if (disp == -1)
      i1 = 1,      i2 = 1;       // non-ghost indices
    else
      i1 = Ni - 1, i2 = Ni - 1;  // non-ghost indices
  }

  // After specifying the limits, move values to fieldSend[]
  int c = 0; // starting index for fieldSend
  for (int i = i1; i <= i2; i++)
    for (int j = j1; j <= j2; j++)
      for(int k = k1; k <= k2; k++)
      {
        fieldSend[c] = phi[t0](i,j,k);
        c = c+1;
      }
  assert( (c<MaxBufLen) && "CopySendBuff: SendBuff larger than expected.");
}

void CopyRecvBuf(Array3D phi[], // Solution old and new
                 int t0,
                 double disp, double dir, double fieldRecv[], int MaxBufLen)
{
  int Ni = phi[t0].sizex()-1;
  int Nj = phi[t0].sizey()-1;
  int Nk = phi[t0].sizez()-1;
  int i1, i2, j1, j2, k1, k2; // left, right limits in appropriate directions
  int start = 0;
  assert ((dir >= 0 && dir <3 ) && "CopySendBuf: Incorrect dir");
  assert((disp == 1 || disp == -1) && "CopySendBuf: Incorrect disp");

  // Specifying the limits of moving fieldRecv[] to phi
  if (dir == 0) // i-j plane
  {
    i1 = start + 1, i2 = Ni - 1; // non-ghost indices
    j1 = start + 1, j2 = Nj - 1; // non-ghost indices
    if (disp == 1)
      k1 = 0,    k2 = 0;     // ghost indices
    else
      k1 = Nk,   k2 = Nk;    // ghost indices
  }
  else if (dir == 1) // i-k plane
  {
    i1 = start + 1, i2 = Ni - 1; // non-ghost indices
    k1 = start + 1, k2 = Nk - 1; // non-ghost indices
    if (disp == 1)
      j1 = 0,    j2 = 0;   // ghost indices
    else
      j1 = Nj,   j2 = Nj;  // ghost indices
  }
  else if (dir == 2) // j-k plane
  {
    j1 = start + 1, j2 = Nj - 1; // non-ghost indices
    k1 = start + 1, k2 = Nk - 1; // non-ghost indices
    if (disp == 1)
      i1 = 0,      i2 = 0;  // ghost indices
    else
      i1 = Ni,     i2 = Ni; // ghost indices
  }

  int c = 0; // starting index for fieldRecv
  for (int i = i1; i <= i2; i++)
    for (int j = j1; j <= j2; j++)
      for(int k = k1; k <= k2; k++)
      {
        phi[t0](i,j,k) = fieldRecv[c];
        c = c+1;
      }
  assert( (c<MaxBufLen) && "CopySendBuff: SendBuff larger than expected.");
}

void Jacobi_sweep(int udim[][p_dim], // local_dim(pts in dimension)
                  Array3D phi[], int t0, int t1,
                  double xmin, double ymin, double zmin,
                  double h,
                  double *maxdelta)
{
  double x, y, z;
  // This is where the bug is.
  for (int i = udim[0][2]; i <= udim[1][2]; i++)
    for (int j = udim[0][1]; j <= udim[1][1]; j++)
      for (int k = udim[0][0]; k <= udim[1][0]; k++)
      {
        // printf("(%d,%d,%d)\n",i,j,k);
        x = xmin + i * h, y = ymin + j * h, z = zmin + k * h;
        phi[t1](i,j,k) = ( h*h * f(x,y,z)
                             + (phi[t0](i+1,j,k) + phi[t0](i-1,j,k))
                             + (phi[t0](i,j+1,k) + phi[t0](i,j-1,k))
                             + (phi[t0](i,j,k+1) + phi[t0](i,j,k-1)) ) / 6.0;
        *maxdelta = fmax(*maxdelta, fabs(phi[t1](i,j,k)-phi[t0](i,j,k)));
      }
  /* The code
    for (int i = udim[0][2]; i <= udim[1][2]; i++)
    for (int j = udim[0][1]; j <= udim[1][1]; j++)
      for (int k = udim[0][0]; k <= udim[1][0]; k++)
      {
        // printf("(%d,%d,%d)\n",i,j,k);
        x = xmin + i * h, y = ymin + j * h, z = zmin + k * h;
        phi[i][j][k][t1] = ( h*h * f(x,y,z)
                             + (phi[i+1][j][k][t0] + phi[i-1][j][k][t0])
                             + (phi[i][j+1][k][t0] + phi[i][j-1][k][t0])
                             + (phi[i][j][k+1][t0] + phi[i][j][k-1][t0]) ) / 6.0;
        *maxdelta = fmax(*maxdelta, fabs(phi[i][j][k][t1]-phi[i][j][k][t0]));
      }
   is working fine, so there's an issue with exchange? Hm...
   */
}


// What about splitting as 1 x 1 x np?
