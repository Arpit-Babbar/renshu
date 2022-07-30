// Solve 3D Poisson's equation with source term f(x) = 1 using Jacobi's method

// Visualizing the steady state equation as time dependent equation where we
// are evolving from t0 to t1, we shall store solution_old in solution(:,t0)
// and solution_new in (:,t1) while running the iterates of Jacobi's method.

// Please see the tutorials on mpi_cart before this

// TODO - Order Ni, Nj, Nk correctly

#include <iostream>
#include <cmath>
#include <mpi.h>

#include "../include/array3d.h"

using namespace std;

// Define physical dimension globally
#define p_dim 3

double f(double x, double y, double z) // Source term
{
  double value = 1.0;
  return value;
}

// Move from solution array to intermediary array fieldSend before MPI_Send
void CopySendBuf(Array3D phi[],                // Solution old and new
                 int t0,
                 double disp, double dir, double fieldSend[], int MaxBufLen);

// Move from intermediate array fieldRecv to solution array after MPI_Recv
void CopyRecvBuf(Array3D phi[],                   // Solution old and new
                 int t0,
                 double disp, double dir, double fieldRecv[], int MaxBufLen);

// One Jacobi iteration
void Jacobi_sweep(int udim[][p_dim], // local_dim(pts in dimension)
                  Array3D phi[], int t0, int t1,
                  double xmin, double ymin, double zmin,
                  double h,
                  double *maxdelta);

int main(int argc, char** argv)
{
  int myid, numprocs, ierr; // rank, size renamed for problem
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  int N = 60;
  int spat_dim[p_dim];  // N x N x N grid
  int proc_dim[p_dim];  // np1 X np2 X np2
  int pbc_check[p_dim]; // 0/1 to indicate open/periodic bc in particular dimension

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
                   0,              // Root rank
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

 int local_dim[p_dim], mycoord[p_dim];

  // mycoord = Cartesian coordinate rank of process
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
    if (mycoord[i] < spat_dim[i] % proc_dim[i])
                     // total unused points
    {
      local_dim[i] += 1; // Adding unused point
    }
    //Total processes for which `if` holds = total unused points
  }

  int Ni = local_dim[2], Nj = local_dim[1], Nk = local_dim[0];
  vector<double> grid_z(Ni+2), grid_y(Nj+2), grid_x(Nk+2);
  printf("For rank %d with coordinates (%d,%d,%d), Nk, Nj, Ni = %d, %d, %d\n",
          myid, mycoord[0], mycoord[1], mycoord[2], Nk, Nj, Ni);
  Array3D phi[2];
  phi[0].resize(Ni+2,Nj+2,Nk+2); phi[1].resize(Ni+2,Nj+2,Nk+2);

  int totmsgsize[p_dim], MaxBufLen=0;

  // j-k plane
  totmsgsize[2] = local_dim[0] * local_dim[1]; // Nk*Nj
  MaxBufLen     = max(totmsgsize[2], MaxBufLen);
  // i-k plane
  totmsgsize[1] = local_dim[0] * local_dim[2]; // Nk*Ni
  MaxBufLen     = max(totmsgsize[1], MaxBufLen);
  // i-j plane
  totmsgsize[0] = local_dim[1] * local_dim[2]; // Nj*Ni
  MaxBufLen     = max(totmsgsize[0], MaxBufLen);
  double *fieldSend = new double [MaxBufLen];
  double *fieldRecv = new double [MaxBufLen];
  // char *filename1 = new char[filenameLength];
  // double fieldSend[MaxBufLen] = {0.0}; // usage: phi -> fieldSend -> MPI_Send
  // double fieldRecv[MaxBufLen] = {0.0}; // usage: MPI_Recv -> fieldRecv -> phi
  cout << "MaxBufLen = "<<MaxBufLen<<endl;
  // left, right physical limits, i.e., values to be updated
  int udim[2][p_dim] = {0};
  int disp = -1;
  for (int dir = 0; dir < p_dim; dir++)
  {
    int source, dest;
    ierr = MPI_Cart_shift(GRID_COMM_WORLD,
                          dir,             // physical direction
                          disp,            // displacement
                          &source,         // output negative dir proc
                          &dest            // output positive dir proc
                          );
    if (dest != MPI_PROC_NULL)   // non-boundary, neighbour on 'left'
      udim[0][dir] = 1;
    else                         // boundary, no neighbour of 'left'
      udim[0][dir] = 2;
    if (source != MPI_PROC_NULL) // non-boundary, neighbour on 'right'
      udim[1][dir] = local_dim[dir];
    else                         // boundary, no neighbour on 'right'
      udim[1][dir] = local_dim[dir] - 1;
  } // udim[][0] -> Nk, udim[][1] -> Nj, udim[][2] -> Ni


  int t0=0, t1=1;      // Indicate solution_old, solution_new in Jacobi
  double maxdelta=0.0; // Diff b/w 2 jacobi iterates to measure convergence

  int tag = 0;         // Unused
  int itermax = 1000;
  double eps = 1e-10;  // Tolerance
  int source, dest;    // Neighbouring processes with which we'd trade.
  int iter = 0;

  while (iter < itermax)
  {
    maxdelta = 0.0;
    for (int disp = -1; disp <= 1; disp = disp + 2) // disp = -1,1 for 2 directions
    {
      for (int dir = 0; dir < 3; dir++) // Looping over all directions.
      {
        MPI_Request req;
        MPI_Status status;
        ierr = MPI_Cart_shift(GRID_COMM_WORLD,
                              dir,
                              disp,
                              &source,
                              &dest
                              );
        if (source != MPI_PROC_NULL) // if source exists
        {
          // Recieve intermediary array fieldRecv
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
        // Copy phi to fieldSend to transfer ghost values
        CopySendBuf(phi, t0,               // solution old
                    disp, dir,             // neighbour indicator
                    fieldSend, MaxBufLen   // intermediary array, size
                    );
        ierr = MPI_Send(fieldSend, totmsgsize[dir], MPI_DOUBLE_PRECISION,
                        dest, tag, GRID_COMM_WORLD);
        }
        if (source != MPI_PROC_NULL)
        {
          ierr = MPI_Wait(&req, &status);
          // Copy fieldRecv to phi to recieve ghost values
          CopyRecvBuf(phi, t0,
                      disp, dir,
                      fieldRecv, MaxBufLen);
        }
      } // dir loop
    } // dis loop
    // Perform Jacobi iteration
    Jacobi_sweep(udim, phi, t0, t1,
                 xmin, ymin, zmin, h, &maxdelta);
    ierr = MPI_Allreduce(MPI_IN_PLACE,
                         &maxdelta,
                         1,           // buffer size
                         MPI_DOUBLE,
                         MPI_MAX,
                         GRID_COMM_WORLD
                        );
    iter += 1;
    if (myid==0)
      printf("iter = %d, eps = %.16f, maxdelta = %.16f\n", iter, eps, maxdelta);
    int tmp = t0; t0 = t1; t1 = tmp; // Swap t0 and t1
    if (maxdelta < eps)
      break;
  }
  ierr = MPI_Finalize();
  printf("ierr = %d \n",ierr);
  // delete[] fieldRecv;
  // delete[] fieldSend;
  // Even though these lines are good practise, they are costing a lot in speed.
  return 0;
}

// Move from solution array to intermediary array fieldSend before MPI_Send
void CopySendBuf(Array3D phi[], int t0,
                 double disp, double dir, double fieldSend[], int MaxBufLen)
{
  int iEnd = phi[t0].sizex()-1;
  int jEnd = phi[t0].sizey()-1;
  int kEnd = phi[t0].sizez()-1;
  int i1, i2, j1, j2, k1, k2; // left, right limits in appropriate directions
  int start = 0;
  assert ((dir >= 0 && dir <3 ) && "CopySendBuf: Incorrect dir");
  assert((disp == 1 || disp == -1) && "CopySendBuf: Incorrect disp");

  // Specifying the limits of moving phi to fieldSend[]
  if (dir == 0) // i-j plane
  {
    i1 = start + 1, i2 = iEnd - 1; // non-ghost indices
    j1 = start + 1, j2 = jEnd - 1; // non-ghost indices
    if (disp == -1)
      k1 = 1,      k2 = 1;         // non-ghost indices
    else
      k1 = kEnd-1, k2 = kEnd-1;    // non-ghost indices
  }
  else if (dir == 1) // i-k plane
  {
    i1 = start + 1, i2 = iEnd - 1; // non-ghost indices
    k1 = start + 1, k2 = kEnd - 1; // non-ghost indices
    if (disp == -1)
      j1 = 1,      j2 = 1;         // non-ghost indices
    else
      j1 = jEnd-1, j2 = jEnd - 1;  // non-ghost indices
  }
  else // j-k plane
  {
    j1 = start + 1, j2 = jEnd - 1;  // non-ghost indices
    k1 = start + 1, k2 = kEnd - 1;  // non-ghost indices
    if (disp == -1)
      i1 = 1,        i2 = 1;        // non-ghost indices
    else
      i1 = iEnd - 1, i2 = iEnd - 1; // non-ghost indices
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
  assert( (c<=MaxBufLen) && "CopySendBuff: SendBuff larger than expected.");
}

void CopyRecvBuf(Array3D phi[], // Solution old and new
                 int t0,
                 double disp, double dir, double fieldRecv[], int MaxBufLen)
{
  int iEnd = phi[t0].sizex()-1;
  int jEnd = phi[t0].sizey()-1;
  int kEnd = phi[t0].sizez()-1;
  int i1, i2, j1, j2, k1, k2; // left, right limits in appropriate directions
  int start = 0;
  assert ((dir >= 0 && dir <3 ) && "CopySendBuf: Incorrect dir");
  assert((disp == 1 || disp == -1) && "CopySendBuf: Incorrect disp");

  // Specifying the limits of moving fieldRecv[] to phi
  if (dir == 0) // i-j plane
  {
    i1 = start + 1, i2 = iEnd - 1; // non-ghost indices
    j1 = start + 1, j2 = jEnd - 1; // non-ghost indices
    if (disp == 1)
      k1 = 0,    k2 = 0;           // ghost indices
    else
      k1 = kEnd,   k2 = kEnd;      // ghost indices
  }
  else if (dir == 1) // i-k plane
  {
    i1 = start + 1, i2 = iEnd - 1; // non-ghost indices
    k1 = start + 1, k2 = kEnd - 1; // non-ghost indices
    if (disp == 1)
      j1 = 0,    j2 = 0;   // ghost indices
    else
      j1 = jEnd,   j2 = jEnd;  // ghost indices
  }
  else // j-k plane
  {
    j1 = start + 1, j2 = jEnd - 1; // non-ghost indices
    k1 = start + 1, k2 = kEnd - 1; // non-ghost indices
    if (disp == 1)
      i1 = 0,      i2 = 0;  // ghost indices
    else
      i1 = iEnd,     i2 = iEnd; // ghost indices
  }

  int c = 0; // starting index for fieldRecv
  for (int i = i1; i <= i2; i++)
    for (int j = j1; j <= j2; j++)
      for(int k = k1; k <= k2; k++)
      {
        phi[t0](i,j,k) = fieldRecv[c];
        c = c+1;
      }
  assert( (c<=MaxBufLen) && "CopySendBuff: SendBuff larger than expected.");
}

void Jacobi_sweep(int udim[][p_dim], // local_dim(pts in dimension)
                  Array3D phi[], int t0, int t1,
                  double xmin, double ymin, double zmin,
                  double h,
                  double *maxdelta)
{
  double x, y, z;
  for (int i = udim[0][2]; i <= udim[1][2]; i++)
    for (int j = udim[0][1]; j <= udim[1][1]; j++)
      for (int k = udim[0][0]; k <= udim[1][0]; k++)
      {
        x = xmin + i * h, y = ymin + j * h, z = zmin + k * h; // This is wrong!
        phi[t1](i,j,k) = ( h*h * f(x,y,z)
                             + (phi[t0](i+1,j,k) + phi[t0](i-1,j,k))
                             + (phi[t0](i,j+1,k) + phi[t0](i,j-1,k))
                             + (phi[t0](i,j,k+1) + phi[t0](i,j,k-1)) ) / 6.0;
        *maxdelta = fmax(*maxdelta, fabs(phi[t1](i,j,k)-phi[t0](i,j,k)));
      }
}
