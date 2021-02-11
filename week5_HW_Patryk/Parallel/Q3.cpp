#include <iostream>
#include <mpi.h>
#include <math.h>  
#include <cmath> 
using namespace std;

int main(int argc, char* argv[]){

int rank;
int size;
int nproc= 6; // number of processors    
MPI_Comm comm;
comm = MPI_COMM_WORLD;  //Define communicator
MPI_Init(NULL,NULL); //Initialise parallel sequence
MPI_Comm_rank(comm, &rank);
MPI_Comm_size(comm, &size); 
    
      
if(size != nproc){ //Aborting if number of processors not equal to size
    cout<<"nproc = "<<nproc<<" does not fit number of processes requested which is "
    <<size<<". Restart with matching number!"<<endl;
    MPI_Abort(comm, 911);
    MPI_Finalize();
    return 0;
   } 
    
int M = nproc*(20-1)+2; // M length intervals
int N = 10000; // N time intervals
int K = (M-1)/nproc +2; // number of grid points to allocate for each process

    
double Usol[M+1];  // stores true solution    
double dx = 1./M;  // dx*M=1
double T = 0.1;  
//double T = atof(argv[1]);  // get final time from input argument    
double dt = T/N;
double dtdx = dt/(dx*dx);
cout<< "\ndx="<<dx<<", dt="<<dt<<", dt/dx²="<< dtdx<<endl;   
    
double prev[M+1], next[M+1]; // used to store numerical solution \in spacial dim at previous \& new timestep
double prev_grid[2][K]; //used to store solutions on grid points
double prev_no_boundary[M-1];
double prev_temp[K-2];    
    
prev[0] = 0;prev[M] = 0; next[0] = 0; next[M] = 0; // apply boundary conditions:
  
for (int j = 1; j < M; j++) {prev[j]= sin(2*M_PI*j*dx)+2*sin(5*M_PI*j*dx)+3*sin(20*M_PI*j*dx);}
        
double  t_start, t_end;
t_start = MPI_Wtime();
    
for (int j = 0; j < K; ++j){prev_grid[0][j] = prev[rank*(K-2)+ j];}   //Initialising grid solution points

// use numerical scheme to propagate \in time dimension and assign part of solution to each processor
for (int i = 1; i < N; i++){
    
   if (rank==0){ //calculate solution on interior points
        for (int j = 1; j < K - 1; ++j){
            prev_grid[1][j] = (dtdx*prev_grid[0][j - 1] + (1 - 2*dtdx)*prev_grid[0][j]+ dtdx*prev_grid[0][j + 1]);}   
        MPI_Ssend(&prev_grid[1][K - 2], 1, MPI_DOUBLE, rank + 1, 1, comm); //exchange solution between processors
        MPI_Recv(&prev_grid[1][K - 1], 1, MPI_DOUBLE, rank + 1, 0, comm, MPI_STATUS_IGNORE);}
                  
                   
   else if (rank == size - 1){               
        for (int j = 1; j < K - 1; ++j){ //calculate solution on inside of the segment
            prev_grid[1][j] = (dtdx*prev_grid[0][j - 1] + (1 - 2*dtdx)*prev_grid[0][j]+ dtdx*prev_grid[0][j + 1]);}     
        MPI_Recv(&prev_grid[1][0], 1, MPI_DOUBLE, rank - 1, 1, comm, MPI_STATUS_IGNORE); //exchange solution between processors
        MPI_Ssend(&prev_grid[1][1], 1, MPI_DOUBLE, rank - 1, 0, comm);}
    
   else{       
       for (int j = 1; j < K - 1; ++j){
            prev_grid[1][j] = (dtdx*prev_grid[0][j - 1] + (1 - 2*dtdx)*prev_grid[0][j]+ dtdx*prev_grid[0][j + 1]);}
   
         MPI_Recv(&prev_grid[1][0], 1, MPI_DOUBLE, rank - 1, 1, comm, MPI_STATUS_IGNORE);
         MPI_Ssend(&prev_grid[1][K - 2], 1, MPI_DOUBLE, rank + 1, 1, comm);
         MPI_Recv(&prev_grid[1][K - 1], 1, MPI_DOUBLE, rank + 1, 0, comm, MPI_STATUS_IGNORE);
         MPI_Ssend(&prev_grid[1][1], 1, MPI_DOUBLE, rank - 1, 0, comm);} 
   
    
      MPI_Barrier(comm); //Wait for all processes
      for (int j = 0; j < K; ++j){prev_grid[0][j] = prev_grid[1][j];}
      MPI_Barrier(comm); //Wait for all processes

   }
  
   //Gather solution
   for (int i = 0; i < K - 2; ++i){prev_temp[i] = prev_grid[1][i + 1];}
   MPI_Barrier(comm);  
   MPI_Gather(&prev_temp, K - 2, MPI_DOUBLE, prev_no_boundary, K - 2, MPI_DOUBLE, 0, comm);
   MPI_Barrier(comm);

   if (rank == 0){
      for (int m = 1; m < M; ++m){prev[m] = prev_no_boundary[m - 1];}   
      t_end = MPI_Wtime();
      cout << "Running time  " << t_end - t_start << endl;
    
      // print out array entries of numerical solution next to true solution
      cout << "\nTrue and numerical values at M="<<M<<" space points at time T="<<T<<":"<<endl;
      cout << "\nTrue values \n"<<endl;
      for(int m=0; m<=M; ++m){
          Usol[m] = exp(-4*M_PI*M_PI*T)*sin(2*M_PI*m*dx) + 2*exp(-25*M_PI*M_PI*T)*sin(5*M_PI*m*dx) + 3*exp(-400*M_PI*M_PI*T)*sin(20*M_PI*M*dx);
          cout << Usol[m] << endl;}
      cout << "\nNumerical solutions \n"<<endl;    
      for(int m=0; m<=M; ++m){
          Usol[m] = exp(-4*M_PI*M_PI*T)*sin(2*M_PI*m*dx) + 2*exp(-25*M_PI*M_PI*T)*sin(5*M_PI*m*dx) + 3*exp(-400*M_PI*M_PI*T)*sin(20*M_PI*M*dx);
          cout << prev[m] << endl;}}
    
    return 0;
}
