#include <iostream>
#include <mpi.h>
#include <math.h>  
#include <cmath> 
using namespace std;

int main(int argc, char* argv[]){

int M = 100; // M length intervals
int N = 10000; // N time intervals
double Usol[M+1];  // stores true solution    
double dx = 1./M;  // dx*M=1
double T = 0.1;  
//double T = atof(argv[1]);  // get final time from input argument    
double dt = T/N;
double dtdx = dt/(dx*dx);
cout<< "\ndx="<<dx<<", dt="<<dt<<", dt/dx²="<< dtdx<<endl;   
    
double prev[M+1], next[M+1]; // used to store numerical solution \in spacial dim at previous \& new timestep
prev[0] = 0;prev[M] = 0; next[0] = 0; next[M] = 0; // apply boundary conditions:
  
  for (int j = 1; j < M; j++) {prev[j]= sin(2*M_PI*j*dx)+2*sin(5*M_PI*j*dx)+3*sin(20*M_PI*j*dx);}
        
  double  t_start, t_end;
  t_start = MPI_Wtime();

  // use numerical scheme to propagate \in time dimension over M+1 space points 
  for (int i = 1; i <= N; i++) {
    for (int j = 1; j < M; j++){next[j] = prev[j] + dtdx * (prev[j - 1] - 2 * prev[j] + prev[j + 1]);}
        for (int j = 1; j < M; j++) {prev[j] = next[j];}
        }

t_end = MPI_Wtime();
cout << "Running time  " << t_end - t_start << endl;
    
// print out array entries of numerical solution next to true solution
cout << "\nTrue and numerical values at M="<<M<<" space points at time T="<<T<<":"<<endl;
cout << "\nTrue values           Numerical solutions\n"<<endl;
for(int m=0; m<=M; ++m){
    Usol[m] = exp(-4*M_PI*M_PI*T)*sin(2*M_PI*m*dx) + 2*exp(-25*M_PI*M_PI*T)*sin(5*M_PI*m*dx) + 3*exp(-400*M_PI*M_PI*T)*sin(20*M_PI*M*dx);
    cout << Usol[m] << "            " << next[m] << endl;
    // note that we did not really need to store the true solution in the array just to print out the values.
}
    
return 0;}
