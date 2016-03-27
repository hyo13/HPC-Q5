//
//  main.cpp
//  HPC-Q3c
//
//  Created by hyo13 on 23/3/16.
//  Copyright (c) 2016 hyo13. All rights reserved.
//

#include <iostream>
#include <vector>
#include <math.h>
#include "mpi.h"
#include "TriMatrix.h"
using namespace std;

int main(int argc, char* argv[]) {
    
    // INPUTS
    double L=1;
    double Nx=20;
    double T=1;
    double Nt=1000;
    double alpha=1;
    
    //calculate minimum input time step for Forward Euler to converge with v = or < 0.5
    double dx=L/Nx;
    int Ntmin=2*alpha*T/(pow(dx,2));
    
    //check that Nt is > or = Ntmax to make sure Forward Euler converges
    if (Nt<Ntmin){
        cout<<"ERROR: Nt is smaller than minimum time step allowed. Program Terminated."<<endl;
        terminate();
    }
    if (fmod(Nx,2)==1){
        cout<<"ERROR: Nx cannot be divided by 2 for parallel programming. Program Terminated."<<endl;
        terminate();
    }
    
    // INITIAL CALCULATIONS
    double dt=T/Nt;
    double v=alpha*dt/pow(dx,2);
    
    // X-POSIION VECTOR
    vector <double> x(Nx+1);
    for (int i=0;i<Nx+1;i++){
        x[i]=0+i*dx;
    }
    
    //INITIAL CONDITION
    vector <double> u0(x.size());
    u0[1]=0;
    u0[x.size()]=0;
    for (int i=1;i<x.size()-1;i++){
        u0[i]=x[i]*(1-x[i]);
    }
    
    //CREATE TRI-DIAGONAL MATRIX
    double arg=0.5;
    TriMatrix ML(-arg*v,x.size());
    TriMatrix MR((1-arg)*v,x.size());
    
    //SETUP FOR PARALLEL PROGRAMMING
    MPI_Init(&argc,&argv);
    int rank;
    int retval;
    retval=MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    //IMPLICIT TIME INTEGRATION
    vector <double> u0new(x.size());
    for (int i=0;i<Nt;i++){
        u0new=ML/(MR.MUL(u0,rank));
        u0=u0new;
        cout<<i<<endl;
    }
    
    //END PARALLEL PROGRAMMING
    MPI_Finalize();
    
    //OUTPUT RESULTS
    for (int i=0;i<x.size();i++){
        cout<<u0new[i]<<endl;
    }
    
    return 0;
}
