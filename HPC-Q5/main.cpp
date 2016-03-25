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
#include "TriMatrix.h"
using namespace std;

int main() {
    
    // INPUTS
    double L=1;
    double Nx=10000;
    double T=0.01;
    double Nt=2000;
    double alpha=0.001;
    
    // INITIAL CALCULATIONS
    double dt=T/Nt;
    double dx=L/Nx;
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
    
    //IMPLICIT TIME INTEGRATION
    double arg=0.5;
    TriMatrix ML(-arg*v,x.size());
    TriMatrix MR((1-arg)*v,x.size());
    vector <double> u0new(x.size());
    for (int i=0;i<Nt;i++){
        u0new=ML/(MR*u0);
        u0=u0new;
    }
    
    //OUTPUT RESULTS
    for (int i=0;i<x.size();i++){
        cout<<u0new[i]<<endl;
    }
    
    return 0;
}
