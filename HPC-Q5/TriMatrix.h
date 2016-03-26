//
//  TriMatrix.h
//  HPC-Q3c
//
//  Created by hyo13 on 23/3/16.
//  Copyright (c) 2016 hyo13. All rights reserved.
//

#ifndef CLASS_TriMatrix
#define CLASS_TriMatrix

#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

class TriMatrix{
    
private:
    vector <double> *diagm, *diagu, *diagl;
    
public:
    TriMatrix(double v,double size){
        diagm = new vector <double> (size);
        diagu = new vector <double> (size-1);
        diagl = new vector <double> (size-1);
        
        //create diagm vector
        (*diagm)[0]=1;
        (*diagm)[size-1]=1;
        for (int i=1;i<size-1;i++){
            (*diagm)[i]=1-2*v;
        }
        
        //create diagu vector
        (*diagu)[0]=0;
        for (int i=1;i<size-1;i++){
            (*diagu)[i]=v;
        }
        
        //create diagl vector
        (*diagl)[size-2]=0;
        for (int i=0;i<size-2;i++){
            (*diagl)[i]=v;
        }
        
    }
    
    //Function: Calculate Multiplication of TriMatrix to Vector X
    vector<double> MUL(vector<double> X, int rank){
        double s=X.size();
        vector <double> B(s);
        //RANK 0 PROCESS
        if (rank==0){
            //generate upper B
            B[0]=(*diagm)[0]*X[0]+(*diagu)[0]*X[1];
            for (int i=1;i<s/2;i++){
                B[i]=(*diagl)[i-1]*X[i-1]+(*diagm)[i]*X[i]+(*diagu)[i]*X[i+1];
                MPI_Send(&B[i],1,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
            }
            for (int i=s/2;i<s-1;i++){
                MPI_Recv(&B[i],1,MPI_DOUBLE,1,1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        //RANK 1 PROCESS
        else {
            //generate lower B
            for (int i=s/2;i<s-1;i++){
                B[i]=(*diagl)[i-1]*X[i-1]+(*diagm)[i]*X[i]+(*diagu)[i]*X[i+1];
                MPI_Send(&B[i],1,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
            }
            B[s-1]=(*diagl)[s-2]*X[s-2]+(*diagm)[s-1]*X[s-1];
            for (int i=1;i<s/2;i++){
                MPI_Recv(&B[i],1,MPI_DOUBLE,0,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
        }
        return B;
    }
    
    //Operator Overload: Matrix-Vector Solve Operation
    vector<double> operator/ (vector<double> B){
        double s=B.size();
        vector<double> M=(*diagm);
        vector<double> U=(*diagu);
        vector<double> L=(*diagl);
        
        //Forward Elimination
        for (int i=1;i<s;i++){
            double m=L[i-1]/M[i-1];
            M[i]=M[i]-m*U[i-1];
            B[i]=B[i]-m*B[i-1];
        }
        //Calculate Xn
        vector<double> X(s);
        X[s-1]=B[s-1]/M[s-1];
        
        //Backward Substitution
        for (int i=s-2;i>=0;i--){
            X[i]=(B[i]-U[i]*X[i+1])/M[i];
        }
        
        return X;
    }
    
};

#endif /* defined(__HPC_Q1__TriMatrix__) */
