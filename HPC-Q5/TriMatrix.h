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
    
    //Operator Overload: Calculate Multiplication of TriMatrix to Vector X
    vector<double> operator* (vector<double> X){
        double s=X.size();
        
        //diagm multiply X
        vector<double> Am(s);
        for (int i=0;i<s;i++){
            Am[i]=(*diagm)[i]*X[i];
        }
        
        //diagu multiply X
        vector<double> Au(s);
        for (int i=0;i<s-1;i++){
            Au[i]=(*diagu)[i]*X[i+1];
        }
        Au[s-1]=0;
        
        //diagl multiply X
        vector<double> Al(s);
        Al[0]=0;
        for (int i=1;i<s;i++){
            Al[i]=(*diagl)[i-1]*X[i-1];
        }
        
        //Superposition of Results
        vector<double> B(s);
        for (int i=0;i<s;i++){
            B[i]=Am[i]+Au[i]+Al[i];
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
