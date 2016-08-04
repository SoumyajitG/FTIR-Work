#include<iostream>
#include<cmath>
#include <stdlib.h>
#include<memory>
#include"mex.h"
#include "C:\Users\Lenovo\Documents\Visual Studio 2012\Projects\eigenpinv\eigenpinv\Eigen\Eigen" 
#include<windows.h>
using namespace std;
using namespace Eigen;
void OMP(double *btmp,double *Atmp,int M,int N,int K,double *result)//
{
    Vec