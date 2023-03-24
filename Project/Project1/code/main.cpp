#include<bits/stdc++.h>
#include "BVP.h"
using namespace std;

class F1: public Function
{
    public:
    double operator()(const double &x,const double &y)const{return exp(y+sin(x));}
    double diff_x(const double &x,const double &y){return cos(x)*exp(y+sin(x));}
    double diff_y(const double &x,const double &y){return exp(y+sin(x));}
    double diff2_x(const double &x,const double &y){return -sin(x)*exp(y+sin(x))+cos(x)*cos(x)*exp(y+sin(x));}
    double diff2_y(const double &x,const double &y){return exp(y+sin(x));}
    double laplace(const double &x,const double &y){return diff2_x(x,y)+diff2_y(x,y);}
}F;
int main()
{
    double h=1.0/16;
    FD_regular FD(F,h,1);
    FD.solve();
    for(int i=0;i<1.0/h;i++)
    {
        for(int j=0;j<1.0/h;j++)printf("%lf ",FD(i*h,j*h));
        puts("");
    }
    return 0;
}