#include<bits/stdc++.h>
#include "BVP.h"
using namespace std;

class F1: public Function
{
    public:
    double operator()(const double &x,const double &y)const{return exp(y+sin(x));}
    double diff_x(const double &x,const double &y)const{return cos(x)*exp(y+sin(x));}
    double diff_y(const double &x,const double &y)const{return exp(y+sin(x));}
    double diff2_x(const double &x,const double &y)const{return -sin(x)*exp(y+sin(x))+cos(x)*cos(x)*exp(y+sin(x));}
    double diff2_y(const double &x,const double &y)const{return exp(y+sin(x));}
}F;
int main()
{
    double h=1.0/128;
    FD_regular FD(F,h,3);
    FD.solve();
    for(int i=0;i<1.0/h;i++)
    {
        for(int j=0;j<1.0/h;j++)printf("%lf ",FD(i*h,j*h));
        puts("");
    }
    cerr<<FD.error_norm(1)<<' '<<FD.error_norm(2)<<' '<<FD.error_norm(-1)<<endl;
    return 0;
}