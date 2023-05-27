#ifndef IVP
#define IVP
#include<bits/stdc++.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
using namespace Eigen;
using namespace std;
const double EPS=1e-11;
void ERROR(string wrongmessage)
{
    cerr<<wrongmessage<<endl;
    exit(-1);
}
class IVPsolver
{
public:
    virtual void solve() = 0;  
};
class LMM : public IVPsolver
{
public:
    vector<double> alpha,beta;
};
class RKM : public IVPsolver
{
public:
    vector<double> RKweights,RKnodes;
    vector<vector<double> > RKMatrix;
};

template<typename T> vector<T> operator +(const vector<T> &a,const vector<T> &b)
{
    vector<T> re;re.clear();
    if(a.size()!=b.size())ERROR("wrong vector size");
    for(int i=0;i<(int)a.size();i++)re.push_back(a[i]+b[i]);
    return re;
}
template<typename T> vector<T> operator -(const vector<T> &a,const vector<T> &b)
{
    vector<T> re;re.clear();
    if(a.size()!=b.size())ERROR("wrong vector size");
    for(int i=0;i<(int)a.size();i++)re.push_back(a[i]-b[i]);
    return re;
}
template<typename T> vector<T> operator *(const T &a,const vector<T> &b)
{
    vector<T> re;re.clear();
    for(int i=0;i<(int)b.size();i++)re.push_back(a*b[i]);
    return re;
}
vector<double> f(const vector<double> &v,double t)
{
    double mu=1.0/81.45;
    vector<double> re;re.resize(6);
    re[0]=v[3],re[1]=v[4],re[2]=v[5];
    re[3]=2 * v[4] + v[0] - (mu * (v[0] + mu - 1)) / pow((pow(v[1], 2) + pow(v[2], 2) + pow(v[0] + mu - 1, 2)), 1.5) - ((1 - mu) * (v[0] + mu)) / pow((pow(v[1], 2) + pow(v[2], 2) + pow(v[0] + mu, 2)), 1.5);
    re[4]=-2 * v[3] + v[1] - (mu * v[1]) / pow((pow(v[1], 2) + pow(v[2], 2) + pow(v[0] + mu - 1, 2)), 1.5) - ((1 - mu) * v[1]) / pow((pow(v[1], 2) + pow(v[2], 2) + pow(v[0] + mu, 2)), 1.5);
    re[5]=-(mu * v[2]) / pow((pow(v[1], 2) + pow(v[2], 2) + pow(v[0] + mu - 1, 2)), 1.5) - ((1 - mu) * v[2]) / pow((pow(v[1], 2) + pow(v[2], 2) + pow(v[0] + mu, 2)), 1.5);
    return re;
}
double norm_inf(vector<double> x)
{
    double re=0;
    for(int i=0;i<(int)x.size();i++)re=re>abs(x[i])?re:abs(x[i]);
    return re;
}
double error_norm_inf(vector<double> uhat,vector<double> u,double e_abs,double e_rel)
{
    double re=0;
    if(uhat.size()!=u.size())ERROR("wrong vector size");
    for(int i=0;i<(int)u.size();i++)
        re=max((fabs(uhat[i]-u[i]))/(e_abs+abs(u[i])*e_rel),re);
    return re;
}
#endif