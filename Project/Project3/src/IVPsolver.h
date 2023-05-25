#ifndef IVP
#define IVP
#include<bits/stdc++.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
using namespace Eigen;
using namespace std;
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
private:
    vector<double> alpha,beta;
};
#endif