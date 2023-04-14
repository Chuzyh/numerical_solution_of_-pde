#include<bits/stdc++.h>
#include "multigrid.h"
#include <jsoncpp/json/json.h>
#include "test.h"
using namespace std;

int main()
{
    int N=4096*2*2*2*2;
    double h=1.0/N;
    vector<double> init;init.resize(N+1);
    Multigrid_Method<1> test(Fun4,bound_conditon(Dirichlet),restriction_operators(injection),interpolation_operators(linear),cycles(V_cycle),stopping_criteria(max_iteration),h,init,5.0);
    test.Vcycle(10,10);
    cout<<test.error_norm(1)<<endl;
    // test.out();
    return 0;
}