#include<bits/stdc++.h>
#include "multigrid.h"
#include <jsoncpp/json/json.h>
#include "test.h"
using namespace std;

int main()
{
    int N=64;
    double h=1.0/N;
    vector<double> init;init.resize(N+1);
    Multigrid_Method<1> test(Fun5,bound_conditon(Dirichlet),restriction_operators(injection),interpolation_operators(linear),cycles(V_cycle),stopping_criteria(max_iteration),h,init,2.0);
    double la;
    for(int i=1;i<=1000;i++)
    {

    test.solve(2,2);
    cout<<test.error_norm(1)<<' '<<la/test.error_norm(1)<<endl;
    
    la=test.error_norm(1);
    // test.out();
    // test.out();
    }
    
    // test.out();
    return 0;
}