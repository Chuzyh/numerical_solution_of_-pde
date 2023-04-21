#include<bits/stdc++.h>
#include "multigrid.h"
#include <jsoncpp/json/json.h>
#include "test.h"
using namespace std;

int main()
{
    int N=128*2*2*2*2  ;
    double h=1.0/N;
    vector<double> init;init.resize(N+1);
    Multigrid_Method<1> test(Fun5,bound_conditon(mixed),restriction_operators(injection),interpolation_operators(linear),cycles(V_cycle),stopping_criteria(rela_accuracy),h,init,1e-8,3);
    Multigrid_Method<1> test2(Fun5,bound_conditon(mixed),restriction_operators(injection),interpolation_operators(linear),cycles(FMG),stopping_criteria(rela_accuracy),h,init,1e-8,3);
    double la;
    test.solve(3,3);
    test2.solve(10,10);
    
    
    cout<<test.error_norm(1)<<' '<<test2.error_norm(1)<<' ' <<test.error_norm(1)/test2.error_norm(1)<<endl;
    cout<<test.residual_norm(1)<<' '<<test2.residual_norm(1)<<' ' <<test.residual_norm(1)/test2.residual_norm(1)<<endl;
    return 0;
}