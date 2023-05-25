#include<bits/stdc++.h>
#include "AdamsBashforth.h"
// #include <jsoncpp/json/json.h>
using namespace std;

const double T1=17.06521656015796;
const double T2=19.14045706162071;

const vector<double> INIT={0.994, 0, 0, 0, -2.0015851063790825224, 0};
const vector<double> INIT2={0.87978, 0, 0, 0, -0.3797, 0};
const int n=2400000;
int main()
{
    AdamsBashforth_solver abm(1,T1,INIT,n);
    abm.solve();
    return 0;
}