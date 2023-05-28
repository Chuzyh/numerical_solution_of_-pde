#include<bits/stdc++.h>
#include "AdamsBashforth.h"
#include "AdamsMoulton.h"
#include "BackDifferFormula.h"
#include "classicalRK.h"
#include "ESDIRK.h"
#include "GaussLegendreRK.h"
#include "FehlbergRK.h"
#include "DormandPrinceRK.h"
// #include <jsoncpp/json/json.h>
using namespace std;

const double T1=17.06521656015796;
const double T2=19.14045706162071;

const vector<double> INIT={0.994, 0, 0, 0, -2.0015851063790825224, 0};
const vector<double> INIT2={0.87978, 0, 0, 0, -0.3797, 0};
const int n=10000;
void fac_init()
{
    register_ABF();
    register_ADM();
    register_BDF();
    register_CRK();
    register_FRK();
    register_DPRK();
    register_ESDIRK();
    register_GaussLegendreRK();
}
int main()
{
    fac_init();
    classFactory &F=classFactory::createFactory();
    auto sol=F.createMethod("classicalRK");
    
    (*sol).set(1,T1,INIT,n);
    (*sol).solve();
    (*sol).set(1,T2,INIT2,n);
    (*sol).solve();
    return 0;
}