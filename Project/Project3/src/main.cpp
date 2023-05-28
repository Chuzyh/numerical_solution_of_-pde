#include<bits/stdc++.h>
#include "AdamsBashforth.h"
#include "AdamsMoulton.h"
#include "BackDifferFormula.h"
#include "classicalRK.h"
#include "ESDIRK.h"
#include "GaussLegendreRK.h"
#include "FehlbergRK.h"
#include "DormandPrinceRK.h"
#include <jsoncpp/json/json.h>
using namespace std;

const double T1=17.06521656015796;
const double T2=19.14045706162071;

const vector<double> INIT={0.994, 0, 0, 0, -2.0015851063790825224, 0};
const vector<double> INIT2={0.87978, 0, 0, 0, -0.3797, 0};
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
void test1()
{
    ifstream ifs("test.json");
    Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj);
    for(int i=0;i<obj["method"].size();i++)
    for(int k=0;k<obj["p"][i].size();k++)
    {
        vector<double> ERR,tmp,ti;
        
        for(int j=0;j<obj["tick"].size();j++)
        {
            cerr<<obj["method"][i].asString()<<' '<<obj["tick"][j].asDouble()<<' '<<obj["p"][i][k].asInt()<<endl;
            
            clock_t st1 = clock();
            classFactory &F=classFactory::createFactory();
            auto sol=F.createMethod(obj["method"][i].asString());
            (*sol).set(obj["p"][i][k].asInt(),T1,INIT,T1/obj["tick"][j].asDouble()*(obj["method"][i].asString()=="ESDIRK"?20:1));
            (*sol).solve();
            clock_t ed2 = clock();
            tmp=(*sol).error(INIT);
            ERR.push_back(max(abs(tmp[0]),abs(tmp[1])));
            ti.push_back((double)(ed2-st1)/CLOCKS_PER_SEC);
            printf("%.18lf %.18lf\n",norm_inf(tmp),(double)(ed2-st1)/CLOCKS_PER_SEC);
            
        }
        freopen(("../data/"+obj["method"][i].asString()+to_string(obj["p"][i][k].asInt())+"test1").c_str(),"w",stdout);
        for(int j=0;j<obj["tick"].size();j++)printf("%.18lf %.18lf %lf\n",ERR[j],ti[j],obj["tick"][j].asDouble());
        fclose(stdout);
    }
    
}
void test2()
{
    ifstream ifs("test.json");
    Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj);
    for(int i=0;i<obj["method"].size();i++)
    for(int k=0;k<obj["p"][i].size();k++)
    {
        vector<double> ERR,tmp,ti;
        
        for(int j=0;j<obj["tick"].size();j++)
        {
            cerr<<obj["method"][i].asString()<<' '<<obj["tick"][j].asDouble()<<' '<<obj["p"][i][k].asInt()<<endl;
            
            clock_t st1 = clock();
            classFactory &F=classFactory::createFactory();
            auto sol=F.createMethod(obj["method"][i].asString());
            (*sol).set(obj["p"][i][k].asInt(),T2,INIT2,T2/obj["tick"][j].asDouble()*(obj["method"][i].asString()=="ESDIRK"?20:1));
            (*sol).solve();
            clock_t ed2 = clock();
            tmp=(*sol).error(INIT2-INIT2);
            sol=F.createMethod("FehlbergRK");
            (*sol).set(5,T2,INIT2,T2/obj["tick"][j].asDouble()*(obj["method"][i].asString()=="ESDIRK"?20:1));
            (*sol).solve();
            tmp=(*sol).error(tmp);
            ERR.push_back(max(abs(tmp[0]),abs(tmp[1])));
            ti.push_back((double)(ed2-st1)/CLOCKS_PER_SEC);
            printf("%.18lf %.18lf\n",norm_inf(tmp),(double)(ed2-st1)/CLOCKS_PER_SEC);
            
        }
        freopen(("../data/"+obj["method"][i].asString()+to_string(obj["p"][i][k].asInt())+"test2").c_str(),"w",stdout);
        for(int j=0;j<obj["tick"].size();j++)printf("%.18lf %.18lf %lf\n",ERR[j],ti[j],obj["tick"][j].asDouble());
        fclose(stdout);
    }
    
}
int main()
{
    fac_init();
    test1();
    test2();
    //  classFactory &F=classFactory::createFactory();
    // auto sol=F.createMethod("AdamsBashforth");
    // (*sol).set(1,T1,INIT,24000);
    // (*sol).solve();
    // sol=F.createMethod("classicalRK");
    // (*sol).set(1,T1,INIT,6000);
    // (*sol).solve();
    // sol=F.createMethod("DormandPrinceRK");
    // (*sol).set(5,-1,INIT,100);
    // (*sol).solve();
    return 0;
}