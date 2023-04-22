#include<bits/stdc++.h>
#include "multigrid.h"
#include <jsoncpp/json/json.h>
#include "test.h"
using namespace std;

/*void test_fun1d(Function & f,int num)
{
    ifstream ifs("test.json");
    Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj);
    freopen(("../data/four_grid_convergence_rate_of_fun"+to_string(num)).c_str(),"w",stdout);
    for(int i1=0;i1<obj["condition"].size();i1++)
    for(int i2=0;i2<obj["restriction_operators"].size();i2++)
    for(int i3=0;i3<obj["interpolation_operators"].size();i3++)
    for(int i4=0;i4<obj["cycles"].size();i4++)
    {
        double la_error=0,la_resi=0;
        bound_conditon BC=bound_conditon_string(obj["condition"][i1].asString());
        restriction_operators RO=restriction_operators_string(obj["restriction_operators"][i2].asString());
        interpolation_operators IO=interpolation_operators_string(obj["interpolation_operators"][i3].asString());
        cycles CY=cycles_string(obj["cycles"][i4].asString());
        cout<<"N  "<<obj["condition"][i1].asString()<<' '<<obj["restriction_operators"][i2].asString()<<' '<<obj["interpolation_operators"][i3].asString()<<' '<<obj["cycles"][i4].asString()<<endl;
        for(int i=0;i<obj["n"].size();i++)
        {
            int N=obj["n"][i].asInt();
            double h=1.0/N;
            vector<double> init;init.resize(N+1);
            
            Multigrid_Method<1> test(f,BC,RO,IO,CY,stopping_criteria(rela_accuracy),h,init,1e-8,3);
            test.solve(5,5,0);
            double now_error=test.error_norm(1);
            double now_resi=test.residual_norm(1);
            
            printf(" %d & %.3e & %.3lf & %.3e & %.3lf\\\\ \\hline \n",N,now_error,la_error/now_error,now_resi,la_resi/now_resi);
            la_error=now_error;la_resi=now_resi;
        }
    }
    
    fclose(stdout);
}
void test_cycle(Function & f,int num)
{
    ifstream ifs("test.json");
    Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj);
    freopen(("../data/cycle_iteration_convergence_rate_of_fun"+to_string(num)).c_str(),"w",stdout);
    for(int i1=0;i1<obj["condition"].size();i1++)
    for(int i2=0;i2<obj["restriction_operators"].size();i2++)
    for(int i3=0;i3<obj["interpolation_operators"].size();i3++)
    for(int i4=0;i4<obj["cycles"].size();i4++)
    {
        double la_error=0,la_resi=0;
        bound_conditon BC=bound_conditon_string(obj["condition"][i1].asString());
        restriction_operators RO=restriction_operators_string(obj["restriction_operators"][i2].asString());
        interpolation_operators IO=interpolation_operators_string(obj["interpolation_operators"][i3].asString());
        cycles CY=cycles_string(obj["cycles"][i4].asString());
        for(int i=0;i<obj["n"].size();i++)
        {
            int N=obj["n"][i].asInt();
            double h=1.0/N;
            vector<double> init;init.resize(N+1);
            cout<<N<<"  "<<obj["condition"][i1].asString()<<' '<<obj["restriction_operators"][i2].asString()<<' '<<obj["interpolation_operators"][i3].asString()<<' '<<obj["cycles"][i4].asString()<<endl;
        
            Multigrid_Method<1> test(f,BC,RO,IO,CY,stopping_criteria(rela_accuracy),h,init,1e-8,3);
            test.solve(5,5,1);
        }
    }
    
    fclose(stdout);
}

void test_fun2d(Function2d & f,int num)
{
    ifstream ifs("test.json");
    Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj);
    freopen(("../data/four_grid_convergence_rate_of_fun2d"+to_string(num)).c_str(),"w",stdout);
    for(int i1=0;i1<obj["condition"].size();i1++)
    for(int i2=0;i2<obj["restriction_operators"].size();i2++)
    for(int i3=0;i3<obj["interpolation_operators"].size();i3++)
    for(int i4=0;i4<obj["cycles"].size();i4++)
    {
        double la_error=0,la_resi=0;
        bound_conditon BC=bound_conditon_string(obj["condition"][i1].asString());
        restriction_operators RO=restriction_operators_string(obj["restriction_operators"][i2].asString());
        interpolation_operators IO=interpolation_operators_string(obj["interpolation_operators"][i3].asString());
        cycles CY=cycles_string(obj["cycles"][i4].asString());
        cout<<"N  "<<obj["condition"][i1].asString()<<' '<<obj["restriction_operators"][i2].asString()<<' '<<obj["interpolation_operators"][i3].asString()<<' '<<obj["cycles"][i4].asString()<<endl;
        for(int i=0;i<obj["n"].size();i++)
        {
            int N=obj["n"][i].asInt();
            double h=1.0/N;
            vector<vector<double> > init;init.resize(N+1);
            for(int i=0;i<=N;i++)init[i].resize(N+1);
            
            Multigrid_Method<2> test(f,BC,RO,IO,CY,stopping_criteria(rela_accuracy),h,init,1e-8,3);
            test.solve(5,5,0);
            double now_error=test.error_norm(1);
            double now_resi=test.residual_norm(1);
            
            printf(" %d & %.3e & %.3lf & %.3e & %.3lf\\\\ \\hline \n",N,now_error,la_error/now_error,now_resi,la_resi/now_resi);
            la_error=now_error;la_resi=now_resi;
        }
    }
    
    fclose(stdout);
}
void test_cycle2d(Function2d & f,int num)
{
    ifstream ifs("test.json");
    Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj);
    freopen(("../data/cycle_iteration_convergence_rate_of_fun2d"+to_string(num)).c_str(),"w",stdout);
    for(int i1=0;i1<obj["condition"].size();i1++)
    for(int i2=0;i2<obj["restriction_operators"].size();i2++)
    for(int i3=0;i3<obj["interpolation_operators"].size();i3++)
    for(int i4=0;i4<obj["cycles"].size();i4++)
    {
        double la_error=0,la_resi=0;
        bound_conditon BC=bound_conditon_string(obj["condition"][i1].asString());
        restriction_operators RO=restriction_operators_string(obj["restriction_operators"][i2].asString());
        interpolation_operators IO=interpolation_operators_string(obj["interpolation_operators"][i3].asString());
        cycles CY=cycles_string(obj["cycles"][i4].asString());
        for(int i=0;i<obj["n"].size();i++)
        {
            int N=obj["n"][i].asInt();
            double h=1.0/N;
            vector<vector<double> > init;init.resize(N+1);
            for(int i=0;i<=N;i++)init[i].resize(N+1);
            cout<<N<<"  "<<obj["condition"][i1].asString()<<' '<<obj["restriction_operators"][i2].asString()<<' '<<obj["interpolation_operators"][i3].asString()<<' '<<obj["cycles"][i4].asString()<<endl;
        
            Multigrid_Method<2> test(f,BC,RO,IO,CY,stopping_criteria(rela_accuracy),h,init,1e-8,3);
            test.solve(5,5,1);
        }
    }
    
    fclose(stdout);
}*/
int main()
{
    /*test_fun1d(Fun5,1);
    test_fun1d(Fun6,2);
    test_fun1d(Fun7,3);
    test_cycle(Fun5,1);
    test_cycle(Fun6,2);
    test_cycle(Fun7,3);*/
    //test_fun2d(Fun,1);
    int N=64  ;
    double h=1.0/N;
    vector<vector<double> > init;init.resize(N+1);
    for(int i=0;i<=N;i++)init[i].resize(N+1);

    Multigrid_Method<2> test(Fun,bound_conditon(Neumann),restriction_operators(injection),interpolation_operators(linear),cycles(V_cycle),stopping_criteria(rela_accuracy),h,init,1e-8,3);
    Multigrid_Method<2> test2(Fun,bound_conditon(Neumann),restriction_operators(injection),interpolation_operators(linear),cycles(FMG),stopping_criteria(rela_accuracy),h,init,1e-8,3);
    double la;
    test.solve(5,5,1);
    test2.solve(5,5,1);
    
    
    cout<<test.error_norm(1)<<' '<<test2.error_norm(1)<<' ' <<test.error_norm(1)/test2.error_norm(1)<<endl;
    cout<<test.residual_norm(1)<<' '<<test2.residual_norm(1)<<' ' <<test.residual_norm(1)/test2.residual_norm(1)<<endl;
    return 0;
}