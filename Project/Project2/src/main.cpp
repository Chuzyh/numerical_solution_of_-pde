#include<bits/stdc++.h>
#include "multigrid.h"
#include <jsoncpp/json/json.h>
#include "test.h"
using namespace std;

void test_fun1d(Function & f,int num,string filename)
{
    ifstream ifs(filename);
    Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj);
    freopen(("../data/four_grid_convergence_rate_of_fun"+to_string(num)).c_str(),"w",stdout);
    double eps=obj["esplion"].asDouble();
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
            
            Multigrid_Method<1> test(f,BC,RO,IO,CY,stopping_criteria(rela_accuracy),h,init,eps,3);
            test.solve(5,5,0);
            double now_error=test.error_norm(1);
            double now_resi=test.residual_norm(1);
            
            printf(" %d & %.3e & %.3lf & %.3e & %.3lf\\\\ \\hline \n",N,now_error,la_error/now_error,now_resi,la_resi/now_resi);
            la_error=now_error;la_resi=now_resi;
        }
    }
    
    fclose(stdout);
}
void test_cycle(Function & f,int num,string filename)
{
    ifstream ifs(filename);
    Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj);
    freopen(("../data/cycle_iteration_convergence_rate_of_fun"+to_string(num)).c_str(),"w",stdout);
    double eps=obj["esplion"].asDouble();
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
        
            Multigrid_Method<1> test(f,BC,RO,IO,CY,stopping_criteria(rela_accuracy),h,init,eps,3);
            test.solve(5,5,1);
        }
    }
    
    fclose(stdout);
}

void test_fun2d(Function2d & f,int num,string filename)
{
    ifstream ifs(filename);
    Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj);
    freopen(("../data/four_grid_convergence_rate_of_fun2d"+to_string(num)).c_str(),"w",stdout);
    double eps=obj["esplion"].asDouble();
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
            
            Multigrid_Method<2> test(f,BC,RO,IO,CY,stopping_criteria(rela_accuracy),h,init,eps,3);
            test.solve(5,5,0);
            double now_error=test.error_norm(1);
            double now_resi=test.residual_norm(1);
            
            printf(" %d & %.3e & %.3lf & %.3e & %.3lf\\\\ \\hline \n",N,now_error,la_error/now_error,now_resi,la_resi/now_resi);
            la_error=now_error;la_resi=now_resi;
        }
    }
    
    fclose(stdout);
}
void test_cycle2d(Function2d & f,int num,string filename)
{
    ifstream ifs(filename);
    Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj);
    freopen(("../data/cycle_iteration_convergence_rate_of_fun2d"+to_string(num)).c_str(),"w",stdout);
    double eps=obj["esplion"].asDouble();
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
        
            Multigrid_Method<2> test(f,BC,RO,IO,CY,stopping_criteria(rela_accuracy),h,init,eps,3);
            test.solve(5,5,1);
        }
    }
    
    fclose(stdout);
}
void test_time2d(Function2d & f,int num,string filename)
{
    ifstream ifs(filename);
    Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj);
    freopen(("../data/cpu_time_of_mutigrid_and_LU_of_fun2d"+to_string(num)).c_str(),"w",stdout);
    double eps=obj["esplion"].asDouble();
    for(int i1=0;i1<obj["condition"].size();i1++)
    for(int i4=0;i4<obj["cycles"].size();i4++)
    {
        double la_error=0,la_resi=0;
        bound_conditon BC=bound_conditon_string(obj["condition"][i1].asString());
        restriction_operators RO=full_weighting;
        interpolation_operators IO=quadratic;
        cycles CY=cycles_string(obj["cycles"][i4].asString());
        cout<<"N"<<"  "<<obj["condition"][i1].asString()<<' '<<obj["cycles"][i4].asString()<<endl;
        
        for(int i=0;i<obj["n"].size();i++)
        {
            int N=obj["n"][i].asInt();
            double h=1.0/N;
            vector<vector<double> > init;init.resize(N+1);
            for(int i=0;i<=N;i++)init[i].resize(N+1);
            
            clock_t st1 = clock();
            Multigrid_Method<2> test(f,BC,RO,IO,CY,stopping_criteria(rela_accuracy),h,init,eps,3);
            test.solve(5,5,0);
            clock_t st2 = clock();
            FD_regular test2(f,h,(int)BC+1);
            test2.solve();
            clock_t ed = clock();

            cout<<N<<" & "<<(double)(st2-st1)/CLOCKS_PER_SEC<<" & "<<(double)(ed-st2)/CLOCKS_PER_SEC<<" & "<<(double)(ed-st2)/(st2-st1)<<"\\\\ \\hline "<<endl;
        }
    }
    
    fclose(stdout);
}
int main()
{
    test_fun1d(Fun5,1,"test.json");
    test_fun1d(Fun6,2,"test.json");
    test_fun1d(Fun7,3,"test.json");
    test_cycle(Fun5,1,"test.json");
    test_cycle(Fun6,2,"test.json");
    test_cycle(Fun7,3,"test.json");


    test_fun2d(Fun,1,"test2.json");
    test_fun2d(Fun2,2,"test2.json");
    test_fun2d(Fun3,3,"test2.json");
    test_cycle2d(Fun,1,"test2.json");
    test_cycle2d(Fun2,2,"test2.json");
    test_cycle2d(Fun3,3,"test2.json");
    
    test_time2d(Fun,1,"test2.json");
    test_time2d(Fun2,2,"test2.json");
    test_time2d(Fun3,3,"test2.json");
    return 0;
}