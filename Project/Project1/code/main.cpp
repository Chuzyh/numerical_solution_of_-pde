#include<bits/stdc++.h>
#include "BVP.h"
#include <jsoncpp/json/json.h>
#include "test.h"
using namespace std;

int main()
{
    ifstream ifs("test.json");
    Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj);
    Circle C((point){obj["Circle center"][0].asDouble(),obj["Circle center"][1].asDouble()},obj["Circle radius"].asDouble());
    for(int i=0;i<obj["n"].size();i++)
    {
        int n=obj["n"][i].asInt();
        double h=1.0/n;
        FD_regular FD1(Fun,h,1);
        FD1.solve();
        freopen(("../data/fun1"+to_string(n)).c_str(),"w",stdout);
        for(int i=0;i<1.0/h;i++)
        {
            for(int j=0;j<1.0/h;j++)printf("%lf ",FD1(i*h,j*h));
            puts("");
        }
        fclose(stdout);

        FD_regular FD2(Fun2,h,2);
        FD2.solve();
        freopen(("../data/fun2"+to_string(n)).c_str(),"w",stdout);
        for(int i=0;i<1.0/h;i++)
        {
            for(int j=0;j<1.0/h;j++)printf("%lf ",FD2(i*h,j*h));
            puts("");
        }
        fclose(stdout);

        FD_regular FD3(Fun3,h,3);
        FD3.solve();
        freopen(("../data/fun3"+to_string(n)).c_str(),"w",stdout);
        for(int i=0;i<1.0/h;i++)
        {
            for(int j=0;j<1.0/h;j++)printf("%lf ",FD3(i*h,j*h));
            puts("");
        }
        fclose(stdout);
    }


    freopen("../data/fun1_regu_points","w",stdout);
    for(int i=0;i<obj["n"].size();i++)
    {
        int n=obj["n"][i].asInt();
        double h=1.0/n;
        for(int k=1;k<=3;k++)
        {
            FD_regular FD1(Fun,h,k);
            FD1.solve();
            double avgerror=0;
            for(int j=0;j<obj["x"].size();j++)
            {
                double X=obj["x"][j].asDouble();
               double Y=obj["y"][j].asDouble();
            
                avgerror+=fabs(FD1(X,Y)-Fun(X,Y));
            }
            avgerror/=obj["x"].size();
            cout<<avgerror<<' ';
        }
        puts("");
        
    }   
    freopen("../data/fun2_regu_points","w",stdout);
    for(int i=0;i<obj["n"].size();i++)
    {
        int n=obj["n"][i].asInt();
        double h=1.0/n;
        for(int k=1;k<=3;k++)
        {
            FD_regular FD1(Fun2,h,k);
            FD1.solve();
            double avgerror=0;
            for(int j=0;j<obj["x"].size();j++)
            {
                double X=obj["x"][j].asDouble();
               double Y=obj["y"][j].asDouble();
            
                avgerror+=fabs(FD1(X,Y)-Fun2(X,Y));
            }
            avgerror/=obj["x"].size();
            cout<<avgerror<<' ';
        }
        puts("");
        
    }   
    freopen("../data/fun3_regu_points","w",stdout);
    for(int i=0;i<obj["n"].size();i++)
    {
        int n=obj["n"][i].asInt();
        double h=1.0/n;
        for(int k=1;k<=3;k++)
        {
            FD_regular FD1(Fun3,h,k);
            FD1.solve();
            double avgerror=0;
            for(int j=0;j<obj["x"].size();j++)
            {
                double X=obj["x"][j].asDouble();
               double Y=obj["y"][j].asDouble();
            
                avgerror+=fabs(FD1(X,Y)-Fun3(X,Y));
            }
            avgerror/=obj["x"].size();
            cout<<avgerror<<' ';
        }
        puts("");
        
    }   
    fclose(stdout);

    freopen("../data/fun1_irregu_points","w",stdout);
    for(int i=0;i<obj["n"].size();i++)
    {
        int n=obj["n"][i].asInt();
        double h=1.0/n;
        for(int k=1;k<=3;k++)
        if(k!=2)
        {
            FD_irregular FD1(Fun,k,h,C);
            FD1.solve();
            double avgerror=0;
            for(int j=0;j<obj["x"].size();j++)
            {
                double X=obj["x"][j].asDouble();
               double Y=obj["y"][j].asDouble();
            
                avgerror+=fabs(FD1((int)(X/h),(int)(Y/h))-Fun(X,Y));
            }
            avgerror/=obj["x"].size();
            cout<<avgerror<<' ';
        }
        puts("");
        
    }   

    freopen("../data/fun2_irregu_points","w",stdout);
    for(int i=0;i<obj["n"].size();i++)
    {
        int n=obj["n"][i].asInt();
        double h=1.0/n;
        for(int k=1;k<=3;k++)
        if(k!=2)
        {
            FD_irregular FD1(Fun2,k,h,C);
            FD1.solve();
            double avgerror=0;
            for(int j=0;j<obj["x"].size();j++)
            {
                double X=obj["x"][j].asDouble();
               double Y=obj["y"][j].asDouble();
            
                avgerror+=fabs(FD1((int)(X/h),(int)(Y/h))-Fun2(X,Y));
            }
            avgerror/=obj["x"].size();
            cout<<avgerror<<' ';
        }
        puts("");
        
    }   

    freopen("../data/fun3_irregu_points","w",stdout);
    for(int i=0;i<obj["n"].size();i++)
    {
        int n=obj["n"][i].asInt();
        double h=1.0/n;
        for(int k=1;k<=3;k++)
        if(k!=2)
        {
            FD_irregular FD1(Fun3,k,h,C);
            FD1.solve();
            double avgerror=0;
            for(int j=0;j<obj["x"].size();j++)
            {
                double X=obj["x"][j].asDouble();
               double Y=obj["y"][j].asDouble();
            
                avgerror+=fabs(FD1((int)(X/h),(int)(Y/h))-Fun3(X,Y));
            }
            avgerror/=obj["x"].size();
            cout<<avgerror<<' ';
        }
        puts("");
        
    }   

    freopen("../data/fun1_regu_error_norm","w",stdout);
    for(int i=0;i<obj["n"].size();i++)
    {
        int n=obj["n"][i].asInt();
        double h=1.0/n;
        
        FD_regular FD(Fun,h,1);
        FD.solve();
        cout<<FD.error_norm(1)<<' '<<FD.error_norm(2)<<' '<<FD.error_norm(-1)<<endl;
    }   
    fclose(stdout);
    freopen("../data/fun2_regu_error_norm","w",stdout);
    for(int i=0;i<obj["n"].size();i++)
    {
        int n=obj["n"][i].asInt();
        double h=1.0/n;
        
        FD_regular FD(Fun2,h,2);
        FD.solve();
        cout<<FD.error_norm(1)<<' '<<FD.error_norm(2)<<' '<<FD.error_norm(-1)<<endl;
    }   
    fclose(stdout);
    freopen("../data/fun3_regu_error_norm","w",stdout);
    for(int i=0;i<obj["n"].size();i++)
    {
        int n=obj["n"][i].asInt();
        double h=1.0/n;
        
        FD_regular FD(Fun3,h,3);
        FD.solve();
        cout<<FD.error_norm(1)<<' '<<FD.error_norm(2)<<' '<<FD.error_norm(-1)<<endl;
    }   
    fclose(stdout);

    freopen("../data/fun1_irregu_error_norm","w",stdout);
    for(int i=0;i<obj["n"].size();i++)
    {
        int n=obj["n"][i].asInt();
        double h=1.0/n;
        
        FD_irregular FD(Fun,3,h,C);
        FD.solve();
        cout<<FD.error_norm(1)<<' '<<FD.error_norm(2)<<' '<<FD.error_norm(-1)<<endl;
    }   
    fclose(stdout);
    freopen("../data/fun2_irregu_error_norm","w",stdout);
    for(int i=0;i<obj["n"].size();i++)
    {
        int n=obj["n"][i].asInt();
        double h=1.0/n;
        
        FD_irregular FD(Fun2,3,h,C);
        FD.solve();
        cout<<FD.error_norm(1)<<' '<<FD.error_norm(2)<<' '<<FD.error_norm(-1)<<endl;
    }   
    fclose(stdout);
    freopen("../data/fun3_irregu_error_norm","w",stdout);
    for(int i=0;i<obj["n"].size();i++)
    {
        int n=obj["n"][i].asInt();
        double h=1.0/n;
        
        FD_irregular FD(Fun3,3,h,C);
        FD.solve();
        cout<<FD.error_norm(1)<<' '<<FD.error_norm(2)<<' '<<FD.error_norm(-1)<<endl;
    }   
    fclose(stdout);
    //cerr<<FD.error_norm(1)<<' '<<FD.error_norm(2)<<' '<<FD.error_norm(-1)<<endl;
    return 0;
}