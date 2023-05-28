#ifndef RK
#define RK
#include "IVPsolver.h"
class classicalRK_solver: public RKM
{

public:
    classicalRK_solver() {}
    static pIVPsolver createIVPsolver() { return pIVPsolver(new classicalRK_solver()); }
    
    void solve()
    {
        RKweights={1.0/6,2.0/6,2.0/6,1.0/6};
        RKnodes={0,1.0/2,1.0/2,1.0};
        RKMatrix={{0},{1.0/2,0,0},{0,1.0/2,0},{0,0,1,0}};

        double tick=Ti/n_steps;
        result.clear();
        result.push_back(init);
        auto now=init;
        for(int i=1;i<=n_steps;i++)
        {
            auto y1=f(now,i*tick);
            auto y2=f(now+0.5*tick*y1,(i+0.5)*tick);
            auto y3=f(now+0.5*tick*y2,(i+0.5)*tick);
            auto y4=f(now+tick*y3,(i+1)*tick);
            now=now+(tick/6)*(y1+2.0*y2+2.0*y3+y4);
            result.push_back(now);
        }
        freopen("RK.data","w",stdout);
        for(int i=0;i<result.size();i++)
            cout<<result[i][0]<<' '<<result[i][1]<<endl;
        fclose(stdout);
    }
    void get_exact()
    {
        n_steps=1e5;
        double tick=Ti/n_steps;
        result.clear();
        result.push_back(init);
        auto now=init;
        for(int i=1;i<=n_steps;i++)
        {
            auto y1=f(now,i*tick);
            auto y2=f(now+0.5*tick*y1,(i+0.5)*tick);
            auto y3=f(now+0.5*tick*y2,(i+0.5)*tick);
            auto y4=f(now+tick*y3,(i+1)*tick);
            auto u1=now+(tick/6)*(y1+2.0*y2+2.0*y3+y4);
            
            tick/=2;
            y1=f(now,i*tick);
            y2=f(now+0.5*tick*y1,(i+0.5)*tick);
            y3=f(now+0.5*tick*y2,(i+0.5)*tick);
            y4=f(now+tick*y3,(i+1)*tick);
            auto u1_2=now+(tick/6)*(y1+2.0*y2+2.0*y3+y4);
            y1=f(u1_2,i*tick);
            y2=f(u1_2+0.5*tick*y1,(i+0.5)*tick);
            y3=f(u1_2+0.5*tick*y2,(i+0.5)*tick);
            y4=f(u1_2+tick*y3,(i+1)*tick);
            auto _u1=u1_2+(tick/6)*(y1+2.0*y2+2.0*y3+y4);
            
            now=_u1+1.0/(pow(2.0,4)-1)*(_u1-u1);
            result.push_back(now);
            tick*=2;
        }
        freopen("RK.data","w",stdout);
        for(int i=0;i<result.size();i++)
            cout<<result[i][0]<<' '<<result[i][1]<<endl;
        fclose(stdout);
    }
};void register_CRK()
    {
        classFactory &F=classFactory::createFactory();
        F.registerProduct("classicalRK",[](){return classicalRK_solver::createIVPsolver();});
    }
#endif