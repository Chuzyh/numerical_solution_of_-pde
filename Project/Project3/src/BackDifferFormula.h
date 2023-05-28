#ifndef BDF
#define BDF
#include "IVPsolver.h"
class BackDifferFormula_solver: public LMM
{
public:
    BackDifferFormula_solver() {}
    static pIVPsolver createIVPsolver() { return pIVPsolver(new BackDifferFormula_solver()); }
    
    void solve()
    {
        double tick=Ti/n_steps;
        result.clear();
        result.push_back(init);
        auto now=init;
        for(int i=1;i<p;i++)
        {
            auto y1=f(now,i*tick);
            auto y2=f(now+0.5*tick*y1,(i+0.5)*tick);
            auto y3=f(now+0.5*tick*y2,(i+0.5)*tick);
            auto y4=f(now+tick*y3,(i+1)*tick);
            now=now+(tick/6)*(y1+2.0*y2+2.0*y3+y4);
            result.push_back(now);
        }
        vector<double> beta2;
        if(p==1)
        {
            beta2={1.0};
            alpha={-1.0};
            beta={1.0};
        }else if(p==2)
        {
            beta2={3.0/2,-1.0/2};
            alpha={-4.0/3,1.0/3};
            beta={2.0/3};
        }else if(p==3)
        {
            beta2={23.0/12,-16.0/12,5.0/12};
            alpha={-18.0/11,9.0/11,-2.0/11};
            beta={6.0/11};
        }else if(p==4)
        {
            beta2={55.0/24,-59.0/24,37.0/24,-9.0/24};
            alpha={-48.0/25,36.0/25,-16.0/25,3.0/25};
            beta={12.0/25};
        }

        for(int i=p;i<=n_steps;i++)
        {
            now=result[i-1];
            for(int j=0;j<p;j++)
                now=now+tick*(beta2[j]*f(result[i-j-1],(i-j-1)*tick));
            result.push_back(now);

            now=now-now;
            for(int j=0;j<p;j++)
                now=now-alpha[j]*result[i-j-1];
            now=now+tick*beta[0]*f(result[i],i*tick);
            result[i]=now;
        }
        freopen("BDF.data","w",stdout);
        for(int i=0;i<result.size();i++)
        {
            cout<<result[i][0]<<' '<<result[i][1]<<endl;
        }
        fclose(stdout);
    }
};
void register_BDF()
{
    classFactory &F=classFactory::createFactory();
    F.registerProduct("BackDifferFormula",[](){return BackDifferFormula_solver::createIVPsolver();});
}
#endif