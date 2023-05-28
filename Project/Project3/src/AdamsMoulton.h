#ifndef ADM
#define ADM
#include "IVPsolver.h"
class AdamsMoulton_solver: public LMM
{
public:
    AdamsMoulton_solver() {}
    static pIVPsolver createIVPsolver() { return pIVPsolver(new AdamsMoulton_solver()); }
    
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
        if(p==2)
        {
            beta2={1};
            beta={1.0/2,1.0/2};
        }else if(p==3)
        {
            beta2={3.0/2,-1.0/2};
            beta={5.0/12,8.0/12,-1.0/12};
        }else if(p==4)
        {
            beta2={23.0/12,-16.0/12,5.0/12};
            beta={9.0/24,19.0/24,-5.0/24,1.0/24};
        }else if(p==5)
        {
            beta2={55.0/24,-59.0/24,37.0/24,-9.0/24};
            beta={251.0/720,646.0/720,-264.0/720,106.0/720,-19.0/720};
        }

        for(int i=p;i<=n_steps;i++)
        {
            now=result[i-1];
            for(int j=0;j<p-1;j++)
                now=now+tick*(beta2[j]*f(result[i-j-1],(i-j-1)*tick));
            result.push_back(now);

            now=result[i-1];
            
            for(int j=0;j<p;j++)
                now=now+tick*(beta[j]*f(result[i-j],(i-j)*tick));
            
            result[i]=now;
        }
        freopen("ADM.data","w",stdout);
        for(int i=0;i<result.size();i++)
        {
            cout<<result[i][0]<<' '<<result[i][1]<<endl;
        }
        fclose(stdout);
    }
};
void register_ADM()
{
    classFactory &F=classFactory::createFactory();
    F.registerProduct("AdamsMoulton",[](){return AdamsMoulton_solver::createIVPsolver();});
}
#endif