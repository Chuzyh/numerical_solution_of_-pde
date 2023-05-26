#ifndef ABM
#define ABM
#include "IVPsolver.h"
class AdamsBashforth_solver: public LMM
{
private:
    int n_steps,p;
    double Ti;
    vector<double> init;
    vector<vector<double>> result;
public:
    AdamsBashforth_solver() {}
    AdamsBashforth_solver(int _p, double time, vector<double> init_val,int n) :p(_p), Ti(time), init(init_val),n_steps(n) {}
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
        if(p==1)
        {
            beta={1};
        }else if(p==2)
        {
            beta={3.0/2,-1.0/2};
        }else if(p==3)
        {
            beta={23.0/12,-16.0/12,5.0/12};
        }else if(p==4)
        {
            beta={55.0/24,-59.0/24,37.0/24,-9.0/24};
        }
        for(int i=p;i<=n_steps;i++)
        {
            now=result[i-1];
            for(int j=0;j<p;j++)
                now=now+tick*(beta[j]*f(result[i-j-1],(i-j-1)*tick));
            result.push_back(now);
        }
        freopen("ABF.data","w",stdout);
        for(int i=0;i<result.size();i++)
        {
            cout<<result[i][0]<<' '<<result[i][1]<<endl;
        }
        fclose(stdout);
    }
};
#endif