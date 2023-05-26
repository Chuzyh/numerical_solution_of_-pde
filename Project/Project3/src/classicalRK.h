#ifndef RK
#define RK
#include "IVPsolver.h"
class classicalRK_solver: public RKM
{
private:
    int n_steps;
    double Ti;
    vector<double> init;
    vector<vector<double>> result;
public:
    classicalRK_solver() {}
    classicalRK_solver(double time, vector<double> init_val,int n) : Ti(time), init(init_val),n_steps(n) {}
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
};
#endif