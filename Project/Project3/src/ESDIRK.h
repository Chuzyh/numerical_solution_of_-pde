#ifndef RK
#define RK
#include "IVPsolver.h"
class ESDIRK_solver: public RKM
{
private:
    int n_steps;
    double Ti;
    vector<double> init;
    vector<vector<double>> result;
public:
    ESDIRK_solver() {}
    ESDIRK_solver(double time, vector<double> init_val,int n) : Ti(time), init(init_val),n_steps(n) {}
    void solve()
    {
        RKweights={82889.0/524892,0,15625.0/83664,69875.0/102672,-2260.0/8211,1.0/4};
        RKnodes={0,1.0/2,83.0/250,31.0/50,17.0/20,1.0};
        RKMatrix={{0},{1.0/4,1.0/4},
        {8611.0/62500,-1743.0/31250,1.0/4},
        {5012029.0/34652500,-654441.0/29225000,174375.0/388108,1.0/4},
        {15267082809.0/155376265600,-71443401.0/120774400,730878875.0/902184768,2285395.0/8070912,1.0/4},
        {82889.0/524892,0,15625.0/83664,69875.0/102672,-2260.0/8211,1.0/4}};

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
        freopen("ESDIRK.data","w",stdout);
        for(int i=0;i<result.size();i++)
            cout<<result[i][0]<<' '<<result[i][1]<<endl;
        fclose(stdout);
    }
};
#endif