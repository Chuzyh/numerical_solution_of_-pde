#ifndef FehlbergRK
#define FehlbergRK
#include "IVPsolver.h"
class FehlbergRK_solver: public RKM
{
private:
    int n_steps,p;
    double Ti;
    vector<double> init;
    vector<vector<double>> result;
public:
    FehlbergRK_solver() {}
    FehlbergRK_solver(int p,double time, vector<double> init_val,int n) : p(p),Ti(time), init(init_val),n_steps(n) {}
    void solve()
    {
        if(p==4)RKweights={25.0/216,0,1408.0/2565,2197.0/4104,-1.0/5,0};
        else RKweights={16.0/135,0,6656.0/12825,28561.0/56430,-9.0/50,2.0/55};
        RKnodes={0,1.0/4,3.0/8,12.0/13,1.0,1.0/2};
        RKMatrix={{0},
        {1.0/4},
        {3.0/32,9.0/32},
        {1932.0/2197,-7200.0/2197,7296.0/2197},
        {439.0/216,-8,3680.0/513,-845.0/4104},
        {-8.0/27,2,-3544.0/2565,1859.0/4104,-11.0/40}};

        double tick=Ti/n_steps;
        result.clear();
        result.push_back(init);
        auto now=init;
        
        for(int i=1;i<=n_steps;i++)
        {
            auto y0=f(now,i*tick);
            auto y1=f(now+tick*(RKMatrix[1][0]*y0),(i+RKnodes[1])*tick);
            auto y2=f(now+tick*(RKMatrix[2][0]*y0+RKMatrix[2][1]*y1),(i+RKnodes[2])*tick);
            auto y3=f(now+tick*(RKMatrix[3][0]*y0+RKMatrix[3][1]*y1+RKMatrix[3][2]*y2),(i+RKnodes[3])*tick);
            auto y4=f(now+tick*(RKMatrix[4][0]*y0+RKMatrix[4][1]*y1+RKMatrix[4][2]*y2+RKMatrix[4][3]*y3),(i+RKnodes[4])*tick);
            if(p==5)
            {
                auto y5=f(now+tick*(RKMatrix[5][0]*y0+RKMatrix[5][1]*y1+RKMatrix[5][2]*y2+RKMatrix[5][3]*y3+RKMatrix[5][4]*y4),(i+RKnodes[5])*tick);
                now=now+tick*RKweights[5]*y5;
            }
            now=now+tick*(RKweights[0]*y0+RKweights[1]*y1+RKweights[2]*y2+RKweights[3]*y3+RKweights[4]*y4);
            
            result.push_back(now);
        }
        freopen("FehlbergRK.data","w",stdout);
        for(int i=0;i<result.size();i++)
            cout<<result[i][0]<<' '<<result[i][1]<<endl;
        fclose(stdout);
    }
};
#endif