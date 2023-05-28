#ifndef ESDIRK
#define ESDIRK
#include "IVPsolver.h"
class ESDIRK_solver: public RKM
{
public:
    ESDIRK_solver() {}
    static pIVPsolver createIVPsolver() { return pIVPsolver(new ESDIRK_solver()); }
    void solve()
    {
        RKweights={82889.0/524892,0,15625.0/83664,69875.0/102672,-2260.0/8211,1.0/4};
        RKnodes={0,1.0/2,83.0/250,31.0/50,17.0/20,1.0};
        RKMatrix={{0},
        {1.0/4,1.0/4},
        {8611.0/62500,-1743.0/31250,1.0/4},
        {5012029.0/34652500,-654441.0/2922500,174375.0/388108,1.0/4},
        {15267082809.0/155376265600,-71443401.0/120774400,730878875.0/902184768,2285395.0/8070912,1.0/4},
        {82889.0/524892,0,15625.0/83664,69875.0/102672,-2260.0/8211,1.0/4}};

        double tick=Ti/n_steps;
        result.clear();
        result.push_back(init);
        auto now=init;
        for(int i=1;i<=n_steps;i++)
        {
            auto y0=f(now,i*tick);
            auto y1=y0;
            while(1)
            {
                auto ne=f(now+tick*(RKMatrix[1][0]*y0+RKMatrix[1][1]*y1),(i+RKnodes[1])*tick);
                if(norm_inf(ne-y1)>EPS)y1=ne;else break;
            }
            auto y2=y1;
            while(1)
            {
                auto ne=f(now+tick*(RKMatrix[2][0]*y0+RKMatrix[2][1]*y1+RKMatrix[2][2]*y2),(i+RKnodes[2])*tick);
                if(norm_inf(ne-y2)>EPS)y2=ne;else break;
            }
            auto y3=y2;
            while(1)
            {
                auto ne=f(now+tick*(RKMatrix[3][0]*y0+RKMatrix[3][1]*y1+RKMatrix[3][2]*y2+RKMatrix[3][3]*y3),(i+RKnodes[3])*tick);
                if(norm_inf(ne-y3)>EPS)y3=ne;else break;
            }
            auto y4=y3;
            while(1)
            {
                auto ne=f(now+tick*(RKMatrix[4][0]*y0+RKMatrix[4][1]*y1+RKMatrix[4][2]*y2+RKMatrix[4][3]*y3+RKMatrix[4][4]*y4),(i+RKnodes[4])*tick);
                if(norm_inf(ne-y4)>EPS)y4=ne;else break;
            }
            auto y5=y4;
            while(1)
            {
                auto ne=f(now+tick*(RKMatrix[5][0]*y0+RKMatrix[5][1]*y1+RKMatrix[5][2]*y2+RKMatrix[5][3]*y3+RKMatrix[5][4]*y4+RKMatrix[5][5]*y5),(i+RKnodes[5])*tick);
                if(norm_inf(ne-y5)>EPS)y5=ne;else break;
            }
            now=now+tick*(RKMatrix[5][0]*y0+RKMatrix[5][1]*y1+RKMatrix[5][2]*y2+RKMatrix[5][3]*y3+RKMatrix[5][4]*y4+RKMatrix[5][5]*y5);
            result.push_back(now);
        }
        freopen("ESDIRK.data","w",stdout);
        for(int i=0;i<result.size();i++)
            cout<<result[i][0]<<' '<<result[i][1]<<endl;
        fclose(stdout);
    }
};
void register_ESDIRK()
{
    classFactory &F=classFactory::createFactory();
    F.registerProduct("ESDIRK",[](){return ESDIRK_solver::createIVPsolver();});
}
#endif