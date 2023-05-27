#ifndef GaussLegendreRK
#define GaussLegendreRK
#include "IVPsolver.h"
class GaussLegendreRK_solver: public RKM
{
private:
    int n_steps,p;
    double Ti;
    vector<double> init;
    vector<vector<double>> result;
public:
    GaussLegendreRK_solver() {}
    GaussLegendreRK_solver(int _p,double time, vector<double> init_val,int n) : p(_p),Ti(time), init(init_val),n_steps(n) {}
    void solve()
    {
        if(p==1)
        {
            RKweights={1};
            RKnodes={1.0/2};
            RKMatrix={{1.0/2}};
        }else if(p==2)
        {
            RKweights={1.0/2,1.0/2};
            RKnodes={(3.0-sqrt(3.0))/6,(3.0+sqrt(3.0))/6};
            RKMatrix={{1.0/4,(3.0-sqrt(3.0)*2)/12},{(3.0+sqrt(3.0)*2)/12,1.0/4}};
        }else if(p==3)
        {
            RKweights={5.0/18,4.0/9,5.0/18};
            RKnodes={(3.0-sqrt(3.0))/6,(3.0+sqrt(3.0))/6};
            RKMatrix={{5.0/36,2.0/9-sqrt(15.0)/15,5.0/36-sqrt(15.0)/30},
            {5.0/36+sqrt(15.0)/24,2.0/9,5.0/36-sqrt(15.0)/24},
            {5.0/36+sqrt(15.0)/30,2.0/9+sqrt(15.0)/15,5.0/36}};
        }else ERROR("not finished");
        

        double tick=Ti/n_steps;
        result.clear();
        result.push_back(init);
        auto now=init;
        if(p==1)
        {
            for(int i=1;i<=n_steps;i++)
            {
                auto y0=now;
                while(1)
                {
                    auto ne=f(now+tick*(RKMatrix[0][0]*y0),(i+RKnodes[0])*tick);
                    if(norm_inf(ne-y0)>EPS)y0=ne;else break;
                }
                now=now+tick*(RKweights[0]*y0);
                result.push_back(now);
            }
        }else if(p==2)
        {
            for(int i=1;i<=n_steps;i++)
            {
                auto y0=now;
                auto y1=now;
                
                while(1)
                {
                    auto ne0=f(now+tick*(RKMatrix[0][0]*y0+RKMatrix[0][1]*y1),(i+RKnodes[0])*tick);
                    auto ne1=f(now+tick*(RKMatrix[1][0]*y0+RKMatrix[1][1]*y1),(i+RKnodes[1])*tick);
                    
                    if(norm_inf(ne1-y1)+norm_inf(ne0-y0)>EPS)y0=ne0,y1=ne1;else break;
                }
                now=now+tick*(RKweights[0]*y0+RKweights[1]*y1);
                result.push_back(now);
            }
        }else if(p==3)
        {
           for(int i=1;i<=n_steps;i++)
            {
                auto y0=now;
                auto y1=now;
                auto y2=now;
                while(1)
                {
                    auto ne0=f(now+tick*(RKMatrix[0][0]*y0+RKMatrix[0][1]*y1+RKMatrix[0][2]*y2),(i+RKnodes[0])*tick);
                    auto ne1=f(now+tick*(RKMatrix[1][0]*y0+RKMatrix[1][1]*y1+RKMatrix[1][2]*y2),(i+RKnodes[1])*tick);
                    auto ne2=f(now+tick*(RKMatrix[2][0]*y0+RKMatrix[2][1]*y1+RKMatrix[2][2]*y2),(i+RKnodes[2])*tick);
                    if(norm_inf(ne1-y1)+norm_inf(ne0-y0)+norm_inf(ne2-y2)>EPS)y0=ne0,y1=ne1,y2=ne2;else break;
                }
                now=now+tick*(RKweights[0]*y0+RKweights[1]*y1+RKweights[2]*y2);
                result.push_back(now);
            } 
        }
        freopen("GaussLegendreRK.data","w",stdout);
        for(int i=0;i<result.size();i++)
            cout<<result[i][0]<<' '<<result[i][1]<<endl;
        fclose(stdout);
    }
};
#endif