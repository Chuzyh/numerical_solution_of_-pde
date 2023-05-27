#ifndef DormandPrinceRK
#define DormandPrinceRK
#include "IVPsolver.h"
class DormandPrinceRK_solver: public RKM
{
private:
    int n_steps;
    double Ti=-1,rho_max=2.0,rho_min=0.5,rho=0.8,q=4.0;
    double e_rel=1e-4,e_abs=1e-8;
    vector<double> init;
    vector<vector<double>> result;
public:
    DormandPrinceRK_solver() {}
    DormandPrinceRK_solver(double time, vector<double> init_val,int n) : Ti(time), init(init_val),n_steps(n) {}
    DormandPrinceRK_solver(vector<double> init_val,int n) : init(init_val),n_steps(n) {}
    void solve()
    {
        RKweights={5179.0/57600,0,7571.0/16695,393.0/640,-92097.0/339200,187.0/2100,1.0/40};
        RKnodes={0,1.0/5,3.0/10,4.0/5,8.0/9,1.0,1.0};
        RKMatrix={{0},
        {1.0/5},
        {3.0/40,9.0/40},
        {44.0/45,-56.0/15,32.0/9},
        {19372.0/6561,-25360.0/2187,64448.0/6561,-212.0/729},
        {9017.0/3168,-355.0/33,46732.0/5247,49.0/176,-5103.0/18656},
        {35.0/384,0,500.0/1113,125.0/192,-2187.0/6784,11.0/84}};

        result.clear();
        result.push_back(init);
        auto now=init;
        
        double tick=Ti>0?Ti/n_steps:0.0001;
        
        for(int i=1;i<=n_steps;i++)
        
        {
            double E_ind=1e9;
            while(1)
            {
                auto y0=f(now,i*tick);
                auto y1=f(now+tick*(RKMatrix[1][0]*y0),(i+RKnodes[1])*tick);
                auto y2=f(now+tick*(RKMatrix[2][0]*y0+RKMatrix[2][1]*y1),(i+RKnodes[2])*tick);
                auto y3=f(now+tick*(RKMatrix[3][0]*y0+RKMatrix[3][1]*y1+RKMatrix[3][2]*y2),(i+RKnodes[3])*tick);
                auto y4=f(now+tick*(RKMatrix[4][0]*y0+RKMatrix[4][1]*y1+RKMatrix[4][2]*y2+RKMatrix[4][3]*y3),(i+RKnodes[4])*tick);
                auto y5=f(now+tick*(RKMatrix[5][0]*y0+RKMatrix[5][1]*y1+RKMatrix[5][2]*y2+RKMatrix[5][3]*y3+RKMatrix[5][4]*y4),(i+RKnodes[5])*tick);
            
                auto ne=now+tick*(RKMatrix[6][0]*y0+RKMatrix[6][1]*y1+RKMatrix[6][2]*y2+RKMatrix[6][3]*y3+RKMatrix[6][4]*y4+RKMatrix[6][5]*y5);
                auto y6=f(ne,(i+RKnodes[6])*tick);
            
                auto ne2=now+tick*(RKweights[0]*y0+RKweights[1]*y1+RKweights[2]*y2+RKweights[3]*y3+RKweights[4]*y4+RKweights[5]*y5+RKweights[6]*y6);

                E_ind=error_norm_inf(ne2,ne,e_abs,e_rel);
                
                if(Ti<0)
                {
                    tick=tick*min(rho_max,max(rho_min,rho*pow(E_ind,-1.0/(q+1))));
                    if(E_ind<=1)
                    {
                        now=ne;
                        break;
                    }
                }else
                {
                    now=ne;
                    break;
                }
                
            }
            
            
            result.push_back(now);
        }
        freopen("DormandPrinceRK.data","w",stdout);
        for(int i=0;i<result.size();i++)
            cout<<result[i][0]<<' '<<result[i][1]<<endl;
        fclose(stdout);
    }
};
#endif