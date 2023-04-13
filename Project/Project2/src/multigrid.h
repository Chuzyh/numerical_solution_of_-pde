#include "BVP.h"


template <int dim>
class Multigrid_Method
{
private:
    enum bound_conditon { Dirichlet, Neumann,mixed} BC;
    enum restriction_operators {full_weighting,injection} RO;
    enum interpolation_operators {linear,quadratic} IO;
    enum cycles {V_cycle,FMG} CY;
    enum stopping_criteria{max_iteration,rela_accuracy} SC;
    int max_iter;double esplion;
    double h;
    int N;
};
template<>
class Multigrid_Method<1>
{
    private:
        enum bound_conditon { Dirichlet, Neumann,mixed} BC;
        enum restriction_operators {full_weighting,injection} RO;
        enum interpolation_operators {linear,quadratic} IO;
        enum cycles {V_cycle,FMG} CY;
        enum stopping_criteria{max_iteration,rela_accuracy} SC;
        int max_iter;double esplion;
        double h,ori_h;
        int N;
        Function & f;
        vector<double> u;
    public:
        Multigrid_Method (Function &f,bound_conditon BC,restriction_operators RO,interpolation_operators IO,cycles CY,stopping_criteria SC,double h,vector<double> u):f(f),BC(BC),RO(RO),IO(IO),CY(CY),SC(SC),h(h),u(u)
        {
            N=(int)(1.0/h+1);ori_h=h;
            if(u.size()!=N){cerr<<"wrong size of initial guess"<<endl;exit(-1);}
        }
        vector<double> restriction(vector<double> I)
        {
            vector<double> re;
            re.clear();
            re.push_back(I[0]);
            if(RO==injection)
            {
                for(int i=2;i<(int)I.size()-1;i+=2)re.push_back(I[i]);
            }
            else
            {
                for(int i=2;i<(int)I.size()-1;i+=2)re.push_back((I[i]*2+I[i-1]+I[i+1])/4);
            } 
            re.push_back(I[(int)I.size()-1]);
            return re;
        }
        vector<double> interpolation(vector<double> I)
        {
            vector<double> re;
            re.clear();
            re.push_back(I[0]);
            if(IO==linear)
            {
                for(int i=1;i<(int)I.size()-1;i++)
                {
                    re.push_back((I[i-1]+I[i])/2);
                    re.push_back(I[i]);
                }
            }
            else
            {
                {cerr<<"quadratic not finished"<<endl;exit(-1);}
            } 
            re.push_back(I[(int)I.size()-1]);
            return re;
        }
        void jacobi()
        {
            vector<double> oldu(u);
            //I[-1]=I[0]-(I[1]-I[0])
            if(BC==Dirichlet)
            {
                for(int i=1;i<(int)u.size()-1;i++)u[i]=(oldu[i]*2-oldu[i-1]-oldu[i+1])-f(h*i)*h;
            }else
            {
                {cerr<<"Neumann not finished"<<endl;exit(-1);}
            }
            for(int i=0;i<(int)u.size();i++)u[i]=u[i]*2/3+oldu[i];
        }
        int can_stop()
        {

        }
        void Vcycle(int v1,int v2)
        {   
            for(int i=1;i<=v1;i++)jacobi();
            if(!can_stop())
            {
                h=h*2;
                u=restriction(u);
                Vcycle(v1,v2);
                u=interpolation(u);
            }
            for(int i=1;i<=v2;i++)jacobi();
        }
};