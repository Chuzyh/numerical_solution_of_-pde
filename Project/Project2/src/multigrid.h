#include "BVP.h"

enum bound_conditon { Dirichlet, Neumann,mixed};
enum restriction_operators {full_weighting,injection} ;
enum interpolation_operators {linear,quadratic};
enum cycles {V_cycle,FMG};
enum stopping_criteria{max_iteration,rela_accuracy};

vector<double> operator -(vector<double> a,vector<double> b)
{
    if(a.size()!=b.size()){cerr<<"wrong size"<<endl;exit(-1);}
    vector<double> re;re.resize(a.size());
    for(int i=0;i<(int)a.size();i++)re[i]=a[i]-b[i];
    return re;
}
vector<double> operator +(vector<double> a,vector<double> b)
{
    if(a.size()!=b.size()){cerr<<"wrong size"<<endl;exit(-1);}
    vector<double> re;re.resize(a.size());
    for(int i=0;i<(int)a.size();i++)re[i]=a[i]+b[i];
    return re;
}
template <int dim>
class Multigrid_Method
{
private:
    bound_conditon BC;
    restriction_operators RO;
    interpolation_operators IO;
    cycles CY;
    stopping_criteria SC;
    int max_iter;double esplion;
    double h;
    int N;
};
template<>
class Multigrid_Method<1>
{
    private:
        bound_conditon BC;
        restriction_operators RO;
        interpolation_operators IO;
        cycles CY;
        stopping_criteria SC;
        int mn_N;double esplion;
        double h,ori_h;
        int N;
        Function & f;
        vector<double> u;
    public:
        Multigrid_Method (Function &f,bound_conditon BC,restriction_operators RO,interpolation_operators IO,cycles CY,stopping_criteria SC,double h,vector<double> u,double st_parm):f(f),BC(BC),RO(RO),IO(IO),CY(CY),SC(SC),h(h),u(u)
        {
            cerr<<"init"<<endl;
            N=(int)(1.0/h+1);ori_h=h;if(SC==max_iteration)mn_N=max((int)(1.0/h/pow(2.0,st_parm)+1),3);else esplion=st_parm;
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
        double quad(double x0,double x1,double x2,double X)
        {
            double c=x0;
            double a=(x2-2*x1+x0)/2;
            double b=x1-a-c;
            return a*X*X+b*X+c;
        }
        vector<double> interpolation(vector<double> I)
        {
            vector<double> re;
            re.clear();
            if(IO==linear||re.size()==2)
            {
                for(int i=0;i<(int)I.size()-1;i++)
                {
                    re.push_back(I[i]);
                    re.push_back((I[i]+I[i+1])/2);
                }
            }
            else
            {
                for(int i=0;i<(int)I.size()-2;i++)
                {
                    re.push_back(I[i]);
                    re.push_back(quad(I[i],I[i+1],I[i+2],0.5));
                }
                re.push_back(I[(int)I.size()-2]);
                re.push_back(quad(I[(int)I.size()-3],I[(int)I.size()-2],I[(int)I.size()-1],1.5));
            } 
            re.push_back(I[(int)I.size()-1]);
            return re;
        }
        vector<double> jacobi(vector<double> u,vector<double> F)
        {
            vector<double> oldu(u);
            double omega=2.0/3;
            if(BC==Dirichlet)
            {
                u[0]=F[0];
                u[u.size()-1]=F[u.size()-1];
                for(int i=1;i<(int)u.size()-1;i++)u[i]=(oldu[i-1]+oldu[i+1])*omega/2+(1.0-omega)*oldu[i]+omega*h*h/2*F[i];
            }else 
            {
                u[0]=(-F[0]*h/1.5+2.0/1.5*oldu[1]-1.0/3*oldu[2])*omega+oldu[0]*(1-omega);
                u[u.size()-1]=F[u.size()-1];
                for(int i=1;i<(int)u.size()-1;i++)u[i]=(oldu[i-1]+oldu[i+1])*omega/2+(1.0-omega)*oldu[i]+omega*h*h/2*F[i];
            }
            return u;
        }
        double accuracy(vector<double> u)
        {
            double re=0;
            for(int i=1;i<(int)u.size()-1;i++)
            {
                re+=fabs(f(h*i)-u[i]);
            }
            return re*h;
        }
        int can_stop(vector<double> u)
        {
            if(SC==max_iteration)return u.size()<=mn_N;
            else 
            {
                return accuracy(u)<esplion;
            }
        }
        void out()
        {
            for(int i=0;i<(int)u.size();i++)
                printf("%lf ",u[i]);
            puts("");
            for(int i=0;i<(int)u.size();i++)
                printf("%lf ",f(h*i));
            puts("");
        }
        void solve(int v1,int v2)
        {
            vector<double> F;F.clear();
            for(int i=0;i<(int)u.size();i++)F.push_back(f.laplace(h*i));
            if(BC==Dirichlet)F[0]=f(0),F[u.size()-1]=f(1);
            else F[0]=f.diff(0),F[u.size()-1]=f(1);
            if(CY==V_cycle)u=Vcycle(v1,v2,u,F);else u=FMG(v1,v2,F);
        }
        vector<double> laplace(vector<double> u)
        {
            vector<double> re;re.resize(u.size());
            if(BC==Dirichlet)re[0]=u[0],re[u.size()-1]=u[u.size()-1];
            else re[0]=-1.5*u[0]/h+2.0*u[1]/h-0.5*u[2]/h,re[u.size()-1]=u[u.size()-1];
            for(int i=1;i<(int)u.size()-1;i++)re[i]=u[i]*2/h/h-(u[i-1]+u[i+1])/h/h;
            return re;
        }
        vector<double> Vcycle(int v1,int v2,vector<double> u,vector<double> F)
        {   
            for(int i=1;i<=v1;i++)
            {
                u=jacobi(u,F);
            }
                
            if(!can_stop(u))
            {
                vector<double> u2;u2.resize(u.size()/2+1);
                for(int i=0;i<u2.size();i++)u2[i]=0;
                vector<double> F2=F;
                
                F2=restriction(F-laplace(u));
                h=h*2;
                u2=Vcycle(v1,v2,u2,F2);
                u2=interpolation(u2);
                u=u+u2;
                h/=2;
            }
            for(int i=1;i<=v2;i++)u=jacobi(u,F);
            return u;
        }
        vector<double> FMG(int v1,int v2,vector<double> F)
        {   
            if(!can_stop(u))
            {
                vector<double> F2=F;
                
                F2=restriction(F-laplace(u));
                h=h*2;
                vector<double> U(u);
                u=restriction(u);
                u=interpolation(FMG(v1,v2,F2))+U;
                h/=2;
                
            }
            u=Vcycle(v1,v2,u,F);
            return u;
        }
        double error_norm(double p)
        {
            double re=0;
            double ave=0;
            for(int i=0;i<u.size();i++)ave+=f(h*i);ave/=u.size();
            for(int i=0;i<u.size();i++)re+=h*pow(fabs(u[i]-f(h*i)),p);
            return pow(re,1.0/p);
        }
        
};