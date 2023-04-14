#include "BVP.h"

enum bound_conditon { Dirichlet, Neumann,mixed};
enum restriction_operators {full_weighting,injection} ;
enum interpolation_operators {linear,quadratic};
enum cycles {V_cycle,FMG};
enum stopping_criteria{max_iteration,rela_accuracy};
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
        int max_iter,iter_time;double esplion;
        double h,ori_h;
        int N;
        Function & f;
        vector<double> u;
    public:
        Multigrid_Method (Function &f,bound_conditon BC,restriction_operators RO,interpolation_operators IO,cycles CY,stopping_criteria SC,double h,vector<double> u,double st_parm):f(f),BC(BC),RO(RO),IO(IO),CY(CY),SC(SC),h(h),u(u)
        {
            cerr<<"init"<<endl;iter_time=0;
            N=(int)(1.0/h+1);ori_h=h;if(SC==max_iteration)max_iter=(int)st_parm;else esplion=st_parm;
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
            if(IO==linear)
            {
                for(int i=0;i<(int)I.size()-1;i++)
                {
                    re.push_back(I[i]);
                    re.push_back((I[i]+I[i+1])/2);
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
                oldu[0]=u[0]=f(0);
                oldu[u.size()-1]=u[u.size()-1]=f(1);
                
                for(int i=1;i<(int)u.size()-1;i++)u[i]=h*h/2*(f.laplace(h*i))+(oldu[i-1]+oldu[i+1])/2;//u[i]=-((oldu[i]*2-oldu[i-1]-oldu[i+1])-f.laplace(h*i)*h*h/2);
            }else
            {
                {cerr<<"Neumann not finished"<<endl;exit(-1);}
            }
            for(int i=1;i<(int)u.size()-1;i++)u[i]=u[i]*1/3+oldu[i]*2/3;
        }
        double accuracy()
        {
            double re=0;
            for(int i=1;i<(int)u.size();i++)
            {
                re+=fabs(f(h*i)-u[i]);
            }
            return re*h;
        }
        int can_stop()
        {
            if(SC==max_iteration)return iter_time>=max_iter;
            else 
            {
                return accuracy()<esplion;
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
        void Vcycle(int v1,int v2)
        {   
            for(int i=1;i<=v1;i++)
                jacobi();
            if(!can_stop())
            {
                h=h*2;
                u=restriction(u);iter_time++;
                Vcycle(v1,v2);
                u=interpolation(u);
                h/=2;
            }else
            {
                int M=u.size();
                SparseMatrix<double> A(M,M);
                vector<Triplet<double> > a;
                MatrixXd b(M,1);
                b=MatrixXd::Zero(M,1);
                a.push_back(Triplet<double>(0,0,1));
                a.push_back(Triplet<double>(M-1,M-1,1));
                b(0,0)=f(0);
                b(M-1,0)=f(1);
                for(int i=1;i<M-1;i++)
                {
                    double h2=h*h;
                    a.push_back(Triplet<double>(i,i,2/h2));
                    a.push_back(Triplet<double>(i,i+1,-1/h2));
                    a.push_back(Triplet<double>(i,i-1,-1/h2));
                    b(i,0)=f.laplace(h*i);
                }
                A.setFromTriplets(a.begin(),a.end());
            
                A.makeCompressed();
                SparseLU<SparseMatrix<double> > Solver;
                Solver.compute(A);
                VectorXd X=Solver.solve(b);
                for(int i=0;i<M;i++)u[i]=X(i);
            }
            for(int i=1;i<=v2;i++)jacobi();
            
        }
        double error_norm(double p)
        {
            double re=0;
            for(int i=0;i<u.size()-1;i++)re+=h*pow(fabs(u[i]-f(h*i)),p);
            return pow(re,1.0/p);
        }
        
};