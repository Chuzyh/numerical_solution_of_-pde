#include "BVP.h"

enum bound_conditon { Dirichlet, Neumann,mixed};
enum restriction_operators {full_weighting,injection} ;
enum interpolation_operators {linear,quadratic};
enum cycles {V_cycle,FMG};
enum stopping_criteria{max_iteration,rela_accuracy};
enum bound_conditon bound_conditon_string(string x)
{
    if(x=="Dirichlet")return bound_conditon(Dirichlet);
    if(x=="Neumann")return bound_conditon(Neumann);
    if(x=="mixed")return bound_conditon(mixed);
    cerr<<"wrong bound_conditon"<<endl;
    exit(-1);
}
enum restriction_operators restriction_operators_string(string x)
{
    if(x=="full_weighting")return restriction_operators(full_weighting);
    if(x=="injection")return restriction_operators(injection);
    cerr<<"wrong restriction_operators"<<endl;
    exit(-1);
}
enum interpolation_operators interpolation_operators_string(string x)
{
    if(x=="linear")return interpolation_operators(linear);
    if(x=="quadratic")return interpolation_operators(quadratic);
    cerr<<"wrong interpolation_operators"<<endl;
    exit(-1);
}
enum cycles cycles_string(string x)
{
    if(x=="V_cycle")return cycles(V_cycle);
    if(x=="FMG")return cycles(FMG);
    cerr<<"wrong cycles"<<x<<endl;
    exit(-1);
}

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
vector<vector<double> > operator -(vector<vector<double> > a,vector<vector<double> > b)
{
    if(a.size()!=b.size()){cerr<<"wrong size"<<endl;exit(-1);}
    vector<vector<double> > re;re.resize(a.size());
    for(int i=0;i<(int)a.size();i++)
    {
        if(a[i].size()!=b[i].size()){cerr<<"wrong size"<<endl;exit(-1);}
        re[i]=a[i]-b[i];
    }
    return re;
}
vector<vector<double> > operator +(vector<vector<double> > a,vector<vector<double> > b)
{
    if(a.size()!=b.size()){cerr<<"wrong size"<<endl;exit(-1);}
    vector<vector<double> > re;re.resize(a.size());
    for(int i=0;i<(int)a.size();i++)
    {
        if(a[i].size()!=b[i].size()){cerr<<"wrong size"<<endl;exit(-1);}
        re[i]=a[i]+b[i];
    }
    return re;
}

vector<vector<double> > operator /(vector<vector<double> > a,double b)
{
    vector<vector<double> > re;re.resize(a.size());
    for(int i=0;i<(int)a.size();i++)
    {
        for(int j=0;j<(int)a[i].size();j++)re[i].push_back(a[i][j]/b);
    }
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
        int mn_N;double st_parm;
        double h,ori_h;
        int N;
        Function & f;
        vector<double> u;
        vector<double> laplace(vector<double> u)
        {
            vector<double> re;re.resize(u.size());
            
            for(int i=1;i<(int)u.size()-1;i++)re[i]=u[i]*2/h/h-(u[i-1]+u[i+1])/h/h;
            if(BC==Dirichlet)re[0]=u[0],re[u.size()-1]=u[u.size()-1];
            else if(BC==Neumann)re[0]=-1.5*u[0]/h+2.0*u[1]/h-0.5*u[2]/h,re[u.size()-1]=1.5*u[u.size()-1]/h-2.0*u[u.size()-2]/h+0.5*u[u.size()-3]/h,re[u.size()/2]=u[u.size()/2];
            else re[0]=-1.5*u[0]/h+2.0*u[1]/h-0.5*u[2]/h,re[u.size()-1]=u[u.size()-1];
            return re;
        }
        vector<double> Vcycle(int v1,int v2,vector<double> u,vector<double> F)
        {   
            for(int i=1;i<=v1;i++)
            {
                u=jacobi(u,F);
            }
                
            if(u.size()>mn_N)
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
            if(u.size()>mn_N)
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
                for(int i=2;i<(int)I.size()-1;i+=2)
                {
                    if(BC==Neumann&&i==I.size()/2)re.push_back(I[i]);else re.push_back((I[i]*2+I[i-1]+I[i+1])/4);
                }
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
            }else if(BC==Neumann)
            {
                //cout<<u.size()<<' '<<N<<endl;
                u[0]=(-F[0]*h/1.5+2.0/1.5*oldu[1]-1.0/3*oldu[2])*omega+oldu[0]*(1-omega);
                u[u.size()-1]=(F[u.size()-1]*h/1.5+2.0/1.5*oldu[u.size()-2]-1.0/3*oldu[u.size()-3])*omega+oldu[u.size()-1]*(1-omega);
                
                for(int i=1;i<(int)u.size()-1;i++)u[i]=(oldu[i-1]+oldu[i+1])*omega/2+(1.0-omega)*oldu[i]+omega*h*h/2*F[i];
                u[u.size()/2]=F[u.size()/2];
            }else 
            {
                u[0]=(-F[0]*h/1.5+2.0/1.5*oldu[1]-1.0/3*oldu[2])*omega+oldu[0]*(1-omega);
                u[u.size()-1]=F[u.size()-1];
                for(int i=1;i<(int)u.size()-1;i++)u[i]=(oldu[i-1]+oldu[i+1])*omega/2+(1.0-omega)*oldu[i]+omega*h*h/2*F[i];
            }
            return u;
        }
    public:
        Multigrid_Method (Function &f,bound_conditon BC,restriction_operators RO,interpolation_operators IO,cycles CY,stopping_criteria SC,double h,vector<double> u,double st_parm,int mn_n):f(f),BC(BC),RO(RO),IO(IO),CY(CY),SC(SC),h(h),u(u),st_parm(st_parm),mn_N(mn_n)
        {
            
            N=(int)(1.0/h+1);ori_h=h;mn_N=max(3,mn_N);
            if(u.size()!=N){cerr<<"wrong size of initial guess"<<endl;exit(-1);}
        }
        
        double accuracy(vector<double> u)
        {
            double re=0;
            for(int i=1;i<(int)u.size()-1;i++)
            {
                re+=fabs((f(h*i)-u[i])/f(h*i));
            }
            return re*h;
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
        void solve(int v1,int v2,int show_detail)
        {
            vector<double> F;F.clear();
            for(int i=0;i<(int)u.size();i++)F.push_back(f.laplace(h*i));

            if(BC==Dirichlet)F[0]=f(0),F[u.size()-1]=f(1);
            else if(BC==Neumann)F[0]=f.diff(0),F[u.size()-1]=f.diff(1),F[u.size()/2]=f(0.5);
            else F[0]=f.diff(0),F[u.size()-1]=f(1);

            double la_error=0,la_resi=0;
            
            if(SC==rela_accuracy)
            {
                for(int i=1;i<=20&&(residual_norm(1)>st_parm);i++)
                {
                    if(CY==V_cycle)u=Vcycle(v1,v2,u,F);else u=FMG(v1,v2,F);
                    double now_error=error_norm(1);
                    double now_resi=residual_norm(1);
                    if(show_detail)printf(" %d & %.3e & %.3lf & %.3e & %.3lf\\\\ \\hline \n",i,now_error,la_error/now_error,now_resi,la_resi/now_resi);
                    la_error=now_error;la_resi=now_resi;
                }
                    
            }else
            {
                for(int i=1;i<=st_parm;i++)
                {
                    if(CY==V_cycle)u=Vcycle(v1,v2,u,F);else u=FMG(v1,v2,F);
                    double now_error=error_norm(1);
                    double now_resi=residual_norm(1);
                    if(show_detail)printf(" %d & %.3e & %.3lf & %.3e & %.3lf\\\\ \\hline \n",i,now_error,la_error/now_error,now_resi,la_resi/now_resi);
                    la_error=now_error;la_resi=now_resi;
                }
                    
            }
            
        }
        
        double error_norm(double p)
        {
            double re=0;
            
            for(int i=0;i<u.size();i++)re+=h*pow(fabs(u[i]-f(h*i)),p);
            return pow(re,1.0/p);
        }
        double residual_norm(double p)
        {
            double re=0;
            vector<double> F;F.clear();
            for(int i=0;i<(int)u.size();i++)F.push_back(f.laplace(h*i));
            if(BC==Dirichlet)F[0]=f(0),F[u.size()-1]=f(1);else
            if(BC==Neumann)F[0]=f.diff(0),F[u.size()-1]=f.diff(1),F[u.size()/2]=f(0.5);
            else F[0]=f.diff(0),F[u.size()-1]=f(1);

            vector<double> res=F-laplace(u);
            for(int i=0;i<u.size();i++)re+=h*pow(fabs(res[i]),p);
            return pow(re,1.0/p);
        }
        
};
template<>
class Multigrid_Method<2>
{
    private:
        const int dir[8][2]={{0,1},{0,-1},{1,0},{-1,0},{1,1},{1,-1},{-1,1},{-1,-1}};

        bound_conditon BC;
        restriction_operators RO;
        interpolation_operators IO;
        cycles CY;
        stopping_criteria SC;
        int mn_N;double st_parm;
        double h,ori_h;
        int N;
        Function2d & f;
        vector<vector<double> > u;
        vector<vector<double> > laplace(vector<vector<double> > u)
        {
            vector<vector<double> > re;re.clear();
            re.resize(u.size());
            for(int i=0;i<(int)u.size();i++)re[i].resize(u.size());

            double h2=h*h;
            for(int i=1;i<(int)u.size()-1;i++)
                for(int j=1;j<(int)u.size()-1;j++)
                    re[i][j]=u[i][j]*4/h2-(u[i][j+1]+u[i][j-1]+u[i+1][j]+u[i-1][j])/h2;
            if(BC==Dirichlet)
            {
                for(int i=0;i<(int)u.size();i++)
                {
                    re[0][i]=u[0][i];
                    re[u.size()-1][i]=u[u.size()-1][i];
                    re[i][0]=u[i][0];
                    re[i][u.size()-1]=u[i][u.size()-1];
                }
            }
            else if(BC==Neumann)
            {
                for(int i=1;i<(int)u.size()-1;i++)
                {
                    re[0][i]=-1.5*u[0][i]/h+2.0*u[1][i]/h-0.5*u[2][i]/h;

                    re[u.size()-1][i]=1.5*u[u.size()-1][i]/h-2.0*u[u.size()-2][i]/h+0.5*u[u.size()-3][i]/h;

                    re[i][0]=-1.5*u[i][0]/h+2.0*u[i][1]/h-0.5*u[i][2]/h;

                    re[i][u.size()-1]=1.5*u[i][u.size()-1]/h-2.0*u[i][u.size()-2]/h+0.5*u[i][u.size()-3]/h;

                }
                re[0][u.size()/2]=u[0][u.size()/2];
            }else 
            {
                for(int i=1;i<(int)u.size()-1;i++)
                {
                    re[0][i]=u[0][i];
                    re[u.size()-1][i]=u[u.size()-1][i];

                    re[i][0]=-1.5*u[i][0]/h+2.0*u[i][1]/h-0.5*u[i][2]/h;

                    re[i][u.size()-1]=1.5*u[i][u.size()-1]/h-2.0*u[i][u.size()-2]/h+0.5*u[i][u.size()-3]/h;

                }
            }
            return re;
        }
        vector<vector<double> > Vcycle(int v1,int v2,vector<vector<double> > u,vector<vector<double> > F)
        {   
            for(int i=1;i<=v1;i++)
            {
                u=jacobi(u,F);
            }
                
            if(u.size()>mn_N)
            {
                vector<vector<double> > u2;u2.resize(u.size()/2+1);
                for(int i=0;i<u2.size();i++)
                {
                    u2[i].resize(u.size()/2+1);
                    for(int j=0;j<u2.size();j++)u2[i][j]=0;
                }
                vector<vector<double> > F2=F;
                
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
        vector<vector<double> > FMG(int v1,int v2,vector<vector<double> > F)
        {   
            if(u.size()>mn_N)
            {
                vector<vector<double> > F2=F;
                
                F2=restriction(F-laplace(u));
                h=h*2;
                vector<vector<double> > U(u);
                u=restriction(u);
                u=interpolation(FMG(v1,v2,F2))+U;
                h/=2;
                
            }
            u=Vcycle(v1,v2,u,F);
            return u;
        }
        vector<vector<double> > restriction(vector<vector<double> > I)
        {
            vector<vector<double> > re;
            int grid_N=I.size();
            re.clear();
            re.resize(grid_N/2+1);
            for(int i=0;i<(int)re.size();i++)re[i].resize(grid_N/2+1);

            if(RO==injection)
            {
                for(int i=0;i<grid_N;i+=2)
                    for(int j=0;j<grid_N;j+=2)
                        re[i/2][j/2]=I[i][j];
            }
            else
            {
                for(int i=0;i<grid_N;i+=2)
                    for(int j=0;j<grid_N;j+=2)
                    {
                        if(i==0||j==0||i==grid_N-1||j==grid_N-1)
                            re[i/2][j/2]=I[i][j];
                        else    
                            re[i/2][j/2]=(I[i][j]*4+I[i][j+1]+I[i][j-1]+I[i+1][j]+I[i-1][j])/8;
                        
                    }
                
            } 
            if(BC==Neumann)re[re.size()/2][re.size()/2]=I[re.size()][re.size()];
            return re;
        }
        int in_grid(int grid_N,int x,int y)
        {
            return !(x<0||y<0||x>=grid_N||y>=grid_N);
        }
        double quad(double x0,double x1,double x2,double X)
        {
            double c=x0;
            double a=(x2-2*x1+x0)/2;
            double b=x1-a-c;
            return a*X*X+b*X+c;
        }
        vector<vector<double> > interpolation(vector<vector<double> > I)
        {
            vector<vector<double> > re;
            int grid_N=I.size();
            re.clear();
            re.resize(grid_N*2-1);
            for(int i=0;i<(int)re.size();i++)re[i].resize(grid_N*2-1);
            for(int i=0;i<(int)re.size();i++)
                for(int j=0;j<(int)re.size();j++)re[i][j]=0;
            if(IO==linear||re.size()==2)
            {
                for(int i=0;i<grid_N;i++)
                    for(int j=0;j<grid_N;j++)
                    {
                        if(i==grid_N/2&&j==grid_N/2)continue;
                        re[i*2][j*2]+=I[i][j];
                        for(int k=0;k<4;k++)
                        if(in_grid(re.size(),i*2+dir[k][0],j*2+dir[k][1]))
                            re[i*2+dir[k][0]][j*2+dir[k][1]]+=I[i][j]*0.5; 
                        for(int k=4;k<8;k++)
                        if(in_grid(re.size(),i*2+dir[k][0],j*2+dir[k][1]))
                            re[i*2+dir[k][0]][j*2+dir[k][1]]+=I[i][j]*0.25;
                        
                    }
            }
            else
            {
                for(int i=0;i<grid_N;i++)
                {
                    for(int j=0;j<grid_N-2;j++)
                    {
                        re[i*2][j*2]=I[i][j];
                        re[i*2][j*2+1]=quad(I[i][j],I[i][j+1],I[i][j+2],0.5);
                    }
                    re[i*2][(grid_N-2)*2]=I[i][grid_N-2];
                    re[i*2][(grid_N-2)*2+1]=(quad(I[i][(int)I.size()-3],I[i][(int)I.size()-2],I[i][(int)I.size()-1],1.5));
                    
                    re[i*2][grid_N*2-2]=I[i][(int)I.size()-1];
                }
                for(int j=0;j<re.size();j++)
                {
                    for(int i=1;i<re.size()-3;i+=2)
                    {
                        re[i][j]=quad(re[i-1][j],re[i+1][j],re[i+3][j],0.5);
                    }
                    re[re.size()-2][j]=(quad(re[re.size()-5][j],re[re.size()-3][j],re[re.size()-1][j],1.5));
                    
                }
            } 
            return re;
        }
        vector<vector<double> > jacobi(vector<vector<double> > u,vector<vector<double> > F)
        {
            vector<vector<double> > oldu(u);
            double omega=2.0/3;
            if(BC==Dirichlet)
            {
                for(int i=0;i<(int)u.size();i++)
                {
                    u[0][i]=F[0][i];
                    u[u.size()-1][i]=F[u.size()-1][i];
                    u[i][0]=F[i][0];
                    u[i][u.size()-1]=F[i][u.size()-1];

                    
                }
                
            }else if(BC==Neumann)
            {
                
                for(int i=1;i<(int)u.size()-1;i++)
                {
                    u[0][i]=(-F[0][i]*h/1.5+2.0/1.5*oldu[1][i]-1.0/3*oldu[2][i])*omega+oldu[0][i]*(1-omega);

                    u[u.size()-1][i]=(F[u.size()-1][i]*h/1.5+2.0/1.5*oldu[u.size()-2][i]-1.0/3*oldu[u.size()-3][i])*omega+oldu[u.size()-1][i]*(1-omega);

                    u[i][0]=(-F[i][0]*h/1.5+2.0/1.5*oldu[i][1]-1.0/3*oldu[i][2])*omega+oldu[i][0]*(1-omega);

                    u[i][u.size()-1]=(F[i][u.size()-1]*h/1.5+2.0/1.5*oldu[i][u.size()-2]-1.0/3*oldu[i][u.size()-3])*omega+oldu[i][u.size()-1]*(1-omega);

                }
                
            }else 
            {
                for(int i=1;i<(int)u.size()-1;i++)
                {
                    u[0][i]=F[0][i];
                    u[u.size()-1][i]=F[u.size()-1][i];
                    u[i][0]=(-F[i][0]*h/1.5+2.0/1.5*oldu[i][1]-1.0/3*oldu[i][2])*omega+oldu[i][0]*(1-omega);

                    u[i][u.size()-1]=(F[i][u.size()-1]*h/1.5+2.0/1.5*oldu[i][u.size()-2]-1.0/3*oldu[i][u.size()-3])*omega+oldu[i][u.size()-1]*(1-omega);

                }
            }
            for(int i=1;i<(int)u.size()-1;i++)
                    for(int j=1;j<(int)u.size()-1;j++)
                        u[i][j]=(oldu[i-1][j]+oldu[i+1][j]+oldu[i][j-1]+oldu[i][j+1])*omega/4+(1.0-omega)*oldu[i][j]+omega*h*h/4*F[i][j];
            if(BC==Neumann)u[0][u.size()/2]=F[0][u.size()/2];
            return u;
        }
    public:
        Multigrid_Method (Function2d &f,bound_conditon BC,restriction_operators RO,interpolation_operators IO,cycles CY,stopping_criteria SC,double h,vector<vector<double> > u,double st_parm,int mn_n):f(f),BC(BC),RO(RO),IO(IO),CY(CY),SC(SC),h(h),u(u),st_parm(st_parm),mn_N(mn_n)
        {
            
            N=(int)(1.0/h+1);ori_h=h;mn_N=max(3,mn_N);
            if(u.size()!=N){cerr<<"wrong size of initial guess"<<endl;exit(-1);}
        }
        
        double accuracy(vector<vector<double> > u)
        {
            double re=0;
            for(int i=1;i<(int)u.size()-1;i++)
            for(int j=1;j<(int)u.size()-1;j++)
            {
                re+=fabs((f(h*i,h*j)-u[i][j])/f(h*i,h*j));
            }
            return re*h*h;
        }
        void out()
        {
            for(int i=0;i<(int)u.size();i++)
            {
                for(int j=0;j<(int)u.size();j++)
                    printf("%lf ",u[i][j]);
                puts("");
            }
            
            puts("");
             for(int i=0;i<(int)u.size();i++)
            {
                for(int j=0;j<(int)u.size();j++)
                    printf("%lf ",f(h*i,h*j));
                puts("");
            }
            puts("");
        }
        void solve(int v1,int v2,int show_detail)
        {
            vector<vector<double> > F;F.clear();
            F.resize(u.size());
            for(int i=0;i<(int)u.size();i++)
                for(int j=0;j<(int)u.size();j++)
                    F[i].push_back(f.laplace(h*i,h*j));

            if(BC==Dirichlet)
            {
                for(int i=0;i<(int)u.size();i++)
                {
                    F[0][i]=f(0,h*i);
                    F[u.size()-1][i]=f(1,h*i);
                    F[i][0]=f(h*i,0);
                    F[i][u.size()-1]=f(h*i,1);
                }
            }
           else if(BC==Neumann)
           {

                F[0][u.size()/2]=f(0,0.5);
                for(int i=1;i<(int)u.size()-1;i++)
                {
                    F[0][i]=f.diff_x(0,h*i);
                    F[u.size()-1][i]=f.diff_x(1,h*i);
                    F[i][0]=f.diff_y(h*i,0);
                    F[i][u.size()-1]=f.diff_y(h*i,1);
                }

                v1*=10,v2*=10;
            }else 
            {
                for(int i=1;i<(int)u.size()-1;i++)
                {
                    F[0][i]=f(0,h*i);
                    F[u.size()-1][i]=f(1,h*i);
                    F[i][0]=f.diff_y(h*i,0);
                    F[i][u.size()-1]=f.diff_y(h*i,1);
                }
            }

            double la_error=0,la_resi=0;
            
            if(SC==rela_accuracy)
            {
                for(int i=1;i<=20&&(residual_norm(1)>st_parm);i++)
                {
                    if(CY==V_cycle)u=Vcycle(v1,v2,u,F);else u=FMG(v1,v2,F);
                    double now_error=error_norm(1);
                    double now_resi=residual_norm(1);
                    if(show_detail)printf(" %d & %.3e & %.3lf & %.3e & %.3lf\\\\ \\hline \n",i,now_error,la_error/now_error,now_resi,la_resi/now_resi);
                    la_error=now_error;la_resi=now_resi;
                }
                    
            }else
            {
                for(int i=1;i<=st_parm;i++)
                {
                    if(CY==V_cycle)u=Vcycle(v1,v2,u,F);else u=FMG(v1,v2,F);
                    double now_error=error_norm(1);
                    double now_resi=residual_norm(1);
                    if(show_detail)printf(" %d & %.3e & %.3lf & %.3e & %.3lf\\\\ \\hline \n",i,now_error,la_error/now_error,now_resi,la_resi/now_resi);
                    la_error=now_error;la_resi=now_resi;
                }
                    
            }
            
        }
        
        
        double residual_norm(double p)
        {
            double re=0;
            vector<vector<double> > F;F.clear();
            F.resize(u.size());
            for(int i=0;i<(int)u.size();i++)
                for(int j=0;j<(int)u.size();j++)
                    F[i].push_back(f.laplace(h*i,h*j));

            if(BC==Dirichlet)
            {
                for(int i=0;i<(int)u.size();i++)
                {
                    F[0][i]=f(0,h*i);
                    F[u.size()-1][i]=f(1,h*i);
                    F[i][0]=f(h*i,0);
                    F[i][u.size()-1]=f(h*i,1);
                }
            }
            else if(BC==Neumann)
            {
                for(int i=1;i<(int)u.size()-1;i++)
                {
                    F[0][i]=f.diff_x(0,h*i);
                    F[u.size()-1][i]=f.diff_x(1,h*i);
                    F[i][0]=f.diff_y(h*i,0);
                    F[i][u.size()-1]=f.diff_y(h*i,1);
                }
                F[0][u.size()/2]=f(0,0.5);
            }else 
            {
                for(int i=1;i<(int)u.size()-1;i++)
                {
                    F[0][i]=f(0,h*i);
                    F[u.size()-1][i]=f(1,h*i);
                    F[i][0]=f.diff_y(h*i,0);
                    F[i][u.size()-1]=f.diff_y(h*i,1);
                }
            }

            vector<vector<double> > res=F-laplace(u);
            for(int i=0;i<u.size();i++)
                for(int j=0;j<u.size();j++)
                {
                    if((i==0&&j==0||i==0&&j==u.size()-1||i==u.size()-1&&j==0||i==u.size()-1&&j==u.size()-1))continue;
                    re+=h*h*pow(fabs(res[i][j]),p);
                }
            return pow(re,1.0/p);
        }
        double error_norm(double p)
        {
            double re=0;
            
            for(int i=0;i<u.size();i++)
                for(int j=0;j<u.size();j++)
                {
                    if((i==0&&j==0||i==0&&j==u.size()-1||i==u.size()-1&&j==0||i==u.size()-1&&j==u.size()-1))continue;
                    re+=h*h*pow(fabs(u[i][j]-f(h*i,h*j)),p);
                }
                    
            return pow(re,1.0/p);
        }
};