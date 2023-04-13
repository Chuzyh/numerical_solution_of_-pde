#ifndef BVP
#define BVP
#include<bits/stdc++.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
using namespace Eigen;
using namespace std;
const double EPS=1e-7;
void ERROR(string wrongmessage)
{
    cerr<<wrongmessage<<endl;
    exit(-1);
}
class Function
{
    private:
        int dim;
    public:
        //two dimension
        virtual double operator () (const double &x,const double &y) const =0;
        virtual double diff_x(const double &x,const double &y) 
        const{
            if(dim!=2)ERROR("wrong dim");
            return ((*this)(x+EPS,y)-(*this)(x-EPS,y))/(EPS*2);
        }
        virtual double diff2_x(const double &x,const double &y) 
        const{
            if(dim!=2)ERROR("wrong dim");
            return ((*this).diff_x(x+EPS,y)-(*this).diff_x(x-EPS,y))/(EPS*2);
        }
        virtual double diff_y(const double &x,const double &y) 
        const{
            if(dim!=2)ERROR("wrong dim");
            return ((*this)(x,y+EPS)-(*this)(x,y-EPS))/(EPS*2);
        }
        virtual double diff2_y(const double &x,const double &y) 
        const{
            if(dim!=2)ERROR("wrong dim");
            return ((*this).diff_y(x,y+EPS)-(*this).diff_y(x,y-EPS))/(EPS*2);
        }
        virtual double laplace(const double &x,const double &y) 
        const{
            if(dim!=2)ERROR("wrong dim");
            return -(*this).diff2_x(x,y)-(*this).diff2_y(x,y);
        }

        //one dimension
        virtual double operator () (const double &x) const =0;
        virtual double diff(const double &x) 
        const{
            if(dim!=1)ERROR("wrong dim");
            return ((*this)(x+EPS)-(*this)(x-EPS))/(EPS*2);
        }
        virtual double laplace(const double &x) 
        const{
            if(dim!=1)ERROR("wrong dim");
            return ((*this).diff(x+EPS)-(*this).diff(x-EPS))/(EPS*2);
        }
        
};
class point
{
    private:
        double x,y;
    public:
        point (double x,double y):x(x),y(y){}
        point operator +(point a){return (point){a.x+(*this).x,a.y+(*this).y};} 
        point operator -(point a){return (point){(*this).x-a.x,(*this).y-a.y};} 
        point operator *(double a){return (point){(*this).x*a,(*this).y*a};} 
        point operator /(double a){return (point){(*this).x/a,(*this).y/a};} 
        double len(){return sqrt((*this).x*(*this).x+(*this).y*(*this).y);}
        void out(){printf("%lf %lf ",(*this).x,(*this).y);}
        double X(){return (*this).x;}
        double Y(){return (*this).y;}
};
class Circle
{
    private:
        point cen;
        double radius;
    public:
        Circle (point cen,double radius):cen(cen),radius(radius){}
        point get_center(){return (*this).cen;}
        double get_radius(){return (*this).radius;}
        int in_square(){return (*this).cen.X()+(*this).radius<1&&(*this).cen.X()-(*this).radius>0&&(*this).cen.Y()-(*this).radius>0&&(*this).cen.Y()-(*this).radius>0;}
        pair<double,double> get_inter(int vertical,double pos)
        {
            if(vertical)
            {
                double Y1=cen.Y()-sqrt(radius*radius-(cen.X()-pos)*(cen.X()-pos));
                double Y2=cen.Y()+sqrt(radius*radius-(cen.X()-pos)*(cen.X()-pos));
                return make_pair(Y1,Y2);
            }else
            {
                double X1=cen.X()-sqrt(radius*radius-(cen.Y()-pos)*(cen.Y()-pos));
                double X2=cen.X()+sqrt(radius*radius-(cen.Y()-pos)*(cen.Y()-pos));
                return make_pair(X1,X2);
            }
        }
};
class FD_Methods
{
public:
    virtual void solve() = 0;  
    virtual double operator ()(double x, double y){return 0.0;}  
    virtual double error_norm(double p){return 0;}  
};

class FD_regular : public FD_Methods
{
private:
    Function & f;
    double h;
    int N,cond;
    VectorXd u;
public:
    FD_regular (Function & f,double h,int cond):f(f),h(h),cond(cond)
    {
        if(!(cond==1||cond==2||cond==3)){cerr<<"invalid condition";exit(1);}
        N=(int)(1.0/h+1);
    }
    void solve()
    {
        if (cond==1) // Dirichlet
        {
            SparseMatrix<double> A(N*N,N*N);
            vector<Triplet<double> > a;
            MatrixXd b(N*N,1);
            b=MatrixXd::Zero(N*N,1);
            // 从左到右，从下到上依次从0开始标号至N*N-1
            for(int x=0;x<N;x++)
            for(int y=0;y<N;y++)
            {
                int num=x+y*N;
                double X=x*h;
                double Y=y*h;

                if(x==0||x==N-1||y==0||y==N-1)
                {
                    a.push_back(Triplet<double>(num,num,1.0));
                    b(num,0)=f(X,Y);
                }
                else
                {
                    double h2=h*h;
                    a.push_back(Triplet<double>(num,num,4/h2));
                    a.push_back(Triplet<double>(num,num+1,-1/h2));
                    a.push_back(Triplet<double>(num,num-1,-1/h2));
                    a.push_back(Triplet<double>(num,num+N,-1/h2));
                    a.push_back(Triplet<double>(num,num-N,-1/h2));
                    b(num,0)=f.laplace(X,Y);
                }
		    }
            A.setFromTriplets(a.begin(),a.end());
            
            A.makeCompressed();
            SparseLU<SparseMatrix<double> > Solver;
            Solver.compute(A);
            u=Solver.solve(b);
        }
        if (cond==2) // Neumann
        {
            SparseMatrix<double> A(N*N,N*N);
            vector<Triplet<double> > a;
            MatrixXd b(N*N,1);
            b=MatrixXd::Zero(N*N,1);

            
            // 从左到右，从下到上依次从0开始标号至N*N-1
            for(int x=0;x<N;x++)
            for(int y=0;y<N;y++)
            {
                int num=x+y*N;
                double X=x*h;
                double Y=y*h;

                if(x==0&&y==0||x==0&&y==N-1||x==N-1&&y==0||x==N-1&&y==N-1)//四个角上的特殊处理一下
                {
                    a.push_back(Triplet<double>(num,num,1.0));
                    b(num,0)=f(X,Y);
                }else if(x==0||x==N-1||y==0||y==N-1)
                {
                    double h2=h*h;
                    if(x==0)
                    {
                        if(y==1)//加一个额外的限制条件让系数矩阵可逆
                        {
                            a.push_back(Triplet<double>(num,num,1.0));
                            b(num,0)=f(X,Y);
                            continue;
                        }
                        a.push_back(Triplet<double>(num,num,-1.5/h));
                        a.push_back(Triplet<double>(num,num+1,2.0/h));
                        a.push_back(Triplet<double>(num,num+2,-0.5/h));
                        b(num,0)=f.diff_x(X,Y);
                    } 
                    if(x==N-1)
                    {
                        a.push_back(Triplet<double>(num,num,1.5/h));
                        a.push_back(Triplet<double>(num,num-1,-2.0/h));
                        a.push_back(Triplet<double>(num,num-2,0.5/h));
                        b(num,0)=f.diff_x(X,Y);
                    } 
                    if(y==0)
                    {
                        a.push_back(Triplet<double>(num,num,-1.5/h));
                        a.push_back(Triplet<double>(num,num+N,2.0/h));
                        a.push_back(Triplet<double>(num,num+N*2,-0.5/h));
                        b(num,0)=f.diff_y(X,Y);
                    }
                    if(y==N-1)
                    {
                        a.push_back(Triplet<double>(num,num,1.5/h));
                        a.push_back(Triplet<double>(num,num-N,-2.0/h));
                        a.push_back(Triplet<double>(num,num-N*2,0.5/h));
                        b(num,0)=f.diff_y(X,Y);
                    }
                }
                else
                {
                    double h2=h*h;
                    a.push_back(Triplet<double>(num,num,4.0/h2));
                    a.push_back(Triplet<double>(num,num+1,-1/h2));
                    a.push_back(Triplet<double>(num,num-1,-1/h2));
                    a.push_back(Triplet<double>(num,num+N,-1/h2));
                    a.push_back(Triplet<double>(num,num-N,-1/h2));
                    b(num,0)=f.laplace(X,Y);
                    //cout<<num<<endl;
                }
		    }
            A.setFromTriplets(a.begin(),a.end());
            A.makeCompressed();
            // cout<<A<<endl;
            SparseLU<SparseMatrix<double> > Solver;
            Solver.compute(A);
            u=Solver.solve(b);
        }
        if (cond==3) // mix
        {
            SparseMatrix<double> A(N*N,N*N);
            vector<Triplet<double> > a;
            MatrixXd b(N*N,1);
            b=MatrixXd::Zero(N*N,1);

            
            // 从左到右，从下到上依次从0开始标号至N*N-1
            for(int x=0;x<N;x++)
            for(int y=0;y<N;y++)
            {
                int num=x+y*N;
                double X=x*h;
                double Y=y*h;

                if(x==0&&y==0||x==0&&y==N-1||x==N-1&&y==0||x==N-1&&y==N-1)//四个角上的特殊处理一下
                {
                    a.push_back(Triplet<double>(num,num,1.0));
                    b(num,0)=f(X,Y);
                }else if(x==0||x==N-1||y==0||y==N-1)
                {
                    double h2=h*h;
                    
                    if(y==0)
                    {
                        //找一条边上选Neumann边值条件
                        a.push_back(Triplet<double>(num,num,-1.5/h));
                        a.push_back(Triplet<double>(num,num+N,2.0/h));
                        a.push_back(Triplet<double>(num,num+N*2,-0.5/h));
                        b(num,0)=f.diff_y(X,Y);continue;
                    }
                    a.push_back(Triplet<double>(num,num,1.0));
                    b(num,0)=f(X,Y);
                }
                else
                {
                    double h2=h*h;
                    a.push_back(Triplet<double>(num,num,4.0/h2));
                    a.push_back(Triplet<double>(num,num+1,-1/h2));
                    a.push_back(Triplet<double>(num,num-1,-1/h2));
                    a.push_back(Triplet<double>(num,num+N,-1/h2));
                    a.push_back(Triplet<double>(num,num-N,-1/h2));
                    b(num,0)=f.laplace(X,Y);
                }
		    }
            A.setFromTriplets(a.begin(),a.end());
            A.makeCompressed();
            SparseLU<SparseMatrix<double> > Solver;
            Solver.compute(A);
            u=Solver.solve(b);
        }
    }
    double operator ()(double x,double y)
    {
        if(fabs(x/h-(int)(x/h))>EPS||fabs(y/h-(int)(y/h))>EPS)
        {
            cerr<<"not a grid point"<<endl;
            assert(0);
        }
        return u[(int)(x/h)+((int)(y/h))*N];
    }
    double error_norm(double p)//-1 表示无穷范数
    {
        double re=0,h2=h*h;
        if(p<0)
        {
            for(int x=0;x<N;x++)
                for(int y=0;y<N;y++)
                    re=max(re,fabs(f(h*x,h*y)-(*this)(h*x,h*y)));
            return re;
        }else
        {
            for(int x=0;x<N;x++)
                for(int y=0;y<N;y++)
                    re+=h2*pow(fabs(f(h*x,h*y)-(*this)(h*x,h*y)),p);
            return pow(re,1.0/p);
        }
    }
};
class FD_irregular : public FD_Methods
{
private:
    Function & f;
    double h;
    int N,cond,realN;
    VectorXd u;
    Circle C;
    map<pair<int,int>,pair<int,int> > label;//first int 表示编号， second int 表示类型（1 区域内点，2 圆中的ghost point，3 正方形边界，4 圆边界）

public:
    FD_irregular (Function & f,int cond,double h,Circle C):f(f),cond(cond),h(h),C(C)
    {
        if(!(cond==1||cond==2||cond==3)){cerr<<"invalid condition";exit(1);}
        N=(int)(1.0/h+1);realN=0;
        if(!C.in_square()){cerr<<"circle is not in the square";exit(1);}
        label.clear();
        get_label();
    }
    void get_label()
    {
        for(int x=0;x<N;x++)
            for(int y=0;y<N;y++)
            {
                double X=x*h;
                double Y=y*h;
                point now=(point){X,Y};
                if(x==0||x==N-1||y==0||y==N-1){label[make_pair(x,y)]=make_pair(0,3);}
                else if(fabs((C.get_center()-now).len()-C.get_radius())<EPS){label[make_pair(x,y)]=make_pair(0,4);}
                else if((C.get_center()-now).len()>C.get_radius()){label[make_pair(x,y)]=make_pair(0,1);}
            }
        for(int x=0;x<N;x++)
            for(int y=0;y<N;y++)
            {
                if(label.count(make_pair(x,y))==0)
                {
                    if(label.count(make_pair(x+1,y)))
                    {
                        label[make_pair(x,y)]=make_pair(0,2);
                        continue;
                    }
                    if(label.count(make_pair(x-1,y)))
                    {
                        label[make_pair(x,y)]=make_pair(0,2);
                        continue;
                    }if(label.count(make_pair(x,y+1)))
                    {
                        label[make_pair(x,y)]=make_pair(0,2);
                        continue;
                    }if(label.count(make_pair(x,y-1)))
                    {
                        label[make_pair(x,y)]=make_pair(0,2);
                        continue;
                    }
                    
                }
            }
        for(map<pair<int,int>,pair<int,int> >::iterator it=label.begin();it!=label.end();it++)
        {
            (*it).second.first=realN++;
        }
    }
    void solve()
    {
        if (cond==1) // Dirichlet
        {
            SparseMatrix<double> A(realN,realN);
            vector<Triplet<double> > a;
            MatrixXd b(realN,1);
            b=MatrixXd::Zero(realN,1);
            for(map<pair<int,int>,pair<int,int> >::iterator it=label.begin();it!=label.end();it++)
            {
                int num=(*it).second.first,type=(*it).second.second,x=(*it).first.first,y=(*it).first.second;
                double X=x*h,Y=y*h,h2=h*h;
                if(type==1)
                {
                    a.push_back(Triplet<double>(num,num,4/h2));
                    a.push_back(Triplet<double>(num,label[make_pair(x-1,y)].first,-1/h2));
                    a.push_back(Triplet<double>(num,label[make_pair(x+1,y)].first,-1/h2));
                    a.push_back(Triplet<double>(num,label[make_pair(x,y-1)].first,-1/h2));
                    a.push_back(Triplet<double>(num,label[make_pair(x,y+1)].first,-1/h2));
                    b(num,0)=f.laplace(X,Y);
                }else if(type==2)
                {
                    if(label.count(make_pair(x-1,y)))
                    {
                            pair<double,double> inter=C.get_inter(0,Y);
                            double X_inter=(inter.first>=X-h&&inter.first<=X)?inter.first:inter.second;
                            double f_inter=f(X_inter,Y);
                            a.push_back(Triplet<double>(num,num,1.0/(X-X_inter)));
                            a.push_back(Triplet<double>(num,label[make_pair(x-1,y)].first,1.0/(X_inter-X+h)));
                            b(num,0)=(1.0/(X-X_inter) + 1.0/(X_inter-X+h))*f_inter;
                            continue;
                    }
                    if(label.count(make_pair(x+1,y)))
                    {
                            pair<double,double> inter=C.get_inter(0,Y);
                            double X_inter=(inter.first>=X&&inter.first<=X+h)?inter.first:inter.second;
                            double f_inter=f(X_inter,Y);
                            a.push_back(Triplet<double>(num,num,1.0/(X_inter-X)));
                            a.push_back(Triplet<double>(num,label[make_pair(x+1,y)].first,1.0/(X+h-X_inter)));
                            b(num,0)=(1.0/(X_inter-X) + 1.0/(X+h-X_inter))*f_inter;
                            continue;
                    }
                    if(label.count(make_pair(x,y-1)))
                    {
                            pair<double,double> inter=C.get_inter(1,X);
                            double Y_inter=(inter.first>=Y-h&&inter.first<=Y)?inter.first:inter.second;
                            double f_inter=f(X,Y_inter);
                            a.push_back(Triplet<double>(num,num,1.0/(Y-Y_inter)));
                            a.push_back(Triplet<double>(num,label[make_pair(x,y-1)].first,1.0/(Y_inter-Y+h)));
                            b(num,0)=(1.0/(Y-Y_inter) + 1.0/(Y_inter-Y+h))*f_inter;
                            continue;
                    }
                    if(label.count(make_pair(x,y+1)))
                    {
                            pair<double,double> inter=C.get_inter(1,X);
                            double Y_inter=(inter.first>=Y&&inter.first<=Y+h)?inter.first:inter.second;
                            double f_inter=f(X,Y_inter);
                            a.push_back(Triplet<double>(num,num,1.0/(Y_inter-Y)));
                            a.push_back(Triplet<double>(num,label[make_pair(x,y+1)].first,1.0/(Y+h-Y_inter)));
                            b(num,0)=(1.0/(Y_inter-Y) + 1.0/(Y+h-Y_inter))*f_inter;
                            continue;
                    }
                }else if(type==3||type==4)
                {
                    a.push_back(Triplet<double>(num,num,1.0));
                    b(num,0)=f(X,Y);
                }
            }

            A.setFromTriplets(a.begin(),a.end());
            
            A.makeCompressed();
            SparseLU<SparseMatrix<double> > Solver;
            Solver.compute(A);
            u=Solver.solve(b);
        }
        if (cond==2) // Neumann
        {
            cerr<<"ddl is coming, write it later";
            exit(-1);
        }
        if (cond==3) // mix 圆上的Neumann太难写了，所以我们假设圆上是Dirichlet 条件，方形边界上是Neumann条件   
        {
            SparseMatrix<double> A(realN,realN);
            vector<Triplet<double> > a;
            MatrixXd b(realN,1);
            b=MatrixXd::Zero(realN,1);
            for(map<pair<int,int>,pair<int,int> >::iterator it=label.begin();it!=label.end();it++)
            {
                int num=(*it).second.first,type=(*it).second.second,x=(*it).first.first,y=(*it).first.second;
                double X=x*h,Y=y*h,h2=h*h;
                if(type==1)
                {
                    a.push_back(Triplet<double>(num,num,4/h2));
                    a.push_back(Triplet<double>(num,label[make_pair(x-1,y)].first,-1/h2));
                    a.push_back(Triplet<double>(num,label[make_pair(x+1,y)].first,-1/h2));
                    a.push_back(Triplet<double>(num,label[make_pair(x,y-1)].first,-1/h2));
                    a.push_back(Triplet<double>(num,label[make_pair(x,y+1)].first,-1/h2));
                    b(num,0)=f.laplace(X,Y);
                }else if(type==2)
                {
                    if(label.count(make_pair(x-1,y)))
                    {
                            pair<double,double> inter=C.get_inter(0,Y);
                            double X_inter=(inter.first>=X-h&&inter.first<=X)?inter.first:inter.second;
                            double f_inter=f(X_inter,Y);
                            a.push_back(Triplet<double>(num,num,1.0/(X-X_inter)));
                            a.push_back(Triplet<double>(num,label[make_pair(x-1,y)].first,1.0/(X_inter-X+h)));
                            b(num,0)=(1.0/(X-X_inter) + 1.0/(X_inter-X+h))*f_inter;
                            continue;
                    }
                    if(label.count(make_pair(x+1,y)))
                    {
                            pair<double,double> inter=C.get_inter(0,Y);
                            double X_inter=(inter.first>=X&&inter.first<=X+h)?inter.first:inter.second;
                            double f_inter=f(X_inter,Y);
                            a.push_back(Triplet<double>(num,num,1.0/(X_inter-X)));
                            a.push_back(Triplet<double>(num,label[make_pair(x+1,y)].first,1.0/(X+h-X_inter)));
                            b(num,0)=(1.0/(X_inter-X) + 1.0/(X+h-X_inter))*f_inter;
                            continue;
                    }
                    if(label.count(make_pair(x,y-1)))
                    {
                            pair<double,double> inter=C.get_inter(1,X);
                            double Y_inter=(inter.first>=Y-h&&inter.first<=Y)?inter.first:inter.second;
                            double f_inter=f(X,Y_inter);
                            a.push_back(Triplet<double>(num,num,1.0/(Y-Y_inter)));
                            a.push_back(Triplet<double>(num,label[make_pair(x,y-1)].first,1.0/(Y_inter-Y+h)));
                            b(num,0)=(1.0/(Y-Y_inter) + 1.0/(Y_inter-Y+h))*f_inter;
                            continue;
                    }
                    if(label.count(make_pair(x,y+1)))
                    {
                            pair<double,double> inter=C.get_inter(1,X);
                            double Y_inter=(inter.first>=Y&&inter.first<=Y+h)?inter.first:inter.second;
                            double f_inter=f(X,Y_inter);
                            a.push_back(Triplet<double>(num,num,1.0/(Y_inter-Y)));
                            a.push_back(Triplet<double>(num,label[make_pair(x,y+1)].first,1.0/(Y+h-Y_inter)));
                            b(num,0)=(1.0/(Y_inter-Y) + 1.0/(Y+h-Y_inter))*f_inter;
                            continue;
                    }
                }else if(type==3)
                {
                    if(x==0&&y==0||x==0&&y==N-1||x==N-1&&y==0||x==N-1&&y==N-1)//四个角上的特殊处理一下
                    {
                        a.push_back(Triplet<double>(num,num,1.0));
                        b(num,0)=f(X,Y);
                        continue;
                    }
                    double h2=h*h;
                    if(x==0)
                    {
                        if(y==1)//加一个额外的限制条件让系数矩阵可逆
                        {
                            a.push_back(Triplet<double>(num,num,1.0));
                            b(num,0)=f(X,Y);
                            continue;
                        }
                        a.push_back(Triplet<double>(num,num,-1.5/h));
                        a.push_back(Triplet<double>(num,label[make_pair(x+1,y)].first,2.0/h));
                        a.push_back(Triplet<double>(num,label[make_pair(x+2,y)].first,-0.5/h));
                        b(num,0)=f.diff_x(X,Y);
                    } 
                    if(x==N-1)
                    {
                        a.push_back(Triplet<double>(num,num,1.5/h));
                        a.push_back(Triplet<double>(num,label[make_pair(x-1,y)].first,-2.0/h));
                        a.push_back(Triplet<double>(num,label[make_pair(x-2,y)].first,0.5/h));
                        b(num,0)=f.diff_x(X,Y);
                    } 
                    if(y==0)
                    {
                        a.push_back(Triplet<double>(num,num,-1.5/h));
                        a.push_back(Triplet<double>(num,label[make_pair(x,y+1)].first,2.0/h));
                        a.push_back(Triplet<double>(num,label[make_pair(x,y+2)].first,-0.5/h));
                        b(num,0)=f.diff_y(X,Y);
                    }
                    if(y==N-1)
                    {
                        a.push_back(Triplet<double>(num,num,1.5/h));
                        a.push_back(Triplet<double>(num,label[make_pair(x,y-1)].first,-2.0/h));
                        a.push_back(Triplet<double>(num,label[make_pair(x,y-2)].first,0.5/h));
                        b(num,0)=f.diff_y(X,Y);
                    }
                }else
                {
                    a.push_back(Triplet<double>(num,num,1.0));
                    b(num,0)=f(X,Y);
                }
            }

            A.setFromTriplets(a.begin(),a.end());
            
            A.makeCompressed();
            SparseLU<SparseMatrix<double> > Solver;
            Solver.compute(A);
            u=Solver.solve(b);
        }
        
    }
    double operator ()(int x,int y)
    {
        return u[label[make_pair(x,y)].first];
    }
    double error_norm(double p)//-1 表示无穷范数
    {
        double re=0,h2=h*h;
        if(p<0)
        {
            for(map<pair<int,int>,pair<int,int> >::iterator it=label.begin();it!=label.end();it++)
            {
                int x=(*it).first.first;
                int y=(*it).first.second;
                re=max(re,fabs(f(h*x,h*y)-(*this)(x,y)));
            }
                
            return re;
        }else
        {
            for(map<pair<int,int>,pair<int,int> >::iterator it=label.begin();it!=label.end();it++)
            {
                int x=(*it).first.first;
                int y=(*it).first.second;
                re+=h2*pow(fabs(f(h*x,h*y)-(*this)(x,y)),p);
            }
            return pow(re,1.0/p);
        }
    }
};
#endif