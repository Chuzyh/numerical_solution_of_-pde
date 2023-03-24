#include<bits/stdc++.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
using namespace Eigen;
using namespace std;
const double EPS=1e-7;
class Function
{
    public:
        virtual double operator () (const double &x,const double &y) const =0;
        virtual double diff_x(const double &x,const double &y) 
        const{
            return ((*this)(x+EPS,y)-(*this)(x-EPS,y))/(EPS*2);
        }
        virtual double diff2_x(const double &x,const double &y) 
        const{
            return ((*this).diff_x(x+EPS,y)-(*this).diff_x(x-EPS,y))/(EPS*2);
        }
        virtual double diff_y(const double &x,const double &y) 
        const{
            return ((*this)(x,y+EPS)-(*this)(x,y-EPS))/(EPS*2);
        }
        virtual double diff2_y(const double &x,const double &y) 
        const{
            return ((*this).diff_y(x,y+EPS)-(*this).diff_y(x,y-EPS))/(EPS*2);
        }
        virtual double laplace(const double &x,const double &y) 
        const{
            return (*this).diff2_x(x,y)+(*this).diff2_y(x,y);
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
        if(!(cond==1||cond==2||cond==3)){cerr<<"invalid condition";assert(0);}
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
};