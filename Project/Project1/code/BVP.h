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
    virtual double result(double x, double y){return 0.0;}  
    virtual double error_norm(double p){return 0;}  
};

class FD_regular : public FD_Methods
{
private:
    Function & f;
    double h;
    int N,cond;

};