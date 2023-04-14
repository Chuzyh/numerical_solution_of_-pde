class F1: public Function2d
{
    public:
    double operator()(const double &x,const double &y)const{return exp(y+sin(x));}
    double diff_x(const double &x,const double &y)const{return cos(x)*exp(y+sin(x));}
    double diff_y(const double &x,const double &y)const{return exp(y+sin(x));}
    double diff2_x(const double &x,const double &y)const{return -sin(x)*exp(y+sin(x))+cos(x)*cos(x)*exp(y+sin(x));}
    double diff2_y(const double &x,const double &y)const{return exp(y+sin(x));}
}Fun;

class F2: public Function2d
{
    public:
    double operator()(const double &x,const double &y)const{return sin(x*3+y*3);}
    double diff_x(const double &x,const double &y)const{return 3*cos(x*3+y*3);}
    double diff_y(const double &x,const double &y)const{return 3*cos(x*3+y*3);}
    double diff2_x(const double &x,const double &y)const{return -9*sin(x*3+y*3);}
    double diff2_y(const double &x,const double &y)const{return -9*sin(x*3+y*3);}
}Fun2;
class F3: public Function2d
{
    public:
    double operator()(const double &x,const double &y)const{return exp(x*x+y*y*y);}
    double diff_x(const double &x,const double &y)const{return 2*x*exp(x*x+y*y*y);}
    double diff_y(const double &x,const double &y)const{return 3*y*y*exp(x*x+y*y*y);}
    double diff2_x(const double &x,const double &y)const{return 2*exp(x*x+y*y*y)+4*x*x*exp(x*x+y*y*y);}
    double diff2_y(const double &x,const double &y)const{return 6*y*exp(x*x+y*y*y)+9*y*y*y*y*exp(x*x+y*y*y);}
}Fun3;

class F4: public Function
{
    public:
    double operator()(const double &x)const{return exp(x);}
    double diff(const double &x)const{return exp(x);}
    double laplace(const double &x)const{return -exp(x);}
}Fun4;