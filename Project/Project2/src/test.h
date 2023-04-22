class F1: public Function2d
{
    public:
    double operator()(const double &x,const double &y)const{return exp(y+sin(x));}
    double diff_x(const double &x,const double &y)const{return cos(x)*exp(y+sin(x));}
    double diff_y(const double &x,const double &y)const{return exp(y+sin(x));}
    double diff2_x(const double &x,const double &y)const{return -sin(x)*exp(y+sin(x))+cos(x)*cos(x)*exp(y+sin(x));}
    double diff2_y(const double &x,const double &y)const{return exp(y+sin(x));}
}Fun;

class F10: public Function2d
{
    public:
    double operator()(const double &x,const double &y)const{return x+y;}
    double diff_x(const double &x,const double &y)const{return 1;}
    double diff_y(const double &x,const double &y)const{return 1;}
    double diff2_x(const double &x,const double &y)const{return 0;}
    double diff2_y(const double &x,const double &y)const{return 0;}
}Fun0;

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
    double operator()(const double &x)const{return x;}
    double diff(const double &x)const{return 1;}
    double laplace(const double &x)const{return  0;}
}Fun4;
class F5: public Function
{
    public:
    double operator()(const double &x)const{return exp(sin(x))-1;}
    double diff(const double &x)const{return exp(sin(x))*cos(x);}
    double laplace(const double &x)const{return  exp(sin(x))*sin(x)-exp(sin(x))*cos(x)*cos(x);}
}Fun5;

class F6: public Function
{
    public:
    double operator()(const double &x)const{return sin(x);}
    double diff(const double &x)const{return cos(x);}
    double laplace(const double &x)const{return  sin(x);}
}Fun6;

class F7: public Function
{
    public:
    double operator()(const double &x)const{return exp(x*x);}
    double diff(const double &x)const{return 2*x*exp(x*x);}
    double laplace(const double &x)const{return  -2*exp(x*x)-4*x*x*exp(x*x);}
}Fun7;
