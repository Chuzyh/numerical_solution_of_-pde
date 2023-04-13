#include "BVP.h"

class Multigrid_Method
{
private:
    enum bound_conditon { Dirichlet, Neumann,mixed} BC;
    enum restriction_operators {full_weighting,injection} RO;
    enum interpolation_operators {linear,quadratic} IO;
    enum cycles {V_cycle,FMG} CY;
    enum stopping_criteria{max_iteration,rela_accuracy} SC;
    int max_iter;double esplion;
    
public:

};