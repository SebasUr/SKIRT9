#ifndef DIFFUSEREEMISSIONWEIGHTS_HPP
#define DIFFUSEREEMISSIONWEIGHTS_HPP
#include "FatalError.hpp"
#include <cmath>
#include <limits>

// Photon-number probabilities, conditional only at an actual scattering event.
struct DiffuseReemissionWeights
{
    double weights[4] = {};
    double total = 0.;

    static DiffuseReemissionWeights make(double a, const double p[5], double q, double f)
    {
        auto valid=[](double v) { return std::isfinite(v) && v>=0. && v<=1.; };
        if(!valid(a) || !valid(q) || !valid(f))
            throw FATALERROR("Invalid diffuse absorption/OTS/two-photon probability");
        for(int i=0;i<5;++i) if(!valid(p[i]))
            throw FATALERROR("Diffuse recombination probability outside [0,1]");
        const double tolerance=64.*std::numeric_limits<double>::epsilon();
        if(p[1]+p[2]+p[3]+p[4]>1.+tolerance)
            throw FATALERROR("Helium branching probabilities sum above one");
        DiffuseReemissionWeights result;
        result.weights[0]=a*p[0]+(1.-a)*p[4]*q*p[0];
        result.weights[1]=(1.-a)*p[1];
        result.weights[2]=(1.-a)*p[2];
        result.weights[3]=(1.-a)*(p[3]+p[4]*(1.-q))*f;
        for(double v:result.weights) result.total+=v;
        if(!std::isfinite(result.total) || result.total<0. || result.total>1.+tolerance)
            throw FATALERROR("Invalid diffuse ionizing survival probability");
        // Only machine-roundoff overshoot is normalized; invalid fits are rejected.
        if(result.total>1.) {
            for(double& v:result.weights) v/=result.total;
            result.total=1.;
        }
        return result;
    }

    int select(double unit) const
    {
        if(!std::isfinite(unit) || unit<0. || unit>=1.)
            throw FATALERROR("Invalid conditional channel uniform deviate");
        if(total<=0.)return -1;
        const double target=unit*total;
        double cumulative=0.;int last=-1;
        for(int i=0;i<4;++i) if(weights[i]>0.) {
            last=i;cumulative+=weights[i];
            if(target<cumulative)return i;
        }
        return last; // roundoff at the final cumulative sum, never a zero-weight channel
    }
};
#endif
