#ifndef REEMISSIONCDF_HPP
#define REEMISSIONCDF_HPP
#include "StoredTable.hpp"
#include "FatalError.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>

// Index-ordered CDF; independent of whether physical wavelengths ascend or descend.
class ReemissionCdf
{
    const StoredTable<3>* _table = nullptr;
    vector<double> _temperatures, _wavelengths;
    vector<vector<size_t>> _prefixIndices;
    vector<bool> _lazySafe;
    size_t _instance = 0;

    struct Bracket { size_t lo, hi; double fraction; };
    Bracket bracket(double temperature) const
    {
        if (!std::isfinite(temperature) || temperature <= 0.)
            throw FATALERROR("Reemission temperature must be finite and positive");
        auto it = std::lower_bound(_temperatures.begin(), _temperatures.end(), temperature);
        if (it == _temperatures.begin()) return {0, 0, 0.};
        if (it == _temperatures.end())
            return {_temperatures.size()-1, _temperatures.size()-1, 0.};
        size_t hi = it - _temperatures.begin();
        double f = (std::log(temperature) - std::log(_temperatures[hi-1]))
                   / (std::log(_temperatures[hi]) - std::log(_temperatures[hi-1]));
        return {hi-1, hi, f};
    }
    double at(size_t type, const Bracket& b, size_t i) const
    {
        double left = _table->valueAtIndices(type, b.lo, i);
        if (b.lo == b.hi) return left;
        double right = _table->valueAtIndices(type, b.hi, i);
        return (1.-b.fraction)*left + b.fraction*right;
    }
    void checkType(int type) const
    {
        if (!_table || type < 0 || static_cast<size_t>(type) >= _lazySafe.size())
            throw FATALERROR("Invalid reemission CDF type or uninitialized sampler");
    }
    double invert(size_t index, double u, double p1, double p2) const
    {
        if (index == _wavelengths.size()) return _wavelengths.back();
        if (index == 0) return _wavelengths[0];
        double lambda1 = _wavelengths[index-1], lambda2 = _wavelengths[index];
        if (p2 != p1)
        {
            double fraction = (u-p1)/(p2-p1);
            return lambda1 + fraction*(lambda2-lambda1);
        }
        return lambda1;
    }

public:
    void initialize(const StoredTable<3>& table)
    {
        // Scope the contract explicitly. No silent reinterpretation of another pack.
        if (table.axisUsesLogInterpolation<0>() || !table.axisUsesLogInterpolation<1>()
            || table.axisUsesLogInterpolation<2>() || table.quantityUsesLogInterpolation())
            throw FATALERROR("Unsupported reemission table interpolation metadata");
        Array types, temperatures, wavelengths;
        table.axisArray<0>(types); table.axisArray<1>(temperatures); table.axisArray<2>(wavelengths);
        if (!types.size() || temperatures.size()<2 || wavelengths.size()<2)
            throw FATALERROR("Empty reemission CDF table");
        _temperatures.resize(temperatures.size()); _wavelengths.resize(wavelengths.size());
        for (size_t t=0; t<temperatures.size(); ++t)
        {
            double v=temperatures[t];
            if (!std::isfinite(v) || v<=0. || (t && v<=temperatures[t-1]))
                throw FATALERROR("Reemission temperatures must be positive and ascending");
            _temperatures[t]=v;
        }
        bool ascending=wavelengths[1]>wavelengths[0];
        for (size_t i=0; i<wavelengths.size(); ++i)
        {
            double v=wavelengths[i];
            if (!std::isfinite(v) || v<=0. || (i && (ascending ? v<=wavelengths[i-1] : v>=wavelengths[i-1])))
                throw FATALERROR("Reemission wavelengths must be positive and strictly monotone");
            _wavelengths[i]=v;
        }
        _table=&table;
        _prefixIndices.assign(types.size(), {}); _lazySafe.assign(types.size(), true);
        for (size_t k=0; k<types.size(); ++k)
        {
            if (types[k]!=static_cast<double>(k)) throw FATALERROR("Unsupported reemission type axis");
            for (size_t t=0; t<temperatures.size(); ++t)
            {
                vector<size_t> indices(wavelengths.size());
                size_t best=0; double maximum=-1.;
                for (size_t i=0; i<wavelengths.size(); ++i)
                {
                    double p=table.valueAtIndices(k,t,i);
                    if (!std::isfinite(p) || p<0.) throw FATALERROR("Invalid raw reemission probability");
                    // Latest maximizer makes the map the identity for monotone rows,
                    // including plateaus. A shared map is sufficient, not necessary.
                    if (p>=maximum) { maximum=p; best=i; }
                    indices[i]=best;
                }
                if (!(maximum>0.)) throw FATALERROR("Zero reemission cumulative probability");
                if (!t) _prefixIndices[k]=indices;
                else if (indices!=_prefixIndices[k]) _lazySafe[k]=false;
            }
        }
        static std::atomic<size_t> nextInstance{1};
        _instance=nextInstance.fetch_add(1, std::memory_order_relaxed);
    }

    const vector<double>& wavelengths() const { return _wavelengths; }
    bool lazySafe(int type) const { checkType(type); return _lazySafe[type]; }
    void materialize(int type, double temperature, vector<double>& probabilities) const
    {
        checkType(type); Bracket b=bracket(temperature);
        probabilities.resize(_wavelengths.size());
        for (size_t i=0; i<probabilities.size(); ++i) probabilities[i]=at(type,b,i);
        for (size_t i=1; i<probabilities.size(); ++i)
            if (probabilities[i]<probabilities[i-1]) probabilities[i]=probabilities[i-1];
        double normalization=1./probabilities.back();
        for (double& p:probabilities) p*=normalization;
    }
    double sampleMaterialized(int type, double temperature, double u) const
    {
        checkType(type);
        struct Last { size_t instance=0; int type=-1; double temperature=0.; vector<double> p; };
        thread_local Last last;
        if (last.instance!=_instance || last.type!=type || last.temperature!=temperature)
        {
            materialize(type,temperature,last.p);
            last.instance=_instance; last.type=type; last.temperature=temperature;
        }
        auto it=std::lower_bound(last.p.begin(),last.p.end(),u);
        size_t i=it-last.p.begin();
        return invert(i,u,i && i<last.p.size() ? last.p[i-1] : 0.,i<last.p.size() ? last.p[i] : 0.);
    }
    double sampleLazy(int type, double temperature, double u) const
    {
        checkType(type);
        if (!_lazySafe[type]) return sampleMaterialized(type,temperature,u);
        Bracket b=bracket(temperature);
        const auto& indices=_prefixIndices[type];
        // A common prefix maximizer at every temperature endpoint also maximizes
        // every nonnegative linear interpolation of those endpoints. If the maps
        // differ, the materialized fallback above preserves max/interpolate order.
        double normalization=1./at(type,b,indices.back());
        auto p=[&](size_t i) { return at(type,b,indices[i])*normalization; };
        size_t left=0, right=_wavelengths.size();
        while (left<right)
        {
            size_t mid=left+(right-left)/2;
            if (p(mid)<u) left=mid+1;
            else right=mid;
        }
        return invert(left,u,left && left<_wavelengths.size() ? p(left-1) : 0.,
                      left<_wavelengths.size() ? p(left) : 0.);
    }
};
#endif
