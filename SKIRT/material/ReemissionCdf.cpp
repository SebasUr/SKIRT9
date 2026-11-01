/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ReemissionCdf.hpp"
#include "FatalError.hpp"
#include <atomic>

////////////////////////////////////////////////////////////////////

void ReemissionCdf::initialize(const StoredTable<3>& table)
{
    if (table.axisUsesLogInterpolation<0>() || !table.axisUsesLogInterpolation<1>()
        || table.axisUsesLogInterpolation<2>() || table.quantityUsesLogInterpolation())
        throw FATALERROR("Unsupported reemission table interpolation metadata");

    Array types, temperatures, wavelengths;
    table.axisArray<0>(types);
    table.axisArray<1>(temperatures);
    table.axisArray<2>(wavelengths);
    if (!types.size() || temperatures.size() < 2 || wavelengths.size() < 2)
        throw FATALERROR("Empty reemission CDF table");

    _temperatures.resize(temperatures.size());
    _wavelengths.resize(wavelengths.size());
    for (size_t t = 0; t < temperatures.size(); ++t)
    {
        double value = temperatures[t];
        if (!std::isfinite(value) || value <= 0. || (t && value <= temperatures[t - 1]))
            throw FATALERROR("Reemission temperatures must be positive and ascending");
        _temperatures[t] = value;
    }
    bool ascending = wavelengths[1] > wavelengths[0];
    for (size_t i = 0; i < wavelengths.size(); ++i)
    {
        double value = wavelengths[i];
        if (!std::isfinite(value) || value <= 0.
            || (i && (ascending ? value <= wavelengths[i - 1] : value >= wavelengths[i - 1])))
            throw FATALERROR("Reemission wavelengths must be positive and strictly monotone");
        _wavelengths[i] = value;
    }

    _table = &table;
    _prefixIndices.assign(types.size(), {});
    _lazySafe.assign(types.size(), true);
    for (size_t k = 0; k < types.size(); ++k)
    {
        if (types[k] != static_cast<double>(k)) throw FATALERROR("Unsupported reemission type axis");
        for (size_t t = 0; t < temperatures.size(); ++t)
        {
            vector<size_t> indices(wavelengths.size());
            size_t best = 0;
            double maximum = -1.;
            for (size_t i = 0; i < wavelengths.size(); ++i)
            {
                double p = table.valueAtIndices(k, t, i);
                if (!std::isfinite(p) || p < 0.) throw FATALERROR("Invalid raw reemission probability");

                // taking the latest maximum makes the map the identity for monotone rows, including plateaus
                if (p >= maximum)
                {
                    maximum = p;
                    best = i;
                }
                indices[i] = best;
            }
            if (!(maximum > 0.)) throw FATALERROR("Zero reemission cumulative probability");
            if (!t)
                _prefixIndices[k] = indices;
            else if (indices != _prefixIndices[k])
                _lazySafe[k] = false;
        }
    }

    static std::atomic<size_t> nextInstance{1};
    _instance = nextInstance.fetch_add(1, std::memory_order_relaxed);
}

////////////////////////////////////////////////////////////////////

void ReemissionCdf::materialize(int type, double temperature, vector<double>& probabilities) const
{
    checkType(type);
    Bracket b = bracket(temperature);
    probabilities.resize(_wavelengths.size());
    for (size_t i = 0; i < probabilities.size(); ++i) probabilities[i] = valueAt(type, b, i);
    for (size_t i = 1; i < probabilities.size(); ++i)
        if (probabilities[i] < probabilities[i - 1]) probabilities[i] = probabilities[i - 1];
    double normalization = 1. / probabilities.back();
    for (double& p : probabilities) p *= normalization;
}

////////////////////////////////////////////////////////////////////

double ReemissionCdf::sample(int type, double temperature, double u) const
{
    checkType(type);
    if (!_lazySafe[type]) return sampleMaterialized(type, temperature, u);

    // a position that holds the running maximum at both bracketing temperatures also holds it for any interpolation
    // with nonnegative weights, so that the running maximum can be taken before interpolating
    Bracket b = bracket(temperature);
    const auto& indices = _prefixIndices[type];
    double normalization = 1. / valueAt(type, b, indices.back());
    auto p = [&](size_t i) { return valueAt(type, b, indices[i]) * normalization; };

    size_t left = 0;
    size_t right = _wavelengths.size();
    while (left < right)
    {
        size_t middle = left + (right - left) / 2;
        if (p(middle) < u)
            left = middle + 1;
        else
            right = middle;
    }
    return invert(left, u, left && left < _wavelengths.size() ? p(left - 1) : 0.,
                  left < _wavelengths.size() ? p(left) : 0.);
}

////////////////////////////////////////////////////////////////////

ReemissionCdf::Bracket ReemissionCdf::bracket(double temperature) const
{
    if (!std::isfinite(temperature) || temperature <= 0.)
        throw FATALERROR("Reemission temperature must be finite and positive");
    auto it = std::lower_bound(_temperatures.begin(), _temperatures.end(), temperature);
    if (it == _temperatures.begin()) return {0, 0, 0.};
    if (it == _temperatures.end()) return {_temperatures.size() - 1, _temperatures.size() - 1, 0.};
    size_t hi = it - _temperatures.begin();
    double fraction = (std::log(temperature) - std::log(_temperatures[hi - 1]))
                      / (std::log(_temperatures[hi]) - std::log(_temperatures[hi - 1]));
    return {hi - 1, hi, fraction};
}

////////////////////////////////////////////////////////////////////

double ReemissionCdf::valueAt(size_t type, const Bracket& bracket, size_t index) const
{
    double left = _table->valueAtIndices(type, bracket.lo, index);
    if (bracket.lo == bracket.hi) return left;
    double right = _table->valueAtIndices(type, bracket.hi, index);
    return (1. - bracket.fraction) * left + bracket.fraction * right;
}

////////////////////////////////////////////////////////////////////

void ReemissionCdf::checkType(int type) const
{
    if (!_table || type < 0 || static_cast<size_t>(type) >= _lazySafe.size())
        throw FATALERROR("Invalid reemission CDF type or uninitialized sampler");
}

////////////////////////////////////////////////////////////////////

double ReemissionCdf::invert(size_t index, double u, double p1, double p2) const
{
    if (index == _wavelengths.size()) return _wavelengths.back();
    if (index == 0) return _wavelengths[0];
    double lambda1 = _wavelengths[index - 1];
    double lambda2 = _wavelengths[index];
    if (p2 != p1)
    {
        double fraction = (u - p1) / (p2 - p1);
        return lambda1 + fraction * (lambda2 - lambda1);
    }
    return lambda1;
}

////////////////////////////////////////////////////////////////////

double ReemissionCdf::sampleMaterialized(int type, double temperature, double u) const
{
    struct Last
    {
        size_t instance{0};
        int type{-1};
        double temperature{0.};
        vector<double> probabilities;
    };
    thread_local Last last;
    if (last.instance != _instance || last.type != type || last.temperature != temperature)
    {
        materialize(type, temperature, last.probabilities);
        last.instance = _instance;
        last.type = type;
        last.temperature = temperature;
    }
    const auto& p = last.probabilities;
    size_t i = std::lower_bound(p.begin(), p.end(), u) - p.begin();
    return invert(i, u, i && i < p.size() ? p[i - 1] : 0., i < p.size() ? p[i] : 0.);
}

////////////////////////////////////////////////////////////////////
