/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "PacketWeightedFieldMean.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // the largest packet count for which every integer can be represented exactly in double precision
    constexpr double maxExactWeight = 9007199254740992.;
}

////////////////////////////////////////////////////////////////////

void PacketWeightedFieldMean::configure(size_t size, int burnIn)
{
    if (!size || burnIn < 0) throw FATALERROR("Invalid field mean configuration");
    vector<double>(size, 0.).swap(_mean);
    _burnIn = burnIn;
    _weight = 0.;
    _nextIteration = 1;
    _samples = 0;
    _finalized = false;
}

////////////////////////////////////////////////////////////////////

bool PacketWeightedFieldMean::blend(Array& field, size_t packets, int iteration, bool finalPass)
{
    if (_mean.empty() || field.size() != _mean.size() || _finalized)
        throw FATALERROR("Invalid field mean lifecycle or shape");
    if (!packets || static_cast<double>(packets) > maxExactWeight) throw FATALERROR("Invalid field mean packet weight");
    if ((!finalPass && iteration != _nextIteration) || (finalPass && iteration != 0))
        throw FATALERROR("Field mean iteration sequence changed");
    for (size_t i = 0; i != field.size(); ++i)
        if (!(field[i] >= 0.) || !std::isfinite(field[i]))
            throw FATALERROR("Nonfinite or negative primary radiation field");

    if (!finalPass)
    {
        ++_nextIteration;
        if (iteration <= _burnIn) return false;
    }

    double incoming = static_cast<double>(packets);
    if (_weight > maxExactWeight - incoming) throw FATALERROR("Field mean packet weight exceeds exact integer range");
    double total = _weight + incoming;
    if (_samples == 0)
    {
        for (size_t i = 0; i != field.size(); ++i) _mean[i] = field[i];
    }
    else
    {
        double oldFraction = _weight / total;
        double newFraction = incoming / total;
        for (size_t i = 0; i != field.size(); ++i)
        {
            _mean[i] = oldFraction * _mean[i] + newFraction * field[i];
            if (!std::isfinite(_mean[i])) throw FATALERROR("Nonfinite primary field mean");
            field[i] = _mean[i];
        }
    }
    _weight = total;
    ++_samples;
    if (finalPass) _finalized = true;
    return true;
}

////////////////////////////////////////////////////////////////////
