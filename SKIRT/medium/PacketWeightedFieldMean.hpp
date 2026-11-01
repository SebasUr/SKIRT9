/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PACKETWEIGHTEDFIELDMEAN_HPP
#define PACKETWEIGHTEDFIELDMEAN_HPP

#include "Array.hpp"
#include "Basics.hpp"

////////////////////////////////////////////////////////////////////

/** The PacketWeightedFieldMean class keeps a running mean over iterations of a radiation field that
    has already been synchronized and normalized, weighted by the number of packets launched in
    each iteration. After the burn-in iterations, each new field is replaced by the running mean,
    so that the mean is used by whatever consumes the field. A final pass is folded into the mean
    in the same way, and ends the lifetime of the object. An instance is owned by a single
    MediumSystem object, and its functions are called serially. */
class PacketWeightedFieldMean
{
public:
    /** This function prepares the object for a field with the given number of values, excluding
        the given number of initial iterations from the mean. */
    void configure(size_t size, int burnIn);

    /** This function folds the given field into the mean and replaces it by the mean, unless the
        iteration is part of the burn-in. Iterations must be numbered consecutively starting with
        one; the final pass has iteration number zero. The function returns true if the field was
        folded into the mean. */
    bool blend(Array& field, size_t packets, int iteration, bool finalPass);

    /** This function returns the accumulated packet weight. */
    double weight() const { return _weight; }

    /** This function returns the number of fields folded into the mean. */
    int samples() const { return _samples; }

    /** This function returns the number of bytes allocated for the mean. */
    size_t payloadBytes() const { return _mean.capacity() * sizeof(double); }

private:
    vector<double> _mean;
    double _weight{0.};
    int _burnIn{0};
    int _nextIteration{1};
    int _samples{0};
    bool _finalized{false};
};

////////////////////////////////////////////////////////////////////

#endif
