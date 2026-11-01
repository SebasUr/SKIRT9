/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef REEMISSIONCDF_HPP
#define REEMISSIONCDF_HPP

#include "StoredTable.hpp"

////////////////////////////////////////////////////////////////////

/** The ReemissionCdf class samples wavelengths from the temperature-dependent cumulative
    distributions in a reemission table of type StoredTable<3>, with axes spectrum type,
    temperature and wavelength. The distributions are interpolated logarithmically in temperature
    and by wavelength index, so the result is independent of the order of the wavelength axis.

    Rebuilding a distribution for each new temperature takes a time proportional to the number of
    wavelengths. If, for a given spectrum type, the positions of the running maxima are the same at
    all tabulated temperatures, the distribution can instead be inverted by binary search on the
    two bracketing temperatures, without building it. These positions are determined once during
    initialization. For the other spectrum types the distribution is built, and the most recent
    one is retained per execution thread. */
class ReemissionCdf
{
public:
    /** This function verifies the contents of the given table and stores a pointer to it. The
        table must remain available for the lifetime of this object. */
    void initialize(const StoredTable<3>& table);

    /** This function returns the wavelength axis of the table, in table order. */
    const vector<double>& wavelengths() const { return _wavelengths; }

    /** This function stores in \em probabilities the normalized cumulative distribution for the
        given spectrum type and temperature. */
    void materialize(int type, double temperature, vector<double>& probabilities) const;

    /** This function returns the wavelength corresponding to the uniform deviate \em u for the
        given spectrum type and temperature. */
    double sample(int type, double temperature, double u) const;

private:
    /** This structure holds the temperature bracket and interpolation fraction. */
    struct Bracket
    {
        size_t lo, hi;
        double fraction;
    };

    /** This function returns the bracket for the given temperature, clamped to the table range. */
    Bracket bracket(double temperature) const;

    /** This function returns the interpolated table value at the given wavelength index. */
    double valueAt(size_t type, const Bracket& bracket, size_t index) const;

    /** This function verifies the spectrum type. */
    void checkType(int type) const;

    /** This function returns the wavelength at the given position in the distribution, using the
        probabilities \em p1 and \em p2 on either side of the position. */
    double invert(size_t index, double u, double p1, double p2) const;

    /** This function samples after building the distribution. */
    double sampleMaterialized(int type, double temperature, double u) const;

    const StoredTable<3>* _table{nullptr};
    vector<double> _temperatures;
    vector<double> _wavelengths;
    vector<vector<size_t>> _prefixIndices;  // for each type, the index of the running maximum at each wavelength
    vector<bool> _lazySafe;                 // for each type, true if these indices hold at all temperatures
    size_t _instance{0};                    // distinguishes objects in the per-thread cache
};

////////////////////////////////////////////////////////////////////

#endif
