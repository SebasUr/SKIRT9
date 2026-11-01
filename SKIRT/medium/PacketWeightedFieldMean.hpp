#ifndef PACKETWEIGHTEDFIELDMEAN_HPP
#define PACKETWEIGHTEDFIELDMEAN_HPP
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

// Running mean of already luminosity-normalized field estimates, not packet sums.
// Owned by one MediumSystem, called serially after its existing MPI reduction.
class PacketWeightedFieldMean
{
    std::vector<double> _mean;
    double _weight{0.};
    int _burnIn{0}, _nextIteration{1}, _samples{0};
    bool _finalized{false};
public:
    void configure(std::size_t size, int burnIn)
    {
        if (!size || burnIn < 0) throw std::invalid_argument("Invalid field-mean configuration");
        std::vector<double>(size, 0.).swap(_mean);
        _burnIn = burnIn; _weight = 0.; _nextIteration = 1; _samples = 0; _finalized = false;
    }
    template<class Values> bool blend(Values& field, std::size_t packets, int iteration, bool finalPass)
    {
        if (_mean.empty() || field.size() != _mean.size() || _finalized)
            throw std::logic_error("Invalid field-mean lifecycle or shape");
        if (!packets || packets > 9007199254740992ULL)
            throw std::invalid_argument("Invalid field-mean packet weight");
        if ((!finalPass && iteration != _nextIteration) || (finalPass && iteration != 0))
            throw std::logic_error("Field-mean iteration sequence changed");
        for (std::size_t i = 0; i != field.size(); ++i)
            if (!(field[i] >= 0.) || !std::isfinite(field[i]))
                throw std::invalid_argument("Nonfinite or negative primary radiation field");
        if (!finalPass)
        {
            ++_nextIteration;
            if (iteration <= _burnIn) return false;
        }
        double incoming = static_cast<double>(packets);
        if (_weight > 9007199254740992. - incoming)
            throw std::overflow_error("Field-mean cumulative packet weight exceeds exact integer range");
        double total = _weight + incoming;
        if (_samples == 0)
        {
            for (std::size_t i = 0; i != field.size(); ++i) _mean[i] = field[i];
        }
        else
        {
            double oldFraction = _weight / total;
            double newFraction = incoming / total;
            for (std::size_t i = 0; i != field.size(); ++i)
            {
                _mean[i] = oldFraction * _mean[i] + newFraction * field[i];
                if (!std::isfinite(_mean[i])) throw std::overflow_error("Nonfinite primary field mean");
                field[i] = _mean[i];
            }
        }
        _weight = total; ++_samples;
        if (finalPass) _finalized = true;
        return true;
    }
    double weight() const { return _weight; }
    int samples() const { return _samples; }
    std::size_t payloadBytes() const { return _mean.capacity() * sizeof(double); }
};
#endif
