/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ITERATIONOPTIONS_HPP
#define ITERATIONOPTIONS_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The IterationOptions class simply offers a number of options for configuring iterations during
    primary and/or secondary emission. These options are relevant only when the simulation has a
    dynamic medium state and/or a dynamic secondary emission.

    If the \em averagePrimaryField flag is enabled, the primary radiation field used to update the
    medium state is replaced, after \em primaryFieldAverageBurnIn iterations, by the running mean
    of the synchronized fields weighted by the number of packets. The final primary emission
    segment is folded into the same mean. The option is available only for a single static medium
    with a DiffuseIonizedGasMix, without dust or secondary iterations. */
class IterationOptions : public SimulationItem
{
    ITEM_CONCRETE(IterationOptions, SimulationItem,
                  "a set of options for configuring iterations during primary and/or secondary emission")

        PROPERTY_BOOL(averagePrimaryField, "average the primary radiation field over iterations after a burn-in")
        ATTRIBUTE_DEFAULT_VALUE(averagePrimaryField, "false")
        ATTRIBUTE_RELEVANT_IF(averagePrimaryField, "IteratePrimary")
        ATTRIBUTE_DISPLAYED_IF(averagePrimaryField, "Level3")

        PROPERTY_INT(primaryFieldAverageBurnIn, "the number of initial primary iterations excluded from averaging")
        ATTRIBUTE_MIN_VALUE(primaryFieldAverageBurnIn, "0")
        ATTRIBUTE_MAX_VALUE(primaryFieldAverageBurnIn, "1000")
        ATTRIBUTE_DEFAULT_VALUE(primaryFieldAverageBurnIn, "5")
        ATTRIBUTE_RELEVANT_IF(primaryFieldAverageBurnIn, "IteratePrimary&averagePrimaryField")
        ATTRIBUTE_DISPLAYED_IF(primaryFieldAverageBurnIn, "Level3")

        PROPERTY_INT(minPrimaryIterations, "the minimum number of iterations during primary emission")
        ATTRIBUTE_MIN_VALUE(minPrimaryIterations, "1")
        ATTRIBUTE_MAX_VALUE(minPrimaryIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(minPrimaryIterations, "1")
        ATTRIBUTE_RELEVANT_IF(minPrimaryIterations, "IteratePrimary")

        PROPERTY_INT(maxPrimaryIterations, "the maximum number of iterations during primary emission")
        ATTRIBUTE_MIN_VALUE(maxPrimaryIterations, "1")
        ATTRIBUTE_MAX_VALUE(maxPrimaryIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(maxPrimaryIterations, "10")
        ATTRIBUTE_RELEVANT_IF(maxPrimaryIterations, "IteratePrimary")

        PROPERTY_INT(minSecondaryIterations, "the minimum number of iterations during secondary emission")
        ATTRIBUTE_MIN_VALUE(minSecondaryIterations, "1")
        ATTRIBUTE_MAX_VALUE(minSecondaryIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(minSecondaryIterations, "1")
        ATTRIBUTE_RELEVANT_IF(minSecondaryIterations, "IterateSecondary")

        PROPERTY_INT(maxSecondaryIterations, "the maximum number of iterations during secondary emission")
        ATTRIBUTE_MIN_VALUE(maxSecondaryIterations, "0")
        ATTRIBUTE_MAX_VALUE(maxSecondaryIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(maxSecondaryIterations, "10")
        ATTRIBUTE_RELEVANT_IF(maxSecondaryIterations, "IterateSecondary")

        PROPERTY_BOOL(includePrimaryEmission, "include primary emission in the secondary emission iterations")
        ATTRIBUTE_DEFAULT_VALUE(includePrimaryEmission, "false")
        ATTRIBUTE_RELEVANT_IF(includePrimaryEmission, "IterateSecondary")
        ATTRIBUTE_DISPLAYED_IF(includePrimaryEmission, "Level3")

        PROPERTY_DOUBLE(primaryIterationPacketsMultiplier,
                        "the multiplier on the number of photon packets launched for each primary emission iteration")
        ATTRIBUTE_MIN_VALUE(primaryIterationPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(primaryIterationPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(primaryIterationPacketsMultiplier, "1")
        ATTRIBUTE_RELEVANT_IF(primaryIterationPacketsMultiplier, "IteratePrimary|includePrimaryEmission")
        ATTRIBUTE_DISPLAYED_IF(primaryIterationPacketsMultiplier, "Level3")

        PROPERTY_DOUBLE(primaryIterationInitialPacketsFraction,
                        "the fraction of primary iteration photon packets to use in the first iteration")
        ATTRIBUTE_MIN_VALUE(primaryIterationInitialPacketsFraction, "]0")
        ATTRIBUTE_MAX_VALUE(primaryIterationInitialPacketsFraction, "1]")
        ATTRIBUTE_DEFAULT_VALUE(primaryIterationInitialPacketsFraction, "1")
        ATTRIBUTE_RELEVANT_IF(primaryIterationInitialPacketsFraction, "IteratePrimary")
        ATTRIBUTE_DISPLAYED_IF(primaryIterationInitialPacketsFraction, "Level3")

        PROPERTY_DOUBLE(primaryIterationPacketsRamp,
                        "the geometric ramp-up factor per primary iteration for the photon packet count")
        ATTRIBUTE_MIN_VALUE(primaryIterationPacketsRamp, "]1")
        ATTRIBUTE_MAX_VALUE(primaryIterationPacketsRamp, "100]")
        ATTRIBUTE_DEFAULT_VALUE(primaryIterationPacketsRamp, "2")
        ATTRIBUTE_RELEVANT_IF(primaryIterationPacketsRamp, "IteratePrimary&primaryIterationInitialPacketsFraction")
        ATTRIBUTE_DISPLAYED_IF(primaryIterationPacketsRamp, "Level3")

        PROPERTY_DOUBLE(secondaryIterationPacketsMultiplier,
                        "the multiplier on the number of photon packets launched for each secondary emission iteration")
        ATTRIBUTE_MIN_VALUE(secondaryIterationPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(secondaryIterationPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(secondaryIterationPacketsMultiplier, "1")
        ATTRIBUTE_RELEVANT_IF(secondaryIterationPacketsMultiplier, "IterateSecondary")
        ATTRIBUTE_DISPLAYED_IF(secondaryIterationPacketsMultiplier, "Level3")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
