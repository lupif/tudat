#ifndef ORBITALSTATEPARAMETERS_H
#define ORBITALSTATEPARAMETERS_H



namespace tudat
{
namespace orbital_elements
{

//! Cartesian elements.
enum CartesianElements
{
    xCartesianPosition = 0,
    yCartesianPosition = 1,
    zCartesianPosition = 2,
    xCartesianVelocity = 3,
    yCartesianVelocity = 4,
    zCartesianVelocity = 5
};

//! Keplerian elements.
enum KeplerianElements
{
    semiMajorAxis = 0,
    eccentricity = 1,
    inclination = 2,
    argumentOfPeriapsis = 3,
    longitudeOfAscendingNode = 4,
    trueAnomaly = 5,
    semiLatusRectum = 5
};

//! Modified equinoctial elements.
enum ModifiedEquinoctialElements
{
    semiParameter = 0,
    fElement = 1,
    gElement = 2,
    hElement = 3,
    kElement = 4,
    trueLongitude = 5
};

//! Spherical orbital state elements
enum SphericalOrbitalStateElements
{
    radius = 0,
    latitude = 1,
    longitude = 2,
    speed = 3,
    flightPath = 4,
    headingAngle = 5
};

//! Unified State Model elements
enum UnifiedStateModelElements
{
    CHodograph = 0,
    Rf1Hodograph = 1,
    Rf2Hodograph = 2,
    epsilon1Quaternion = 3,
    epsilon2Quaternion = 4,
    epsilon3Quaternion = 5,
    etaQuaternion = 6
};

//! Cartesian acceleration components
enum CartesianAccelerationElements
{
    xCartesianAcceleration,
    yCartesianAcceleration,
    zCartesianAcceleration
};

//! Acceleration components in CSN frame for orbital elements.
enum CSNAccelerationElements
{
    cAcceleration,
    sAcceleration,
    nAcceleration
};

}

}

#endif // ORBITALSTATEPARAMETERS_H
