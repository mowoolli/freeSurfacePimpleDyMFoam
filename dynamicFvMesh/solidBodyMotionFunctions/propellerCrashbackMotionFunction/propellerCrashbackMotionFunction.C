/*---------------------------------------------------------------------------*\
| Portions 2015 Copyright (C) Marc O. Woolliscroft                            |
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Motion function for the description of a rotating appendage from a body-fixed
    reference frame. User inputs are origin, axis, frequency, and a linear ramp time.

\*---------------------------------------------------------------------------*/

#include "propellerCrashbackMotionFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include <string>

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(propellerCrashbackMotionFunction, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        propellerCrashbackMotionFunction,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::propellerCrashbackMotionFunction::propellerCrashbackMotionFunction
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    origin_(SBMFCoeffs_.lookup("origin")),
    axis_(SBMFCoeffs_.lookup("axis"))
//    omega_(DataEntry<scalar>::New("omega", SBMFCoeffs_))
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::propellerCrashbackMotionFunction::~propellerCrashbackMotionFunction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::propellerCrashbackMotionFunction::transformation() const
{
    scalar t = time_.value();
    
    scalar omegaMax = omega_;

    // Rotation origin (in radians)
    vector omega
    (
        omegaMax*axis_.x(),
        omegaMax*axis_.y(),
        omegaMax*axis_.z()
//        t*degToRad(radialVelocity_.x()),
//        t*degToRad(radialVelocity_.y()),
//        t*degToRad(radialVelocity_.z())
    );
/*
    scalar magOmega = mag(omega);
    quaternion R(omega/magOmega, magOmega);
    septernion TR(septernion(origin_)*R*septernion(-origin_));
*/
    vector eulerAngles(0, 0, 0);

    if (t < trFullAhead_ )
    {
        eulerAngles = omega*t;
//        eulerAngles = omega_->integrate(0, t)*axis_;
        Info << "eulerAngles = " << eulerAngles << endl;
    }
    else if ((t > trFullAhead_ || t == trFullAhead_) && (t < (trFullAhead_ + trSlowDown_)))
    {
        eulerAngles = 0.5*omega*trSlowDown_/pi*(cos(pi*trFullAhead_/trSlowDown_)*sin(pi*t/trSlowDown_) - sin(pi*trFullAhead_/trSlowDown_)*cos(pi*t/trSlowDown_)) + 0.5*omega*t + 0.5*omega*trFullAhead_;
        Info << "eulerAngles = " << eulerAngles << endl;
        Info << "Slowing down to 0 rpm..." << endl;
    }
    else if ((t == (trFullAhead_ + trSlowDown_) || t > (trFullAhead_ + trSlowDown_)) && (t < (trFullAhead_ + 2*trSlowDown_)))
    {
        eulerAngles = 0.5*omega*trSlowDown_/pi*(cos(pi*(trFullAhead_+trSlowDown_)/trSlowDown_)*sin(pi*t/trSlowDown_)-sin(pi*(trFullAhead_+trSlowDown_)/trSlowDown_)*cos(pi*t/trSlowDown_)) - 0.5*omega*t + omega*trFullAhead_ + 0.5*omega*(trFullAhead_+2*trSlowDown_);
        Info << "eulerAngles = " << eulerAngles << endl;
        Info << "Speeding up to full negative rpm..." << endl;
    }
    else
    {
        eulerAngles = -1*omega*t + 2*omega*trFullAhead_ + 2*omega*trSlowDown_;
        Info << "eulerAngles = " << eulerAngles << endl;
//        eulerAngles = omega_->integrate(0, t)*axis_;
    }

    quaternion R(eulerAngles/mag(eulerAngles), mag(eulerAngles));
    septernion TR(septernion(origin_)*R*septernion(-origin_));

    Info<< "solidBodyMotionFunctions::propellerCrashbackMotionFunction::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::propellerCrashbackMotionFunction::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("trFullAhead") >> trFullAhead_;
    SBMFCoeffs_.lookup("trSlowDown") >> trSlowDown_;
    SBMFCoeffs_.lookup("omega") >> omega_;

/*
    omega_.reset
    (
        DataEntry<scalar>::New("omega", SBMFCoeffs_).ptr()
    );
*/
    return true;
}

// ************************************************************************* //
