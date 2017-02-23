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

#include "propellerMotionFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include <string>

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(propellerMotionFunction, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        propellerMotionFunction,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::propellerMotionFunction::propellerMotionFunction
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    origin_(SBMFCoeffs_.lookup("origin")),
    axis_(SBMFCoeffs_.lookup("axis")),
    omega_(DataEntry<scalar>::New("omega", SBMFCoeffs_))
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::propellerMotionFunction::~propellerMotionFunction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::propellerMotionFunction::transformation() const
{
    scalar t = time_.value();

    vector eulerAngles(0, 0, 0);

    if (t < trp_ )
    {
        eulerAngles = omega_->integrate(0, t)*axis_*t/trp_*0.5;
    }
    else
    {
        eulerAngles = omega_->integrate(0, t)*axis_;
    }

    quaternion R(eulerAngles/mag(eulerAngles), mag(eulerAngles));
    septernion TR(septernion(origin_)*R*septernion(-origin_));

    Info<< "solidBodyMotionFunctions::propellerMotionFunction::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::propellerMotionFunction::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("trp") >> trp_;

    omega_.reset
    (
        DataEntry<scalar>::New("omega", SBMFCoeffs_).ptr()
    );

    return true;
}

// ************************************************************************* //
