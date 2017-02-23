/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "linearForcesMotion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(linearForcesMotion, 0);
    addToRunTimeSelectionTable
    (
        PiroSolidBodyMotionFunction,
        linearForcesMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::linearForcesMotion::linearForcesMotion
(
    const dictionary& SBMFCoeffs,
//    const Time& runTime
    const fvMesh& mesh
)
:
//    solidBodyMotionFunction(SBMFCoeffs, runTime)
    PiroSolidBodyMotionFunction(SBMFCoeffs, mesh)
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::linearForcesMotion::~linearForcesMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::linearForcesMotion::transformation()
{
    dictionary forcesDict = SBMFCoeffs_.subDict("forces");
//    forcesDict.remove("CofR");
//    forcesDict.add("CofR",(CofG_+xOld_));

    forces f
    (
        "forces",
        mesh_,
//        SBMFCoeffs_.subDict("forces")
        forcesDict
    );

    scalar t = time_.value();

    // Translation of centre of gravity with constant velocity
    const vector displacement = velocity_*t;

    quaternion R(0, 0, 0);
    septernion TR(septernion(displacement)*R);

    Info<< "solidBodyMotionFunctions::linearForcesMotion::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::linearForcesMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    PiroSolidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("velocity") >> velocity_;

    return true;
}

// ************************************************************************* //
