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

#include "surfaceRegionMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "rigidBody6DoFMotion.H"
#include "PiroSolidBodyMotionFunction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(surfaceRegionMotion, 0);
    addToRunTimeSelectionTable
    (
        PiroSolidBodyMotionFunction,
        surfaceRegionMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::surfaceRegionMotion::surfaceRegionMotion
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

Foam::solidBodyMotionFunctions::surfaceRegionMotion::~surfaceRegionMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::surfaceRegionMotion::transformation()
{
    const dictionary& SBMFCoeffs = PiroSolidBodyMotionFunction::SBMFCoeffs_;
    const fvMesh& mesh = PiroSolidBodyMotionFunction::mesh_;

    scalar t = time_.value();
 
    // bill add
    scalar scalarTest;
    IOdictionary septDict(IOobject("septernion.dat", time_.constant(), mesh_, IOobject::MUST_READ, IOobject::AUTO_WRITE));
    septDict.readIfPresent("scalarTest", scalarTest);
    Info << "bill scalarTest = " << scalarTest << endl;
    // bill done

    // Translation of centre of gravity with constant velocity
    const vector displacement = dummyVector_*t;

//    rigidBody6DoFMotion TRorig(SBMFCoeffs, mesh);
    rigidBody6DoFMotion* TRcopy = new rigidBody6DoFMotion(SBMFCoeffs, mesh);

//    Info << "TRorig = " << TRorig.polyMeshTransformation() << endl;
//    Info << "TRcopy = " << TRcopy->Foam::solidBodyMotionFunctions::rigidBody6DoFMotion::polyMeshTransformation() << endl;

    quaternion R(0, 0, 0);
    septernion TR(septernion(displacement)*R);

//    septernion TRsurface = TRorig.polyMeshTransformation();
//    septernion TRsurface = TRcopy->Foam::solidBodyMotionFunctions::rigidBody6DoFMotion::polyMeshTransformation();
    scalar TRsurface = TRcopy->Foam::solidBodyMotionFunctions::rigidBody6DoFMotion::polyMeshTransformation();
//    scalar TRsurface = TRorig.Foam::solidBodyMotionFunctions::rigidBody6DoFMotion::polyMeshTransformation();

//    Info << "attempt to get scalarTest = " << Foam::solidBodyMotionFunctions::rigidBody6DoFMotion::scalarTest << endl;

    Info<< "solidBodyMotionFunctions::surfaceRegionMotion::transformation(): "
        << "Time = " << t << " transformation: " << TRsurface << endl;

    Info << "TR that is returned in surface = " << TR << endl;

//    delete TRcopy;

    return TR;
//    return TRsurface;
}


bool Foam::solidBodyMotionFunctions::surfaceRegionMotion::read
(
    const dictionary& SBMFCoeffs
)
{
//    PiroSolidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("dummyVector") >> dummyVector_;

    return true;
}


// ************************************************************************* //
