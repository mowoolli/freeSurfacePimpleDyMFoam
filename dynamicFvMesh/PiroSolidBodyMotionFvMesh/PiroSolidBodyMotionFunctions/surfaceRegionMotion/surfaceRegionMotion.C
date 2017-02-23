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
#include "IOReferencer.H"

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
    scalar t = time_.value();

    // Find the septernion associated with the main mesh (region0)    
    const fvMesh& pMesh = time_.db().parent().lookupObject<fvMesh>("region0");
    IOdictionary septDict = pMesh.lookupObject<IOdictionary>("septDict");

    // Initialize the septernion here in the surface mesh to be the same as region0
    septernion mySept = septDict.lookup("TRpassed");

    Info<< "solidBodyMotionFunctions::surfaceRegionMotion::transformation(): "
        << "Time = " << t << " transformation: " << mySept << endl;

    return mySept;
}


bool Foam::solidBodyMotionFunctions::surfaceRegionMotion::read
(
    const dictionary& SBMFCoeffs
)
{

    return true;
}


// ************************************************************************* //
