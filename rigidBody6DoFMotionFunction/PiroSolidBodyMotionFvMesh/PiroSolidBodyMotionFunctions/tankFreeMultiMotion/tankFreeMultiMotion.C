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

#include "tankFreeMultiMotion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(tankFreeMultiMotion, 0);
    addToRunTimeSelectionTable
    (
        PiroSolidBodyMotionFunction,
        tankFreeMultiMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::tankFreeMultiMotion::tankFreeMultiMotion
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

Foam::solidBodyMotionFunctions::tankFreeMultiMotion::~tankFreeMultiMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::tankFreeMultiMotion::transformation()
{
    scalar t = time_.value();

    septernion TR = SBMFs_[0].transformation();

    for (label i = 1; i < SBMFs_.size(); i++)
    {
        TR *= SBMFs_[i].transformation();
    }

    Info<< "solidBodyMotionFunctions::tankFreeMultiMotion::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::tankFreeMultiMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    PiroSolidBodyMotionFunction::read(SBMFCoeffs);

    label i = 0;
    SBMFs_.setSize(SBMFCoeffs_.size());

    forAllConstIter(IDLList<entry>, SBMFCoeffs_, iter)
    {
        if (iter().isDict())
        {
            SBMFs_.set
            (
                i,
//                solidBodyMotionFunction::New(iter().dict(), time_)
                PiroSolidBodyMotionFunction::New(iter().dict(), mesh_)
            );

            Info<< "Constructed SBMF " << i << " : "
                << iter().keyword() << " of type "
                << SBMFs_[i].type() << endl;

            i++;
        }
    }
    SBMFs_.setSize(i);

    return true;
}

// ************************************************************************* //
