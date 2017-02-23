/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "PiroSolidBodyMotionFunction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
//defineTypeNameAndDebug(Foam::PiroSolidBodyMotionFunction, 0);
defineTypeNameAndDebug(PiroSolidBodyMotionFunction, 0);
//defineRunTimeSelectionTable(Foam::PiroSolidBodyMotionFunction, dictionary);
defineRunTimeSelectionTable(PiroSolidBodyMotionFunction, dictionary);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PiroSolidBodyMotionFunction::PiroSolidBodyMotionFunction
(
    const dictionary& SBMFCoeffs,
    //const Time& runTime
    const fvMesh& mesh
)
:
    SBMFCoeffs_
    (
        SBMFCoeffs.subDict
        (
            word(SBMFCoeffs.lookup("PiroSolidBodyMotionFunction")) + "Coeffs"
        )
    ),
    mesh_(mesh),
    time_(mesh.time())
    //time_(runTime)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PiroSolidBodyMotionFunction::~PiroSolidBodyMotionFunction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::PiroSolidBodyMotionFunction::read(const dictionary& SBMFCoeffs)
{
    SBMFCoeffs_ = SBMFCoeffs.subDict(type() + "Coeffs");

    return true;
}


// ************************************************************************* //
