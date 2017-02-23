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

Class
    Foam::solidBodyMotionFunctions::surfaceRegionMotion

Description
    SolidBodyMotionFunction governing the motion of the surface region. The motion
    function simply retrieves the septernion that is calculated in the zigZagMotionSolver
    of turningCircleMotionSolver from the main region (region0).

SourceFiles
    surfaceRegionMotion.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceRegionMotion_H
#define surfaceRegionMotion_H

#include "PiroSolidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"
#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class surfaceRegionMotion Declaration
\*---------------------------------------------------------------------------*/

class surfaceRegionMotion
:
    public PiroSolidBodyMotionFunction
{
    // Private data

    // Private Member Functions

        //- Disallow copy construct
        surfaceRegionMotion(const surfaceRegionMotion&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceRegionMotion&);


public:

    //- Runtime type information
    TypeName("surfaceRegionMotion");


    // Constructors

        //- Construct from components
        surfaceRegionMotion
        (
            const dictionary& SBMFCoeffs,
//            const Time& runTime
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~surfaceRegionMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
//        virtual septernion transformation(const septernion& TR);
        virtual septernion transformation();

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
