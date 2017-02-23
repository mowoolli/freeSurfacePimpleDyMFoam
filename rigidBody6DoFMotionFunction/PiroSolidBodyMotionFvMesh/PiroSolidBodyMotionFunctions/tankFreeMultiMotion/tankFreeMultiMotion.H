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
    Foam::solidBodyMotionFunctions::tankFreeMultiMotion

Description
    Combination of SolidBodyMotionFvMesh 6DoF motion functions.

SourceFiles
    tankFreeMultiMotion.C

\*---------------------------------------------------------------------------*/

#ifndef tankFreeMultiMotion_H
#define tankFreeMultiMotion_H

#include "PiroSolidBodyMotionFunction.H"
//#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class tankFreeMultiMotion Declaration
\*---------------------------------------------------------------------------*/

class tankFreeMultiMotion
:
    public PiroSolidBodyMotionFunction
{
    // Private data

        //- Motions to combine
        PtrList<PiroSolidBodyMotionFunction> SBMFs_;


    // Private Member Functions

        //- Disallow copy construct
        tankFreeMultiMotion(const tankFreeMultiMotion&);

        //- Disallow default bitwise assignment
        void operator=(const tankFreeMultiMotion&);


public:

    //- Runtime type information
    TypeName("tankFreeMultiMotion");


    // Constructors

        //- Construct from components
        tankFreeMultiMotion
        (
            const dictionary& SBMFCoeffs,
//            const Time& runTime
            const fvMesh& mesh
        );
/*
        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new tankFreeMultiMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }
*/

    //- Destructor
    virtual ~tankFreeMultiMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
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
