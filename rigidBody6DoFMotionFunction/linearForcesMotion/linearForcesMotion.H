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
    Foam::solidBodyMotionFunctions::linearForcesMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function. Constant velocity displacement.

SourceFiles
    linearForcesMotion.C

\*---------------------------------------------------------------------------*/

#ifndef linearForcesMotion_H
#define linearForcesMotion_H

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
                          Class linearForcesMotion Declaration
\*---------------------------------------------------------------------------*/

class linearForcesMotion
:
    public PiroSolidBodyMotionFunction
{
    // Private data

        //- Linear velocity
        vector velocity_;

//        objectRegistry& obr_;

    // Private Member Functions

        //- Disallow copy construct
        linearForcesMotion(const linearForcesMotion&);

        //- Disallow default bitwise assignment
        void operator=(const linearForcesMotion&);


public:

    //- Runtime type information
    TypeName("linearForcesMotion");


    // Constructors

        //- Construct from components
        linearForcesMotion
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
                new linearForcesMotion
                (
                    SBMFCoeffs_,
//                    time_
                    mesh_
                )
            );
        }
*/

    //- Destructor
    virtual ~linearForcesMotion();


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
