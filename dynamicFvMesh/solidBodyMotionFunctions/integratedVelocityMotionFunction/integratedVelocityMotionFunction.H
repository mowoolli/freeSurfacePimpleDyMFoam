/*---------------------------------------------------------------------------*\
| Portions 2015 Copyright (C) Marc O. Woolliscroft                            |
\*---------------------------------------------------------------------------*/
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

Class
    Foam::solidBodyMotionFunctions::integratedVelocityMotionFunction

Description
    For constant-speed PMM maneuvers. Surge velocity is varied to coincide
    with sway velocity, resulting in a constant magnitude of velocity. Very
    similar to tankFixedFvMesh. However, this being a 
    solidBodyMotionFunction and not a dynamicFvMesh, it is made for appended
    hullforms using the propellerMotionFunction for rotating propellers.
    Therefore, it uses the multiSolidBodyMotionFvMesh as the dynamicFvMesh option.

SourceFiles
    integratedVelocityMotionFunction.C

\*---------------------------------------------------------------------------*/

#ifndef integratedVelocityMotionFunction_H
#define integratedVelocityMotionFunction_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "forces.H"
#include "uniformDimensionedFields.H"
#include "Vector2D.H"
#include <string>
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class integratedVelocityMotionFunction Declaration
\*---------------------------------------------------------------------------*/

class integratedVelocityMotionFunction
:
    public solidBodyMotionFunction 
{
    // Private data

        //- time to reach desired speed
        scalar tr_;
        
        //- desired speed
        vector Ufinal_;

        //- amplitude of sway motion
        vector swayAmp_;

        //- amplitude of yaw motion
        vector yawAmp_;

        //- phase of yaw motion
        scalar yawPhase_;

        //- frequency of sway and yaw motions
        scalar omega_;

        //- center of gravity
        vector CofG_;

        //- position file ptr
        autoPtr<OFstream> positionFilePtr_;

        //- velocity file ptr
        autoPtr<OFstream> velocityFilePtr_;

        //- acceleration file ptr
        autoPtr<OFstream> accelFilePtr_;


    // Private Member Functions

        //- Disallow copy construct
        integratedVelocityMotionFunction(const integratedVelocityMotionFunction&);

        //- Disallow default bitwise assignment
        void operator=(const integratedVelocityMotionFunction&);

        //- If the accel file has not been created create it
        void makeFiles();


public:

    //- Runtime type information
    TypeName("integratedVelocityMotionFunction");


    // Constructors

        //- Construct from components
        integratedVelocityMotionFunction
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new integratedVelocityMotionFunction 
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }

    // Destructor

        virtual ~integratedVelocityMotionFunction();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
