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
    Solve for 6DoF rigid body motion based on forces from forces object.
    Able to constrain any degree of freedom, and add springs and dampers
    to any degree of freedom, though no coupling between degrees of freedom
    is done.

See Also
    SDA (Ship design analysis) for 3DoF motion.

SourceFiles
    integratedVelocityMotionFunction.C

\*---------------------------------------------------------------------------*/

#ifndef integratedVelocityMotionFunction_H
#define integratedVelocityMotionFunction_H

//#include "PiroSolidBodyMotionFunction.H"
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

/*
        //- mass of body
        scalar mass_;
        
        //- estimate of added mass of body (diagonal components only)
        vector m_a_;
        
        //- center of gravity
        vector CofG_;
        
        //- moments of inertia about center of gravity
        symmTensor I_;
        
        //- estimate of added intertia of body (diagonal components only)
        vector I_a_;
        
        //- gravity
        uniformDimensionedVectorField g_;
        
        //- new displacement vector
        vector xNew_;
        
        //- new velocity vector
        vector vNew_;
        
        //- new rotation vector
        vector thetaNew_;
        
        //- new rotational velocity vector
        vector omegaNew_;

        //- old displacement vector
        vector xOld_;
        
        //- old velocity vector
        vector vOld_;
        
        //- old rotation vector
        vector thetaOld_;
        
        //- old rotational velocity vector
        vector omegaOld_;
        
        //- vector of free translational dof - 1 if free, 0 if fixed
        vector translationFree_;
        
        //- vector of free rotational dof - 1 if free, 0 if fixed
        vector rotationFree_;

        //- vector of translational spring coefficients
        vector k_;
        
        //- vector of rotational spring coefficients
        vector K_;
        
        //- vector of translational damping coefficients
        vector c_;
        
        //- vector of rotational damping coefficients
        vector C_;
        
        //- relaxation parameter (should be between 0 and 1)
        scalar beta_;
        
        //- current time index - to determine if new time level
        label curTimeIndex_;
*/

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
            //const fvMesh& mesh
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
