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
    Foam::rigidBody6DoFMotionFvMesh

Description
    Solid-body motion of the mesh specified by 6DoF rigid body motion based 
    on forces from forces object.  Able to constrain any degree of freedom, 
    and add springs and dampers to any degree of freedom, though no coupling 
    between degrees of freedom is done.

SourceFiles
    rigidBody6DoFMotionFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef rigidBody6DoFMotionFvMesh_H
#define rigidBody6DoFMotionFvMesh_H

#include "dynamicFvMesh.H"
#include "primitiveFields.H"
#include "dictionary.H"
#include "pointIOField.H"
#include "forces.H"
#include "uniformDimensionedFields.H"
#include "Vector2D.H"
#include <string>
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class rigidBody6DoFMotionFvMesh Declaration
\*---------------------------------------------------------------------------*/

class rigidBody6DoFMotionFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Dictionary of motion control parameters
        dictionary dynamicMeshCoeffs_;

        //- The reference points which are transformed
        pointIOField undisplacedPoints_;

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
        vector x_;
        
        //- new velocity vector
        vector v_;
        
        //- new rotation vector
        vector theta_;
        
        //- new rotational velocity vector
        vector omega_;

        //- old displacement vector
        vector x0_;
        
        //- old velocity vector
        vector v0_;
        
        //- old rotation vector
        vector theta0_;
        
        //- old rotational velocity vector
        vector omega0_;

        //- old-old displacement vector
        vector x00_;
        
        //- old-old velocity vector
        vector v00_;
        
        //- old-old rotation vector
        vector theta00_;
        
        //- old-old rotational velocity vector
        vector omega00_;
        
        //- current timestep
        scalar dt_;
        
        //- old timestep
        scalar dt0_;
        
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
        
        //- order of accuracy for time derivative
        label order_;
        
        //- current time index - to determine if new time level
        label curTimeIndex_;

        //- position file ptr
        autoPtr<OFstream> positionFilePtr_;

        //- velocity file ptr
        autoPtr<OFstream> velocityFilePtr_;

        //- acceleration file ptr
        autoPtr<OFstream> accelFilePtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        rigidBody6DoFMotionFvMesh(const rigidBody6DoFMotionFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const rigidBody6DoFMotionFvMesh&);

        //- If the accel file has not been created create it
        void makeFiles();


public:

    //- Runtime type information
    TypeName("rigidBody6DoFMotionFvMesh");


    // Constructors

        //- Construct from IOobject
        rigidBody6DoFMotionFvMesh(const IOobject& io);


    // Destructor

        ~rigidBody6DoFMotionFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update properties from given dictionary
        virtual bool read(const dictionary& coeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
