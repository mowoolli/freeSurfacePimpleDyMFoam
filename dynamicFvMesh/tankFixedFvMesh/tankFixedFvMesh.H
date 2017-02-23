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
    Foam::tankFixedFvMesh

Description
    Solid-body motion of the mesh specified by 3DoF rigid body motion in the 
    horizontal plane. Ideal for replicating captive motion PMM maneuvers.

SourceFiles
    tankFixedFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef tankFixedFvMesh_H
#define tankFixedFvMesh_H

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
                           Class tankFixedFvMesh Declaration
\*---------------------------------------------------------------------------*/

class tankFixedFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Dictionary of motion control parameters
        dictionary dynamicMeshCoeffs_;

        //- The reference points which are transformed
        pointIOField undisplacedPoints_;

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

        //- Disallow default bitwise copy construct
        tankFixedFvMesh(const tankFixedFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const tankFixedFvMesh&);

        //- If the accel file has not been created create it
        void makeFiles();


public:

    //- Runtime type information
    TypeName("tankFixedFvMesh");


    // Constructors

        //- Construct from IOobject
        tankFixedFvMesh(const IOobject& io);


    // Destructor

        ~tankFixedFvMesh();


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
