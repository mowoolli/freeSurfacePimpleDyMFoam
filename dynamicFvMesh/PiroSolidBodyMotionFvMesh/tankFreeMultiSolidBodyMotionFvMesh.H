/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::tankFreeMultiSolidBodyMotionFvMesh

Description
    Solid-body motion of the mesh specified by a run-time selectable
    motion function.

SourceFiles
    tankFreeMultiSolidBodyMotionFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef tankFreeMultiSolidBodyMotionFvMesh_H
#define tankFreeMultiSolidBodyMotionFvMesh_H

#include "dynamicFvMesh.H"
#include "dictionary.H"
#include "pointIOField.H"
#include "PiroSolidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class tankFreeMultiSolidBodyMotionFvMesh Declaration
\*---------------------------------------------------------------------------*/

class tankFreeMultiSolidBodyMotionFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Dictionary of motion control parameters
        dictionary dynamicMeshCoeffs_;

        //- The motion control function
        PtrList<PiroSolidBodyMotionFunction> SBMFs_;

        //- The reference points which are transformed
        pointIOField undisplacedPoints_;

        //- Specified cellZones
        labelList zoneIDs_;

        //- Points to move per cellZone
        labelListList pointIDs_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tankFreeMultiSolidBodyMotionFvMesh(const tankFreeMultiSolidBodyMotionFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const tankFreeMultiSolidBodyMotionFvMesh&);


public:

    //- Runtime type information
    TypeName("tankFreeMultiSolidBodyMotionFvMesh");


    // Constructors

        //- Construct from IOobject
        tankFreeMultiSolidBodyMotionFvMesh(const IOobject& io);


    //- Destructor
    ~tankFreeMultiSolidBodyMotionFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
