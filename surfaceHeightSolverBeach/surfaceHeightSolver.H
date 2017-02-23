/*---------------------------------------------------------------------------*\
Copyright (C) 2010 Engys Ltd.
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
    Foam::surfaceHeightSolver

Description
    Class for storing and calculating free surface height.

SourceFiles
    surfaceHeightSolver.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceHeightSolver_H
#define surfaceHeightSolver_H

#include "autoPtr.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class surfaceHeightSolver Declaration
\*---------------------------------------------------------------------------*/

class surfaceHeightSolver
:
    public IOdictionary
{
    // Private data
    
        //- reference to volume mesh
        const fvMesh& mesh_;
        
        //- volumetric height field
        volScalarField zeta_;

        //- surface mesh
        autoPtr<dynamicFvMesh> sMeshPtr_;

        dynamicFvMesh& sMesh_;

        //- full surface velocity
//        volVectorField Ufree_;

        //- surface height field
        volScalarField eta_;
        
        //- height field max displacement
        autoPtr<scalar> etaLimit_;
        
        //- solution switch
        Switch executeSolver_;

        //- surface advection velocity 
        vector Usurface_;

        //- ramp time
//        scalar tr_;

        //- surface final velocity
//        vector Ufinal_;

        // - beach coefficients
            
            //- numerical beach coefficients
            Switch beachActive_;
            
            //- list of fa patch names around which waves are of interest, 
            //- e.g. names of ship hull patches
            wordList featureNames_;
            
            //- distance from mid point of feature to mid point of beach
            //- as a multiple of feature size
            scalar beachDistance_;
            
            //- width of beach as a multiple of feature size
            scalar beachWidth_;
            
            //- derived beach properties
                //- feature mid-point
                vector featureCentroid_;

                //- feature size 
                scalar featureSize_;
        
    
    // Private Member Functions

        //- calculate beach properties
        void initBeach();

        //- apply beach to height field
        void beach();
        
        //- Disallow default bitwise copy construct and assignment
        surfaceHeightSolver(const surfaceHeightSolver&);
        void operator=(const surfaceHeightSolver&);

protected:

        static label solutionIndex_;

public:

    //- Runtime type information
    TypeName("surfaceHeightSolver");

    // Constructors

        //- Construct from components
        surfaceHeightSolver
        (
            const fvMesh& mesh,
            const word& regionName
        );

    //- Destructor
    virtual ~surfaceHeightSolver()
    {}


    // Member Functions

        //- update the height field
        bool correct(const surfaceScalarField& phi, const volVectorField& U);
        
        //- Access
            //- return height field for patch
            const scalarField& height(label) const;
                  
        //- Read the solver input dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
