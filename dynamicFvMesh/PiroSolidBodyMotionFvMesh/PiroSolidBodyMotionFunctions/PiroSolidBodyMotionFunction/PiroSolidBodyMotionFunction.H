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

Namespace
    Foam::PiroSolidBodyMotionFunctions

Description
    Namespace for solid-body motions


Class
    Foam::PiroSolidBodyMotionFunction

Description
    Base class for defining solid-body motions

SourceFiles
    PiroSolidBodyMotionFunction.C
    newDynamicFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef PiroSolidBodyMotionFunction_H
#define PiroSolidBodyMotionFunction_H

#include "Time.H"
#include "fvMesh.H"
#include "dictionary.H"
#include "septernion.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PiroSolidBodyMotionFunction Declaration
\*---------------------------------------------------------------------------*/

class PiroSolidBodyMotionFunction
{
protected:

    // Protected data

        dictionary SBMFCoeffs_;
        const fvMesh& mesh_;
        const Time& time_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        PiroSolidBodyMotionFunction(const PiroSolidBodyMotionFunction&);

        //- Disallow default bitwise assignment
        void operator=(const PiroSolidBodyMotionFunction&);


public:
//    const fvMesh& mesh_;

    //- Runtime type information
    TypeName("PiroSolidBodyMotionFunction");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            PiroSolidBodyMotionFunction,
            dictionary,
            //(const dictionary& SBMFCoeffs, const Time& runTime),
            //(SBMFCoeffs, runTime)
            (const dictionary& SBMFCoeffs, const fvMesh& mesh),
            (SBMFCoeffs, mesh)
        );


    // Constructors

        //- Construct from the SBMFCoeffs dictionary and Time
        PiroSolidBodyMotionFunction
        (
            const dictionary& SBMFCoeffs,
            //const Time& runTime
            const fvMesh& mesh
        );


    // Selectors

        //- Select constructed from the SBMFCoeffs dictionary and Time
        static autoPtr<PiroSolidBodyMotionFunction> New
        (
            const dictionary& SBMFCoeffs,
            //const Time& runTime
            const fvMesh& mesh
        );


    // Destructor

        virtual ~PiroSolidBodyMotionFunction();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() = 0;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
