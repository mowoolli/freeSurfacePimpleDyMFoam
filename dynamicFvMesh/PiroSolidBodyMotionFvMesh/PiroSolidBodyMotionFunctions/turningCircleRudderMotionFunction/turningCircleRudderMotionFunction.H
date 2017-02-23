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
    Foam::solidBodyMotionFunctions::turningCircleRudderMotionFunction

Description
    SolidBodyMotionFvMesh 6DoF motion function.

    The rotation is defined by an origin and axis of rotation and an angular
    speed.

SourceFiles
    turningCircleRudderMotionFunction.C

\*---------------------------------------------------------------------------*/

#ifndef turningCircleRudderMotionFunction_H
#define turningCircleRudderMotionFunction_H

#include "PiroSolidBodyMotionFunction.H"
//#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"
#include "DataEntry.H"
#include "autoPtr.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class turningCircleRudderMotionFunction Declaration
\*---------------------------------------------------------------------------*/

class turningCircleRudderMotionFunction
:
    public PiroSolidBodyMotionFunction
{
    // Private data

        //- Origin of the axis
        const vector origin_;

        //- Axis vector
        const vector axis_;

        //- Angular velocity (helm rate of rudder) in rad/sec
        autoPtr<DataEntry<scalar> > omega_;

//        IOobject psi;

//        IOobject eulerAngles;

        //- Ramp time. Should be greater than or equal to ramp time of zigZagMotionSolver
        scalar trr_;

//        scalar minHeadingAngle_;
        
//        scalar maxHeadingAngle_;
        
        vector referenceHeadingAngle_;

//        scalar minRudderAngle_;

        scalar maxRudderAngle_;

        //- position file ptr
        autoPtr<OFstream> headingAngleFilePtr_;

        //- velocity file ptr
        autoPtr<OFstream> rudderAngleFilePtr_;
        

    // Private Member Functions

        //- Disallow copy construct
        turningCircleRudderMotionFunction(const turningCircleRudderMotionFunction&);

        //- Disallow default bitwise assignment
        void operator=(const turningCircleRudderMotionFunction&);

        //- If the headingAngle file has not been created create it
        void makeFiles();


public:

    //- Runtime type information
    TypeName("turningCircleRudderMotionFunction");


    // Constructors

        //- Construct from components
        turningCircleRudderMotionFunction
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
                new turningCircleRudderMotionFunction
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }
*/

    //- Destructor
    virtual ~turningCircleRudderMotionFunction();


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
