/*---------------------------------------------------------------------------*\
| Portions 2015 Copyright (C) Marc O. Woolliscroft                            |
\*---------------------------------------------------------------------------*/
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
    Foam::solidBodyMotionFunctions::propellerCrashbackMotionFunction

Description
    Motion function used to describe the rotation of an appendage from
    a body-fixed frame of reference. The rotation is defined by an origin
    and axis of rotation and an angular speed.

SourceFiles
    propellerCrashbackMotionFunction.C

\*---------------------------------------------------------------------------*/

#ifndef propellerCrashbackMotionFunction_H
#define propellerCrashbackMotionFunction_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"
#include "DataEntry.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class propellerCrashbackMotionFunction Declaration
\*---------------------------------------------------------------------------*/

class propellerCrashbackMotionFunction
:
    public solidBodyMotionFunction
{
    // Private data

        //- Origin of the axis
        const vector origin_;

        //- Axis vector
        const vector axis_;

        //- Angular velocity (rad/sec)
//        autoPtr<DataEntry<scalar> > omega_;

        //- Time that the propeller rotates at full RPM
        scalar trFullAhead_;

        //- Time for prop to go from full ahead RPM to zero
        //- (Also time for going from zero to full reverse)
        scalar trSlowDown_;

        //- Angular velocity in rad/sec
        scalar omega_;

    // Private Member Functions

        //- Disallow copy construct
        propellerCrashbackMotionFunction(const propellerCrashbackMotionFunction&);

        //- Disallow default bitwise assignment
        void operator=(const propellerCrashbackMotionFunction&);


public:

    //- Runtime type information
    TypeName("propellerCrashbackMotionFunction");


    // Constructors

        //- Construct from components
        propellerCrashbackMotionFunction
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new propellerCrashbackMotionFunction
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~propellerCrashbackMotionFunction();


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
