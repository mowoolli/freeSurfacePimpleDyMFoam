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

\*---------------------------------------------------------------------------*/

#include "integratedVelocityMotionFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "Tuple2.H"
#include "IFstream.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(integratedVelocityMotionFunction, 0);
    addToRunTimeSelectionTable(solidBodyMotionFunction, integratedVelocityMotionFunction, dictionary);
};
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::solidBodyMotionFunctions::integratedVelocityMotionFunction::makeFiles()
{
    // Create the position file if not already created
    if (positionFilePtr_.empty())
    {
        if (debug)
        {
            Info<< "Creating position file." << endl;
        }

        // File update
        if (Pstream::master())
        {
            fileName positionDir;
            word startTimeName =
                time_.timeName(time_.startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                positionDir = time_.path()/".."/"RBmotion";
            }
            else
            {
                positionDir = time_.path()/"RBmotion";
            }

            // Create directory if does not exist.
            mkDir(positionDir);

            // Open new file at start up
            positionFilePtr_.reset(new OFstream(positionDir/"position.dat"));
        }
    }

    // Create the velocity file if not already created
    if (velocityFilePtr_.empty())
    {
        if (debug)
        {
            Info<< "Creating velocity file." << endl;
        }

        // File update
        if (Pstream::master())
        {
            fileName velocityDir;
            word startTimeName =
                time_.timeName(time_.startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                velocityDir = time_.path()/".."/"RBmotion";
            }
            else
            {
                velocityDir = time_.path()/"RBmotion";
            }

            // Create directory if does not exist.
            mkDir(velocityDir);

            // Open new file at start up
            velocityFilePtr_.reset(new OFstream(velocityDir/"velocity.dat"));
        }
    }

    // Create the accel file if not already created
    if (accelFilePtr_.empty())
    {
        if (debug)
        {
            Info<< "Creating acceleration file." << endl;
        }

        // File update
        if (Pstream::master())
        {
            fileName accelDir;
            word startTimeName =
                time_.timeName(time_.startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                accelDir = time_.path()/".."/"RBmotion";
            }
            else
            {
                accelDir = time_.path()/"RBmotion";
            }

            // Create directory if does not exist.
            mkDir(accelDir);

            // Open new file at start up
            accelFilePtr_.reset(new OFstream(accelDir/"accel.dat"));
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::integratedVelocityMotionFunction::integratedVelocityMotionFunction
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
    //const fvMesh& mesh
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    //solidBodyMotionFunction(SBMFCoeffs, mesh),
    positionFilePtr_(NULL),
    velocityFilePtr_(NULL),
    accelFilePtr_(NULL)
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::integratedVelocityMotionFunction::~integratedVelocityMotionFunction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion Foam::solidBodyMotionFunctions::integratedVelocityMotionFunction::transformation() const
{
//    scalar dt = time_.deltaT().value();
    scalar t = time_.value();

    vector rotation(0, 0, 0);

    vector displacement(0, 0, 0);

    int N = 1000;
    int i;

    scalar yawMag;
    yawMag = (yawAmp_ & (yawAmp_/mag(yawAmp_)))*pi/180;

    if (t < tr_ )
    {
        for(i = 0; i <= N; i++)
        {
            if (i == 0)
            {
                displacement += t/N*0.5*(Ufinal_*0.5*(1 - cos(pi*(t-t)/(tr_ + VSMALL)))*cos(-yawMag*cos(omega_*(t-t)))) \
                             + t/tr_*(swayAmp_/mag(swayAmp_))*t/N*0.5*(mag(Ufinal_)*sin(-yawMag*cos(omega_*(t-t))));
            }
            else if (i == N)
            {
                displacement += t/N*0.5*(Ufinal_*0.5*(1 - cos(pi*t/(tr_ + VSMALL)))*cos(-yawMag*cos(omega_*t))) \
                             + t/tr_*(swayAmp_/mag(swayAmp_))*t/N*0.5*(mag(Ufinal_)*sin(-yawMag*cos(omega_*t)));
            }
            else
            {
                displacement += t/N*(Ufinal_*0.5*(1 - cos(pi*(t-(N-i)*t/N)/(tr_ + VSMALL)))*cos(-yawMag*cos(omega_*(t-(N-i)*t/N)))) \
                             + t/tr_*(swayAmp_/mag(swayAmp_))*t/N*(mag(Ufinal_)*sin(-yawMag*cos(omega_*(t-(N-i)*t/N))));
            }
        }

        rotation = t/tr_*yawAmp_*sin(omega_*t + yawPhase_)*pi/180.;
    }
    else
    {
        for(i = 0; i <= N; i++)
        {
            if (i == 0)
            {
                displacement += tr_/N*0.5*(Ufinal_*0.5*(1 - cos(pi*(tr_-(N-i)*tr_/N)/(tr_ + VSMALL)))*cos(-yawMag*cos(omega_*(tr_-(N-i)*tr_/N)))) \
                             + (swayAmp_/mag(swayAmp_))*tr_/N*0.5*(mag(Ufinal_)*sin(-yawMag*cos(omega_*(tr_-(N-i)*tr_/N)))) \
                             + (t-tr_)/N*0.5*(Ufinal_*cos(-yawMag*cos(omega_*(tr_)))) \
                             + (swayAmp_/mag(swayAmp_))*(t-tr_)/N*0.5*(mag(Ufinal_)*sin(-yawMag*cos(omega_*(tr_))));
            }
            else if (i == N)
            {
                displacement += tr_/N*0.5*(Ufinal_*0.5*(1 - cos(pi*(tr_-(N-i)*tr_/N)/(tr_ + VSMALL)))*cos(-yawMag*cos(omega_*(tr_-(N-i)*tr_/N)))) \
                             + (swayAmp_/mag(swayAmp_))*tr_/N*0.5*(mag(Ufinal_)*sin(-yawMag*cos(omega_*(tr_-(N-i)*tr_/N)))) \
                             + (t-tr_)/N*0.5*(Ufinal_*cos(-yawMag*cos(omega_*t))) \
                             + (swayAmp_/mag(swayAmp_))*(t-tr_)/N*0.5*(mag(Ufinal_)*sin(-yawMag*cos(omega_*t)));
            }
            else
            {
                displacement += tr_/N*(Ufinal_*0.5*(1 - cos(pi*(tr_-(N-i)*tr_/N)/(tr_ + VSMALL)))*cos(-yawMag*cos(omega_*(tr_-(N-i)*tr_/N)))) \
                             + (swayAmp_/mag(swayAmp_))*tr_/N*(mag(Ufinal_)*sin(-yawMag*cos(omega_*(tr_-(N-i)*tr_/N)))) \
                             + (t-tr_)/N*(Ufinal_*cos(-yawMag*cos(omega_*(t-(N-i)*(t-tr_)/N)))) \
                             + (swayAmp_/mag(swayAmp_))*(t-tr_)/N*(mag(Ufinal_)*sin(-yawMag*cos(omega_*(t-(N-i)*(t-tr_)/N))));
            }
        }

        rotation = yawAmp_*sin(omega_*t + yawPhase_)*pi/180.;
    }

    quaternion R(rotation.x(), rotation.y(), rotation.z());
    septernion TR(septernion(CofG_+displacement)*R*septernion(-CofG_));

    Info<< "solidBodyMotionFunctions::integratedVelocityMotionFunction::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::integratedVelocityMotionFunction::read(const dictionary& SBMFCoeffs)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("tr") >> tr_;
    SBMFCoeffs_.lookup("Ufinal") >> Ufinal_;
    SBMFCoeffs_.lookup("swayAmp") >> swayAmp_;
    SBMFCoeffs_.lookup("yawAmp") >> yawAmp_;
    SBMFCoeffs_.lookup("yawPhase") >> yawPhase_;
    SBMFCoeffs_.lookup("omega") >> omega_;
    SBMFCoeffs_.lookup("CofG") >> CofG_;
    
    return true;
}


// ************************************************************************* //
