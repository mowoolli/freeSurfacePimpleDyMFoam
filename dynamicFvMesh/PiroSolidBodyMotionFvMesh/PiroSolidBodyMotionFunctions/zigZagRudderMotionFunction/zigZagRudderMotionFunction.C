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

\*---------------------------------------------------------------------------*/

#include "zigZagRudderMotionFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "scalar.H"
#include <string>
//#include "OFstream.H"
#include "IFstream.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(zigZagRudderMotionFunction, 0);
    addToRunTimeSelectionTable
    (
        PiroSolidBodyMotionFunction,
        zigZagRudderMotionFunction,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::solidBodyMotionFunctions::zigZagRudderMotionFunction::makeFiles()
{
    // Create the heading angle file if not already created
    if (headingAngleFilePtr_.empty())
    {
        if (debug)
        {
            Info<< "Creating heading angle file." << endl;
        }

        // File update
        if (Pstream::master())
        {
            fileName headingAngleDir;
            word startTimeName =
                time_.timeName(time_.startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                headingAngleDir = time_.path()/".."/"overshootData";
            }
            else
            {
                headingAngleDir = time_.path()/"overshootData";
            }

            // Create directory if does not exist.
            mkDir(headingAngleDir);

            // Open new file at start up
            headingAngleFilePtr_.reset(new OFstream(headingAngleDir/"headingAngle.dat"));
        }
    }

    // Create the rudder angle file if not already created
    if (rudderAngleFilePtr_.empty())
    {
        if (debug)
        {
            Info<< "Creating rudder angle file." << endl;
        }

        // File update
        if (Pstream::master())
        {
            fileName rudderAngleDir;
            word startTimeName =
                time_.timeName(time_.startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                rudderAngleDir = time_.path()/".."/"overshootData";
            }
            else
            {
                rudderAngleDir = time_.path()/"overshootData";
            }

            // Create directory if does not exist.
            mkDir(rudderAngleDir);

            // Open new file at start up
            rudderAngleFilePtr_.reset(new OFstream(rudderAngleDir/"rudderAngle.dat"));
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::zigZagRudderMotionFunction::zigZagRudderMotionFunction
(
    const dictionary& SBMFCoeffs,
//    const Time& runTime
    const fvMesh& mesh
)
:
//    solidBodyMotionFunction(SBMFCoeffs, runTime),
    PiroSolidBodyMotionFunction(SBMFCoeffs, mesh),
    origin_(SBMFCoeffs_.lookup("origin")),
    axis_(SBMFCoeffs_.lookup("axis")),
    omega_(DataEntry<scalar>::New("omega", SBMFCoeffs_)),
    headingAngleFilePtr_(NULL),
    rudderAngleFilePtr_(NULL)

{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::zigZagRudderMotionFunction::~zigZagRudderMotionFunction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::zigZagRudderMotionFunction::transformation()
{
//    OFstream outputFile("overshootData.dat");

    scalar t = time_.value();
    scalar dt = time_.deltaT().value();

    vector eulerAngles(0, 0, 0);
//    vector eulerAngles;
    vector ihat(1, 0, 0);
    vector jhat(0, 1, 0);
    vector khat(0, 0, 1);
    scalar psiMin = minHeadingAngle_*(pi/180);
    scalar psiMax = maxHeadingAngle_*(pi/180);
    scalar deltaMin = minRudderAngle_*(pi/180);
    scalar deltaMax = maxRudderAngle_*(pi/180);

    static scalar intStartTime;
    static scalar intStartTimeForYawCheck;
//    scalar intEndTime;

    static bool posToNeg;
    bool loopPop = false;

    //patch ID whose normal matches the heading angle
    const label patchID = mesh_.boundaryMesh().findPatchID("inlet");
    const vectorField & headingAngleField = mesh_.boundary()[patchID].nf();
    vector headingAngle = gSum(headingAngleField)/mag(gSum(headingAngleField));
    scalar psi = Foam::acos((headingAngle & referenceHeadingAngle_)/(mag(headingAngle)*mag(referenceHeadingAngle_)))*Foam::sign(headingAngle & jhat);

    Info << "gSum(headingAngleField) = " << gSum(headingAngleField) << endl;
//    Info << "mag(gSum(headingAngleField)) = " << mag(gSum(headingAngleField)) << endl;
//    Info << "headingAngle = " << headingAngle << endl;
//    Info << "psi = " << psi << endl;

    quaternion R;

    if ((t < trr_) || (t == trr_))
    {
        eulerAngles.x() *= 0;
        eulerAngles.y() *= 0;
        eulerAngles.z() *= 0;
        intStartTime = trr_;
//        quaternion R(eulerAngles/mag(eulerAngles + headingAngle), mag(eulerAngles));
        R = (eulerAngles/mag(eulerAngles + headingAngle), 1); // dummy value to initialize quaternion with zero rotation
    }
    else
    {
        Info << "eulerAngles = " << eulerAngles << endl;
        while (psiMin < psi && psi < psiMax && loopPop == false)
        {
            Info << "omega_->integrate(intStartTime,t) = " << omega_->integrate(intStartTime,t) << endl;
//            if ((eulerAngles.x() == 0) && (eulerAngles.y() == 0) && (eulerAngles.z() == 0))
            if ((t < trr_ + dt) || (t == trr_ + dt))
            {
                Info << "statement 1..." << endl;
                intStartTime = trr_;
//                intStartTimeForYawCheck = t;
                eulerAngles = -omega_->integrate(intStartTime,t)*axis_; // initial decrease in rudder angle after ramp
//                intEndTime = t;
                posToNeg = false;
                loopPop = true;
            }
            else if ((deltaMin < -omega_->integrate(intStartTime,t) && omega_->integrate(intStartTime,t) < deltaMax) && posToNeg)
            {
//                intEndTime = t;
                eulerAngles = omega_->integrate(intStartTime,t)*axis_; // keep increasing rudder angle
                loopPop = true;
                Info << "statement 2..." << endl; // WORK
            }
            else if ((omega_->integrate(intStartTime,t) > deltaMax || omega_->integrate(intStartTime,t) == deltaMax) && posToNeg)
            {
                eulerAngles = deltaMax*axis_; // maintain rudder angle at maximum deflection // POSSIBLE ERROR
                intStartTimeForYawCheck = t;
                loopPop = true;
                Info << "statement 3..." << endl;
            }
            else if ((deltaMax > omega_->integrate(intStartTime,t) && -omega_->integrate(intStartTime,t) > deltaMin) && !posToNeg)
            {
//                intEndTime = t;
                eulerAngles = -omega_->integrate(intStartTime,t)*axis_; // keep decreasing rudder angle
                loopPop = true;
                Info << "statement 4..." << endl; // WORK
            }
            else if ((-omega_->integrate(intStartTime,t) < deltaMin || -omega_->integrate(intStartTime,t) == deltaMin) && !posToNeg)
            {
                eulerAngles = deltaMin*axis_; // maintain rudder angle at minimum deflection // POSSIBLE ERROR
                intStartTimeForYawCheck = t;
                loopPop = true;
                Info << "intStartTime = " << intStartTime << endl;
                Info << "intStartTimeForYawCheck = " << intStartTimeForYawCheck << endl;
                Info << "statement 5..." << endl;
            }
        }
        while ((psi > psiMax || psi == psiMax) && loopPop == false)
        {
//            if (-omega_->integrate(intStartTime,t) < deltaMin || -omega_->integrate(intStartTime,t) == deltaMin)
            if (intStartTimeForYawCheck == t - dt)
            {
                Info << "intStartTime = " << intStartTime << endl;
                intStartTime = t;
                eulerAngles = deltaMin*axis_+omega_->integrate(intStartTimeForYawCheck,t)*axis_; // initiate increase in rudder angle
                posToNeg = true;
                loopPop = true;
                Info << "intStartTime = " << intStartTime << endl;
                Info << "intStartTimeForYawCheck = " << intStartTimeForYawCheck << endl;
                Info << "statement 6..." << endl;
            }
            else if((deltaMin+omega_->integrate(intStartTimeForYawCheck,t) > deltaMin) && (deltaMin+omega_->integrate(intStartTimeForYawCheck,t) < deltaMax))
            {
//                intEndTime = t;
                eulerAngles = deltaMin*axis_+omega_->integrate(intStartTimeForYawCheck,t)*axis_; // keep increasing rudder angle
                loopPop = true;
                Info << "intStartTime = " << intStartTime << endl;
                Info << "intStartTimeForYawCheck = " << intStartTimeForYawCheck << endl;
                Info << "statement 7..." << endl;
            }
            else
            {
                eulerAngles = deltaMax*axis_; // maintain rudder angle at maximum deflection
                loopPop = true;
//                intStartTime = t;
                Info << "intStartTime = " << intStartTime << endl;
                Info << "intStartTimeForYawCheck = " << intStartTimeForYawCheck << endl;
                Info << "statement 8..." << endl;
            }
        }
        while ((psi < psiMin || psi == psiMin) && loopPop == false)
        {
//            if (omega_->integrate(intStartTime,t) > deltaMax || omega_->integrate(intStartTime,t) == deltaMax)
            if (intStartTimeForYawCheck == t - dt)
            {
                intStartTime = t;
                eulerAngles = deltaMax*axis_-omega_->integrate(intStartTimeForYawCheck,t)*axis_; // initiate decrease in rudder angle
                posToNeg = false;
                loopPop = true;
                Info << "statement 9..." << endl;
            }
            else if ((deltaMax-omega_->integrate(intStartTimeForYawCheck,t) < deltaMax) && (deltaMax-omega_->integrate(intStartTimeForYawCheck,t) > deltaMin))
            {
//                intEndTime = t;
                eulerAngles = deltaMax*axis_-omega_->integrate(intStartTimeForYawCheck,t)*axis_; // keep decreasing rudder angle
                loopPop = true;
                Info << "intStartTime = " << intStartTime << endl;
                Info << "intStartTimeForYawCheck = " << intStartTimeForYawCheck << endl;
                Info << "statement 10..." << endl;
            }
            else
            {
                eulerAngles = deltaMin*axis_; // maintain rudder angle at minimum deflection
                loopPop = true;
//                intStartTime = t;
                Info << "intStartTime = " << intStartTime << endl;
                Info << "intStartTimeForYawCheck = " << intStartTimeForYawCheck << endl;
                Info << "statement 11..." << endl;
            }
        }
//        Info << "eulerAngles before R calculation = " << eulerAngles << endl;
//        Info << "mag(eulerAngles) = " << mag(eulerAngles) << endl;
        scalar theta = eulerAngles.z()/axis_.z();
//        quaternion R(eulerAngles/mag(eulerAngles), mag(eulerAngles));
//        R = (eulerAngles/mag(eulerAngles), mag(eulerAngles));
        R.v() = axis_.x()*Foam::sin(theta/2)*ihat + axis_.y()*Foam::sin(theta/2)*jhat + axis_.z()*Foam::sin(theta/2)*khat;
//        R.w() = Foam::cos(theta/2);
//        R.v() = eulerAngles/mag(eulerAngles);
        R.w() = Foam::cos(mag(eulerAngles)/2);
//        Info << "R = " << R << endl;
    }

//    outputFile << time_.value() << "  " << psi << "  " << headingAngle << "  " << (eulerAngles.z()/axis_.z())*(180/pi) << endl;

//    Info << "eulerAngles = " << eulerAngles << endl;
//    Info << "R = " << R << endl;
//    quaternion R(eulerAngles/mag(eulerAngles), mag(eulerAngles));
    septernion TR(septernion(origin_)*R*septernion(-origin_));
    Info << "TR = " << TR << endl;

    Info<< "solidBodyMotionFunctions::zigZagRudderMotionFunction::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    makeFiles();

    if (Pstream::master())
    {
        headingAngleFilePtr_() << t << tab << psi << tab << headingAngle << endl;
        rudderAngleFilePtr_() << t << tab << (eulerAngles.z()/axis_.z())*(180/pi) << endl;
    }

    return TR;
}


bool Foam::solidBodyMotionFunctions::zigZagRudderMotionFunction::read
(
    const dictionary& SBMFCoeffs
)
{
    PiroSolidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("trr") >> trr_;
    SBMFCoeffs_.lookup("minHeadingAngle") >> minHeadingAngle_;
    SBMFCoeffs_.lookup("maxHeadingAngle") >> maxHeadingAngle_;
    SBMFCoeffs_.lookup("referenceHeadingAngle") >> referenceHeadingAngle_;
    SBMFCoeffs_.lookup("minRudderAngle") >> minRudderAngle_;
    SBMFCoeffs_.lookup("maxRudderAngle") >> maxRudderAngle_;

    omega_.reset
    (
        DataEntry<scalar>::New("omega", SBMFCoeffs_).ptr()
    );

    return true;
}

// ************************************************************************* //
