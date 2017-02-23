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

#include "rigidBody3DoFMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "Tuple2.H"
#include "IFstream.H"
#include "surfaceRegionMotion.H"
#include "OFstream.H"
#include "IOReferencer.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(rigidBody3DoFMotionSolver, 0);
    addToRunTimeSelectionTable(PiroSolidBodyMotionFunction, rigidBody3DoFMotionSolver, dictionary);
};
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::solidBodyMotionFunctions::rigidBody3DoFMotionSolver::makeFiles()
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

Foam::solidBodyMotionFunctions::rigidBody3DoFMotionSolver::rigidBody3DoFMotionSolver
(
    const dictionary& SBMFCoeffs,
    //const Time& runTime
    const fvMesh& mesh
)
:
    //PiroSolidBodyMotionFunction(SBMFCoeffs, runTime),
    PiroSolidBodyMotionFunction(SBMFCoeffs, mesh),
    g_
    (
        IOobject
        (
            "g",
            //runTime.constant(),
            //runTime,
            mesh.time().constant(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    curTimeIndex_(-1),
    positionFilePtr_(NULL),
    velocityFilePtr_(NULL),
    accelFilePtr_(NULL)
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::rigidBody3DoFMotionSolver::~rigidBody3DoFMotionSolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion Foam::solidBodyMotionFunctions::rigidBody3DoFMotionSolver::transformation()
{
    scalar dt = time_.deltaT().value();
    scalar t = time_.value();

    vector accel_est = ( vNew_ - vOld_ ) / dt;
    vector alpha_est = ( omegaNew_ - omegaOld_ ) / dt;
    
    if (curTimeIndex_ != time_.timeIndex())
    {
        xOld_ = xNew_;
        vOld_ = vNew_;
        thetaOld_ = thetaNew_;
        omegaOld_ = omegaNew_;
    }
    
    dictionary forcesDict = SBMFCoeffs_.subDict("forces");
    forcesDict.remove("CofR");
    forcesDict.add("CofR",(CofG_+xOld_));
    
    forces f
    (
        "forces",
        mesh_,
        //SBMFCoeffs_.subDict("forces")
        forcesDict
    );
    
//    forces::forcesMoments fm = f.calcForcesMoment();
    f.calcForcesMoment();
//    vector F = fm.first().first() + fm.first().second();  // fluid force
    vector F = f.forceEff();
    Info << "f.forceEff() = " << F << endl;
//    vector M = (fm.second().first() + fm.second().second());  // fluid moment about CofG
    vector M = f.momentEff();
    Info << "f.momentEff() = " << M << endl;
    
    // rotate into body coordinate system
    quaternion rot(thetaNew_.x(), thetaNew_.y(), thetaNew_.z());
    M = rot.invTransform(M);
    
    // gravitational force
    F += mass_*g_.value();
    
    // springs and dampers
    F.x() -= (k_.x()*xOld_.x()+c_.x()*vOld_.x());
    F.y() -= (k_.y()*xOld_.y()+c_.y()*vOld_.y());
    F.z() -= (k_.z()*xOld_.z()+c_.z()*vOld_.z());
    M.x() -= (K_.x()*thetaOld_.x()+C_.x()*omegaOld_.x());
    M.y() -= (K_.y()*thetaOld_.y()+C_.y()*omegaOld_.y());
    M.z() -= (K_.z()*thetaOld_.z()+C_.z()*omegaOld_.z());
    
    // estimates of added mass
    F.x() += m_a_.x()*accel_est.x();
    F.y() += m_a_.y()*accel_est.y();
    F.z() += m_a_.z()*accel_est.z();    
    M.x() += I_a_.x()*alpha_est.x();
    M.y() += I_a_.y()*alpha_est.y();
    M.z() += I_a_.z()*alpha_est.z();
    
    // accelerations
    vector acceleration(vector::zero);
    vector alpha = inv( I_ + symmTensor(I_a_.x(), 0, 0, I_a_.y(), 0, I_a_.z() ) ) & M;

    if (t < tr_)
    {
//        vector acceleration(vector::zero);
        acceleration.x() = F.x() / ( mass_ + m_a_.x() );
        acceleration.y() = F.y() / ( mass_ + m_a_.y() );
        acceleration.z() = F.z() / ( mass_ + m_a_.z() );
        acceleration.x() *= translationFree_.x();
        acceleration.y() *= 0;
        acceleration.z() *= 0;
//        vector alpha = inv( I_ + symmTensor(I_a_.x(), 0, 0, I_a_.y(), 0, I_a_.z() ) ) & M;
        alpha.x() *= rotationFree_.x();
        alpha.y() *= rotationFree_.y();
        alpha.z() *= 0; 
    }
    else
    {
//        vector acceleration(vector::zero);
        acceleration.x() = F.x() / ( mass_ + m_a_.x() );
        acceleration.y() = F.y() / ( mass_ + m_a_.y() );
        acceleration.z() = F.z() / ( mass_ + m_a_.z() );
        acceleration.x() *= translationFree_.x();
        acceleration.y() *= translationFree_.y();
        acceleration.z() *= translationFree_.z();
//        vector alpha = inv( I_ + symmTensor(I_a_.x(), 0, 0, I_a_.y(), 0, I_a_.z() ) ) & M;
        alpha.x() *= rotationFree_.x();    
        alpha.y() *= rotationFree_.y();    
        alpha.z() *= rotationFree_.z();
    }

    // velocities updated explicitly
    vector vTilde = vOld_ + acceleration*dt;
    vector omegaTilde = omegaOld_ + alpha*dt;
    
    // positions updated implicity
    vector xTilde = xOld_ + vTilde*dt;
    vector thetaTilde = thetaOld_ + omegaTilde*dt;
    
    // relax updates
    vNew_ = beta_*vTilde + (1 - beta_)*vNew_;
    omegaNew_ = beta_*omegaTilde + (1 - beta_)*omegaNew_;
    xNew_ = beta_*xTilde + (1 - beta_)*xNew_;
    thetaNew_ = beta_*thetaTilde + (1 - beta_)*thetaNew_;
    
    quaternion R(thetaNew_.x(), thetaNew_.y(), thetaNew_.z());
    septernion TR(septernion(CofG_ + xNew_)*R*septernion(-CofG_));
    
//    scalar t = time_.value();
    Info<< "solidBodyMotionFunctions::rigidBody3DoFMotionSolver::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    // Reference the septernion dictionary from createFields.H
    IOdictionary& septDict = const_cast<IOdictionary&>(mesh_.lookupObject<IOdictionary>("septDict"));

    // Set the value to the septernion (TR) calculated above so the surface mesh region can access and use the same septernion
    septDict.set("TRpassed",TR);

    makeFiles();

    if (Pstream::master())
    {
        positionFilePtr_() << t << tab << xNew_ << tab << thetaNew_ << endl;
        velocityFilePtr_() << t << tab << vNew_ << tab << omegaNew_ << endl;
//        accelFilePtr_() << t << tab << acceleration << tab << alpha << endl;
        accelFilePtr_() << t << tab << (vNew_-vOld_)/dt << tab << (omegaNew_-omegaOld_)/dt << endl;
 
    if (curTimeIndex_ != time_.timeIndex())
        curTimeIndex_ = time_.timeIndex();
    }

    return TR;
}

bool Foam::solidBodyMotionFunctions::rigidBody3DoFMotionSolver::read(const dictionary& SBMFCoeffs)
{
    PiroSolidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("tr") >> tr_;
    SBMFCoeffs_.lookup("mass") >> mass_;
    SBMFCoeffs_.lookup("CofG") >> CofG_;
    SBMFCoeffs_.lookup("I") >> I_;
    SBMFCoeffs_.lookup("x0") >> xOld_;
    SBMFCoeffs_.lookup("v0") >> vOld_;
    SBMFCoeffs_.lookup("theta0") >> thetaOld_;
    SBMFCoeffs_.lookup("omega0") >> omegaOld_;
    xNew_ = xOld_;
    vNew_ = vOld_;
    thetaNew_ = thetaOld_;
    omegaNew_ = omegaOld_;
    
    m_a_ = SBMFCoeffs_.lookupOrDefault<vector>("m_a", vector::zero);
    I_a_ = SBMFCoeffs_.lookupOrDefault<vector>("I_a", vector::zero);
    
    SBMFCoeffs_.lookup("translationFree") >> translationFree_;
    SBMFCoeffs_.lookup("rotationFree") >> rotationFree_;
    
    SBMFCoeffs_.lookup("k") >> k_;
    SBMFCoeffs_.lookup("K_t") >> K_;
    SBMFCoeffs_.lookup("c") >> c_;
    SBMFCoeffs_.lookup("C_t") >> C_;
    
    SBMFCoeffs_.lookup("relaxation") >> beta_;
    
    return true;
}


// ************************************************************************* //
