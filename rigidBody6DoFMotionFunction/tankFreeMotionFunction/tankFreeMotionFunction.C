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

#include "tankFreeMotionFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "transformField.H"
#include "Tuple2.H"
#include "IFstream.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(tankFreeMotionFunction, 0);
    addToRunTimeSelectionTable(PiroSolidBodyMotionFunction, tankFreeMotionFunction, dictionary);
};
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::solidBodyMotionFunctions::tankFreeMotionFunction::makeFiles()
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

Foam::solidBodyMotionFunctions::tankFreeMotionFunction::tankFreeMotionFunction
(
    const dictionary& SBMFCoeffs,
//    const Time& runTime
    const fvMesh& mesh
)
:
//    solidBodyMotionFunction(SBMFCoeffs, runTime),
    PiroSolidBodyMotionFunction(SBMFCoeffs, mesh),
/*
    dynamicMeshCoeffs_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
//                io.time().constant(),
//                time_.constant(),
//                time_,
                mesh.time().constant(),
                mesh.time(),
//                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
*/
//    patches_(wordReList(dynamicMeshCoeffs_.lookup("patches"))),
//    patchSet_(mesh.boundaryMesh().patchSet(patches_)),
//    rhoInf_(1.0),
//    rhoName_(dynamicMeshCoeffs_.lookupOrDefault<word>("rhoName", "rho")),
    g_
    (
        IOobject
        (
            "g",
//            time_.constant(),
//            time_,
            mesh.time().constant(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

//    functionObjectLibs_(runTime, 0),
    dt_(time_.deltaT().value()),
    dt0_(time_.deltaT().value()),
    curTimeIndex_(-1),
    positionFilePtr_(NULL),
    velocityFilePtr_(NULL),
    accelFilePtr_(NULL)
{
    read(SBMFCoeffs);

//    if (rhoName_ == "rhoInf")
//    {
//        rhoInf_ = readScalar(dynamicMeshCoeffs_.lookup("rhoInf"));
//    }
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::tankFreeMotionFunction::~tankFreeMotionFunction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion Foam::solidBodyMotionFunctions::tankFreeMotionFunction::transformation()
{
//    static bool hasWarned = false;
    
    vector accel_est, alpha_est;
    
    if (order_ <= 1)
    {
        accel_est = ( v_ - v0_ ) / dt_;
        alpha_est = ( omega_ - omega0_ ) / dt_;
    }
    else
    {
        scalar a0 = (dt0_+2.0*dt_)/dt_/(dt0_+dt_);
        scalar a1 = -(dt0_+dt_)/dt0_/dt_;
        scalar a2 = dt_/dt0_/(dt0_+dt_);
        
        accel_est = a0*v_ + a1*v0_ + a2*v00_;
        alpha_est = a0*omega_ + a1*omega0_ + a2*omega00_;
    }
    
    // if new time-step, update old position and velocity
    if (curTimeIndex_ != time_.timeIndex())
    {
        x00_ = x0_;
        v00_ = v0_;
        theta00_ = theta0_;
        omega00_ = omega0_;
        x0_ = x_;
        v0_ = v_;
        theta0_ = theta_;
        omega0_ = omega_;
        dt0_ = dt_;
    }

    dt_ = time_.deltaT().value();
    
    // generate forces dictionary and calculate forces and moments
//    dictionary forcesDict = dynamicMeshCoeffs_.subDict("forces");
    dictionary forcesDict = SBMFCoeffs_.subDict("forces");
//    forcesDict.add("type", forces::typeName);
//    forcesDict.add("type", type_);
//    forcesDict.add("functionObjectLibs", functionObjectLibs_);
//    forcesDict.add("patches", patches_);
//    forcesDict.add("pName", p_);
//    forcesDict.add("UName", U_);
//    forcesDict.add("rhoName", rhoInf_);
//    forcesDict.add("rhoInf", rhoInf_);
//    forcesDict.add("unknown", unknown_);
//    forcesDict.add("cS", cS_);
    forcesDict.remove("CofR");
    forcesDict.add("CofR",(CofG_+x0_));
    
    forces f
    (
        "forces",
//        type_,
//        functionObjectLibs_,
//        patches_,
//        pName_,
//        UName_,
//        rhoName_,
//        rhoInf_,
//        unknown_,
//        cS_
//        this->db(),
//        *this,
        mesh_,
        forcesDict
    );
    
    f.calcForcesMoment();

    vector F = f.forceEff();  // fluid force
    Info << "f.forceEff() = " << F << endl;
    vector M = f.momentEff();  // fluid moment about CofG
    Info << "f.momentEff() = " << M << endl;

    PstreamBuffers buf(Pstream::defaultCommsType);
    if (Pstream::master())
    {
        // rotate moments into body coordinate system
        quaternion rot(theta_.x(), theta_.y(), theta_.z());
        M = rot.invTransform(M);
        
        // apply gravitational force
        F += mass_*g_.value();
        
        // apply spring and damper forces
        F.x() -= (k_.x()*x0_.x()+c_.x()*v0_.x());
        F.y() -= (k_.y()*x0_.y()+c_.y()*v0_.y());
        F.z() -= (k_.z()*x0_.z()+c_.z()*v0_.z());
        M.x() -= (K_.x()*theta0_.x()+C_.x()*omega0_.x());
        M.y() -= (K_.y()*theta0_.y()+C_.y()*omega0_.y());
        M.z() -= (K_.z()*theta0_.z()+C_.z()*omega0_.z());
        
        // apply estimates of added mass force
        F.x() += m_a_.x()*accel_est.x();
        F.y() += m_a_.y()*accel_est.y();
        F.z() += m_a_.z()*accel_est.z();    
        M.x() += I_a_.x()*alpha_est.x();
        M.y() += I_a_.y()*alpha_est.y();
        M.z() += I_a_.z()*alpha_est.z();
        
        // calculate new accelerations
        vector acceleration(vector::zero);
        acceleration.x() = F.x() / ( mass_ + m_a_.x() );
        acceleration.y() = F.y() / ( mass_ + m_a_.y() );
        acceleration.z() = F.z() / ( mass_ + m_a_.z() );
        acceleration.x() *= translationFree_.x();
        acceleration.y() *= translationFree_.y();
        acceleration.z() *= translationFree_.z();
        vector alpha = inv( I_ + symmTensor(I_a_.x(), 0, 0, I_a_.y(), 0, I_a_.z() ) ) & M;
        alpha.x() *= rotationFree_.x();    
        alpha.y() *= rotationFree_.y();    
        alpha.z() *= rotationFree_.z();
        
        vector vTilde, omegaTilde, xTilde, thetaTilde;
        if (order_ <= 1)
        {
            // velocities updated explicitly
            vTilde = v0_ + acceleration*dt_;
            omegaTilde = omega0_ + alpha*dt_;
            
            // positions updated implicity
            xTilde = x0_ + vTilde*dt_;
            thetaTilde = theta0_ + omegaTilde*dt_;
        }
        else
        {
            scalar a0 = (dt0_+2.0*dt_)/dt_/(dt0_+dt_);
            scalar a1 = -(dt0_+dt_)/dt0_/dt_;
            scalar a2 = dt_/dt0_/(dt0_+dt_);
            
            // velocities updated explicitly
            vTilde = (-a1*v0_ - a2*v00_ + acceleration)/a0;
            omegaTilde = (-a1*omega0_ - a2*omega00_ + alpha)/a0;
            
            // positions updated implicity
            xTilde = (-a1*x0_ - a2*x00_ + vTilde)/a0;
            thetaTilde = (-a1*theta0_ - a2*theta00_ + omegaTilde)/a0;
        }
        
        // explicitly relax updates
        v_ = beta_*vTilde + (1 - beta_)*v_;
        omega_ = beta_*omegaTilde + (1 - beta_)*omega_;
        x_ = beta_*xTilde + (1 - beta_)*x_;
        theta_ = beta_*thetaTilde + (1 - beta_)*theta_;
        
        // create the recording files if not already done
        makeFiles();
    
        // write out position, velocity, and acceleration
        scalar t = time_.value();
        positionFilePtr_() << t << tab << x_ << tab << theta_ << endl;
        velocityFilePtr_() << t << tab << v_ << tab << omega_ << endl;
        accelFilePtr_() << t << tab << (v_-v0_)/dt_ << tab << (omega_-omega0_)/dt_ << endl;
        
        // scatter output to different processors
        for (int ii = 0; ii < Pstream::nProcs(); ii++)
        {
            if (ii != Pstream::myProcNo())
            {
                UOPstream toNeighbor(ii, buf);
                vectorField out(2);
                out[0] = x_;
                out[1] = theta_;
                toNeighbor << out;
            }
        }
        buf.finishedSends();
    }
    else
    {
        buf.finishedSends();
        UIPstream fromMaster(Pstream::masterNo(), buf);
        vectorField in(2, vector::zero);
        
        fromMaster >> in;
        x_ = in[0];
        theta_ = in[1];
    }
    
    // update current time index
    if (curTimeIndex_ != time_.timeIndex())
        curTimeIndex_ = time_.timeIndex();

    scalar t = time_.value();


    if (t < tr_ )
    {
        x_[0] = Ufinal_[0]*0.5*(t - tr_/(pi + VSMALL)*sin(pi*t/(tr_ + VSMALL)));
    }
    else
    {
        x_[0] = Ufinal_[0]*(t - tr_*0.5);
    }

    // generate transformation
    quaternion R(theta_.x(), theta_.y(), theta_.z());
    septernion TR(septernion(CofG_ + x_)*R*septernion(-CofG_));
    
    Info<< "tankFreeFvMesh::update(): "
        << "Time = " << t << " transformation: " << TR << endl;
    
    // transform the mesh
//    fvMesh::movePoints
//    (
//        transform(TR, undisplacedPoints_)
//    );

    // update velocity boundary conditions
/*
    if (foundObject<volVectorField>("U"))
    {
        const_cast<volVectorField&>(lookupObject<volVectorField>("U"))
            .correctBoundaryConditions();
    }
    else if (!hasWarned)
    {
        hasWarned = true;

        WarningIn("tankFreeFvMesh::update()")
            << "Did not find volVectorField U."
            << " Not updating U boundary conditions." << endl;
    }
*/
    Info<< "solidBodyMotionFunctions::tankFreeMotionFunction::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::tankFreeMotionFunction::read(const dictionary& SBMFCoeffs)
{
    PiroSolidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("tr") >> tr_;
    SBMFCoeffs_.lookup("Ufinal") >> Ufinal_;
    SBMFCoeffs_.lookup("mass") >> mass_;
    SBMFCoeffs_.lookup("CofG") >> CofG_;
    SBMFCoeffs_.lookup("I") >> I_;
    SBMFCoeffs_.lookup("x0") >> x0_;
    SBMFCoeffs_.lookup("v0") >> v0_;
    SBMFCoeffs_.lookup("theta0") >> theta0_;
    SBMFCoeffs_.lookup("omega0") >> omega0_;
    x_ = x0_;
    v_ = v0_;
    theta_ = theta0_;
    omega_ = omega0_;
    x00_ = SBMFCoeffs_.lookupOrDefault<vector>("x00", x0_);
    v00_ = SBMFCoeffs_.lookupOrDefault<vector>("v00", v0_);
    theta00_ = SBMFCoeffs_.lookupOrDefault<vector>("theta00", theta0_);
    omega00_ = SBMFCoeffs_.lookupOrDefault<vector>("omega00", omega0_);
    
    m_a_ = SBMFCoeffs_.lookupOrDefault<vector>("m_a", vector::zero);
    I_a_ = SBMFCoeffs_.lookupOrDefault<vector>("I_a", vector::zero);
    
    SBMFCoeffs_.lookup("translationFree") >> translationFree_;
    SBMFCoeffs_.lookup("rotationFree") >> rotationFree_;
    
    SBMFCoeffs_.lookup("k") >> k_;
    SBMFCoeffs_.lookup("K_t") >> K_;
    SBMFCoeffs_.lookup("c") >> c_;
    SBMFCoeffs_.lookup("C_t") >> C_;
    
    SBMFCoeffs_.lookup("relaxation") >> beta_;
    
    order_ = SBMFCoeffs_.lookupOrDefault<label>("order", 1);

    
    return true;
}


// ************************************************************************* //
