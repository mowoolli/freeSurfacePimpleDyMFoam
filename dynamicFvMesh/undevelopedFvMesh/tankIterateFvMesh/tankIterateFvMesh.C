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

#include "tankIterateFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "transformField.H"
#include "Tuple2.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tankIterateFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, tankIterateFvMesh, IOobject);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
// Create files to record body motion
void Foam::tankIterateFvMesh::makeFiles()
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
                time().timeName(time().startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                positionDir = time().path()/".."/"RBmotion";
            }
            else
            {
                positionDir = time().path()/"RBmotion";
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
                time().timeName(time().startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                velocityDir = time().path()/".."/"RBmotion";
            }
            else
            {
                velocityDir = time().path()/"RBmotion";
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
                time().timeName(time().startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                accelDir = time().path()/".."/"RBmotion";
            }
            else
            {
                accelDir = time().path()/"RBmotion";
            }

            // Create directory if does not exist.
            mkDir(accelDir);

            // Open new file at start up
            accelFilePtr_.reset(new OFstream(accelDir/"accel.dat"));
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tankIterateFvMesh::tankIterateFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                io.time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    undisplacedPoints_
    (
        IOobject
        (
            "points",
            io.time().constant(),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    g_
    (
        IOobject
        (
            "g",
            time().constant(),
            time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    dt_(time().deltaT().value()),
    dt0_(time().deltaT().value()),
    curTimeIndex_(-1),
    positionFilePtr_(NULL),
    velocityFilePtr_(NULL),
    accelFilePtr_(NULL)
{
    read(dynamicMeshCoeffs_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tankIterateFvMesh::~tankIterateFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::tankIterateFvMesh::update()
{
    static bool hasWarned = false;
    
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
    if (curTimeIndex_ != time().timeIndex())
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

    dt_ = time().deltaT().value();
    
    // generate forces dictionary and calculate forces and moments
    dictionary forcesDict = dynamicMeshCoeffs_.subDict("forces");
    forcesDict.remove("CofR");
    forcesDict.add("CofR",(CofG_+x0_));
    
    forces f
    (
        "forces",
        *this,
        forcesDict
    );
    
    f.calcForcesMoment();

    vector F = f.forceEff();  // fluid force
    vector M = f.momentEff();  // fluid moment about CofG
    
    PstreamBuffers buf(Pstream::defaultCommsType);

    scalar t = time().value();
#    include "business.H"
    // update current time index
    if (curTimeIndex_ != time().timeIndex())
        curTimeIndex_ = time().timeIndex();

/*
    if (t < 1.5*tr_ )
    {
        x_          *= 0.0;
        theta_      *= 0.0;
        x0_         *= 0.0;
        theta0_     *= 0.0;
        x00_        *= 0.0;
        theta00_    *= 0.0;
        v_          *= 0.0;
        omega_      *= 0.0;
        v0_         *= 0.0;
        omega0_     *= 0.0;
        v00_        *= 0.0;
        omega00_    *= 0.0;
    }
*/

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
    
    Info<< "tankIterateFvMesh::update(): "
        << "Time = " << t << " transformation: " << TR << endl;
    
    // transform the mesh
    fvMesh::movePoints
    (
        transform(TR, undisplacedPoints_)
    );

    // update velocity boundary conditions
    if (foundObject<volVectorField>("U"))
    {
        const_cast<volVectorField&>(lookupObject<volVectorField>("U"))
            .correctBoundaryConditions();
    }
    else if (!hasWarned)
    {
        hasWarned = true;

        WarningIn("tankIterateFvMesh::update()")
            << "Did not find volVectorField U."
            << " Not updating U boundary conditions." << endl;
    }

    return true;
}

// read in rigid body properties from a dictionary
bool Foam::tankIterateFvMesh::read(const dictionary& coeffs)
{
    coeffs.lookup("tr") >> tr_;
    coeffs.lookup("Ufinal") >> Ufinal_;
    coeffs.lookup("mass") >> mass_;
    coeffs.lookup("CofG") >> CofG_;
    coeffs.lookup("I") >> I_;
    coeffs.lookup("x0") >> x0_;
    coeffs.lookup("v0") >> v0_;
    coeffs.lookup("theta0") >> theta0_;
    coeffs.lookup("omega0") >> omega0_;
    x_ = x0_;
    v_ = v0_;
    theta_ = theta0_;
    omega_ = omega0_;
    x00_ = coeffs.lookupOrDefault<vector>("x00", x0_);
    v00_ = coeffs.lookupOrDefault<vector>("v00", v0_);
    theta00_ = coeffs.lookupOrDefault<vector>("theta00", theta0_);
    omega00_ = coeffs.lookupOrDefault<vector>("omega00", omega0_);
    
    m_a_ = coeffs.lookupOrDefault<vector>("m_a", vector::zero);
    I_a_ = coeffs.lookupOrDefault<vector>("I_a", vector::zero);
    
    coeffs.lookup("translationFree") >> translationFree_;
    coeffs.lookup("rotationFree") >> rotationFree_;
    
    coeffs.lookup("k") >> k_;
    coeffs.lookup("K_t") >> K_;
    coeffs.lookup("c") >> c_;
    coeffs.lookup("C_t") >> C_;
    
    coeffs.lookup("relaxation") >> beta_;
    
    order_ = coeffs.lookupOrDefault<label>("order", 1);
    
    return true;
}


// ************************************************************************* //
