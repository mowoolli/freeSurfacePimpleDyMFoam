/*---------------------------------------------------------------------------*\
| Portions 2015 Copyright (C) Marc O. Woolliscroft                            |
\*---------------------------------------------------------------------------*/
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

Description
    A dynamicFvMesh for replicating PMM maneuvers in the horizontal plane.
    Motions of surge, sway, and yaw are prescribed by the user, and a ramp
    option is available to avoid impulsive starts.

\*---------------------------------------------------------------------------*/

#include "tankFixedFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "transformField.H"
#include "Tuple2.H"
#include "mathematicalConstants.H"
#include "fvMesh.H"

using namespace Foam::constant::mathematical;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tankFixedFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, tankFixedFvMesh, IOobject);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
// Create files to record body motion
void Foam::tankFixedFvMesh::makeFiles()
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

Foam::tankFixedFvMesh::tankFixedFvMesh(const IOobject& io)
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
    positionFilePtr_(NULL),
    velocityFilePtr_(NULL),
    accelFilePtr_(NULL)
{
    read(dynamicMeshCoeffs_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tankFixedFvMesh::~tankFixedFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::tankFixedFvMesh::update()
{
    static bool hasWarned = false;

    // generate transformation

    scalar t = time().value();

    vector displacement(0, 0, 0);

    vector rotation(0, 0, 0);

    if (t < tr_ )
    {
        displacement = Ufinal_*0.5*(t - tr_/(pi + VSMALL)*sin(pi*t/(tr_ + VSMALL))) + t/tr_*swayAmp_*sin(omega_*t);
        rotation = t/tr_*yawAmp_*sin(omega_*t + yawPhase_)*pi/180;
    }
    else
    {
        displacement = Ufinal_*(t - tr_*0.5) + swayAmp_*sin(omega_*t);
        rotation = yawAmp_*sin(omega_*t + yawPhase_)*pi/180;
    }

    quaternion R(rotation.x(), rotation.y(), rotation.z());
    septernion TR(septernion(CofG_+displacement)*R*septernion(-CofG_));
   
    Info<< "tankFixedFvMesh::update(): "
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

        WarningIn("tankFixedFvMesh::update()")
            << "Did not find volVectorField U."
            << " Not updating U boundary conditions." << endl;
    }

    return true;
}

// read in rigid body properties from a dictionary
bool Foam::tankFixedFvMesh::read(const dictionary& coeffs)
{
    coeffs.lookup("tr") >> tr_;
    coeffs.lookup("Ufinal") >> Ufinal_;
    coeffs.lookup("swayAmp") >> swayAmp_;
    coeffs.lookup("yawAmp") >> yawAmp_;
    coeffs.lookup("yawPhase") >> yawPhase_;
    coeffs.lookup("omega") >> omega_;
    coeffs.lookup("CofG") >> CofG_;
    return true;
}


// ************************************************************************* //
