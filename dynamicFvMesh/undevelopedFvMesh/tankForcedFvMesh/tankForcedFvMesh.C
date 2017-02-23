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

#include "tankForcedFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "transformField.H"
#include "Tuple2.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tankForcedFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, tankForcedFvMesh, IOobject);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
// Create files to record body motion
void Foam::tankForcedFvMesh::makeFiles()
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

Foam::tankForcedFvMesh::tankForcedFvMesh(const IOobject& io)
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

Foam::tankForcedFvMesh::~tankForcedFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::tankForcedFvMesh::update()
{
    static bool hasWarned = false;

    // generate transformation

    scalar t = time().value();

    vector displacement(0, 0, 0);

    vector eulerAngles(0, 0, 0);

    if (t < tr_ )
    {
        displacement = Ufinal_*0.5*(t - tr_/(pi + VSMALL)*sin(pi*t/(tr_ + VSMALL)));
    }
    else
    {
        displacement = Ufinal_*(t - tr_*0.5);

    	eulerAngles = forcedAmplitude_*sin(forcedOmega_*t);

    	// Convert the rotational motion from deg to rad
    	eulerAngles *= pi/180.0;
    }
    
    //quaternion R(theta_.x(), theta_.y(), theta_.z());
    //quaternion R(0, 0, 0);
    quaternion R(eulerAngles.x(), eulerAngles.y(), eulerAngles.z());
    septernion TR(septernion(CofG_ + displacement)*R*septernion(-CofG_));
    
    Info<< "tankForcedFvMesh::update(): "
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

        WarningIn("tankForcedFvMesh::update()")
            << "Did not find volVectorField U."
            << " Not updating U boundary conditions." << endl;
    }

    return true;
}

// read in rigid body properties from a dictionary
bool Foam::tankForcedFvMesh::read(const dictionary& coeffs)
{
    coeffs.lookup("tr") >> tr_;
    coeffs.lookup("Ufinal") >> Ufinal_;
    coeffs.lookup("forcedAmplitude") >> forcedAmplitude_;
    coeffs.lookup("forcedOmega") >> forcedOmega_;
    coeffs.lookup("CofG") >> CofG_;
    return true;
}


// ************************************************************************* //
