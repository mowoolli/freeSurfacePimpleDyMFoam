/*---------------------------------------------------------------------------*\
Copyright (C) 2010 Engys Ltd.
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

#include "surfaceHeightSolver.H"
#include "mappedPatchBase.H"
#include "bound.H"
#include "fvcMeshPhi.H"
#include "fvMesh.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(surfaceHeightSolver, 0);

label surfaceHeightSolver::solutionIndex_ = -1;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void surfaceHeightSolver::initBeach()
{
    if (beachActive_)
    {
    }
}

void surfaceHeightSolver::beach()
{
    scalar t = time().value();

    scalar displacement(0);

    if (t < 9 )
    {
        displacement = -2.556*0.5*(t - 9/(pi + VSMALL)*sin(pi*t/(9 + VSMALL)));
    }
    else
    {
        displacement = -2.556*(t - 9*0.5);
    }

    if (beachActive_)
    {
        scalarField distFunc
        (
            (
                (beachDistance_ - displacement)*5
                - mag(sMesh_.C().internalField().component(0) - 0)
            )
            *(6/(5*beachWidth_))
        );

        eta_.internalField() *= 0.5*(1 + tanh(distFunc));
        eta_.correctBoundaryConditions();
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surfaceHeightSolver::surfaceHeightSolver
(
    const fvMesh& mesh,
    const word& region
)
:
    IOdictionary
    (
        IOobject
        (
            "freeSurfaceProperties",
            mesh.time().constant(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),

    zeta_
    (
        IOobject
        (
            "zeta",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zeta", dimLength, 0.0) 
    ),

    sMeshPtr_
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                region,
                mesh.time().timeName(),
                mesh.time(),
                IOobject::MUST_READ
            )
        )
    ),
    sMesh_(sMeshPtr_()),

    eta_
    (
        IOobject
        (
            "eta",
            mesh.time().timeName(),
            sMesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sMesh_
    ),
    etaLimit_(),
    executeSolver_(this->lookupOrDefault<Switch>("solve", true)),
    Usurface_(this->lookup("Ufs")),
    beachActive_(this->lookupOrDefault<Switch>("beach", false)),
    featureNames_
    (
        this->lookupOrDefault<wordList>("beachPatches", wordList(0))
    ),
    beachDistance_(this->lookupOrDefault<scalar>("beachDistance", 1.5)),
    beachWidth_(this->lookupOrDefault<scalar>("beachWidth", 0.5)),
    featureCentroid_(vector::zero),
    featureSize_(0.0)
{
    if (this->found("etaLimit"))
    {
        etaLimit_.reset(new scalar(readScalar(this->lookup("etaLimit"))));
    }

    initBeach();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool surfaceHeightSolver::correct(const surfaceScalarField& phi, const volVectorField& U)
{
    if (solutionIndex_ == this->db().time().timeIndex() || !executeSolver_)
//    if (!executeSolver_)
    {
        return false;
    }
    else
    {
           sMesh_.update();
           dimensionedScalar omega_("omega", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0);
           volVectorField Vel_(sMesh_.C()*omega_);
            surfaceScalarField phis
        (
            IOobject
            (
                "phis",
                sMesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            linearInterpolate(Vel_) & sMesh_.Sf()
        );

        phis *= 0.0;

        if (phis.mesh().moving())
        {
            fvc::makeRelative(phis,Vel_);
        }

        //cell height field used to make eta relative to cell height
        volVectorField cellHeight
        (
            IOobject
            (
                "Hc",
                mesh_.time().timeName(),
                sMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sMesh_,
            dimensionedVector
            (
                "Hc", dimless, vector::zero
            ),
            zeroGradientFvPatchVectorField::typeName
        );

        //experimental non-linear velocity
        volVectorField Usv
        (
            IOobject
            (
                "Usv",
                mesh_.time().timeName(),
                sMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sMesh_,
            dimensionedVector
            (
                "Usv", dimVelocity, vector::zero
            ),
            zeroGradientFvPatchScalarField::typeName
        );
        
        //correct phis boundaries from mappedPatches
        //if (false)
        {
            // Since we're inside initEvaluate/evaluate there might be processor
            // comms underway. Change the tag we use.
            int oldTag = UPstream::msgType();
            UPstream::msgType() = oldTag + 1;

            // loop over all surface boundary patches
            forAll(sMesh_.boundary(), pi)
            {
                if (isA<mappedPatchBase>(sMesh_.boundary()[pi].patch()))
                {
            
                    // Get the scheduling information from the mappedPatchBase
                    const mappedPatchBase& mpp =
                        refCast<const mappedPatchBase>
                        (sMesh_.boundary()[pi].patch());
                        
                    const label samplePatchI = mpp.samplePolyPatch().index();

                    vectorField Unbr =
                        vectorField(U.boundaryField()[samplePatchI]);

                    // Swap to obtain full local values
                    mpp.distribute(Unbr);
                
                    //correct phis
                    phis.boundaryField()[pi]
                        = (Unbr & sMesh_.boundary()[pi].Sf());

                    //corret eta base patch
                    const labelList patchCells
                        (sMesh_.boundary()[pi].faceCells());
                    const scalarField Vol(sMesh_.V().field());

                    forAll(patchCells, fi)
                    {
                        label pci = patchCells[fi];
                        //eta_.boundaryField()[pi][fi] 
                        //    = Vol[pci]/ sMesh_.boundary()[pi].magSf()[fi];
                        //etaSource[pci] = Un[fi];
                        Usv[pci] = Unbr[fi];
                    }
                }
            }

            // Restore tag
            UPstream::msgType() = oldTag;
        }

        label nCorr = 1;
        for(label corr = 0; corr < nCorr; corr++)
        {

	        fvScalarMatrix etaEqn
	        (
//	            fvm::div(phis,eta_) + fvm::SuSp(fvc::div(phis), eta_)
                  fvm::ddt(eta_) + fvm::div(phis,eta_)
	        );

	        etaEqn.relax();
            etaEqn.solve();
        }

        //- limit maximum height deviation
        scalar maxMagEta(gMax(mag(eta_.internalField())));

        if (etaLimit_ < maxMagEta)
        {
            dimensionedScalar etaMin
                ("etaMmin", eta_.dimensions(), -etaLimit_);
            dimensionedScalar etaMax("etaMax", eta_.dimensions(), etaLimit_);
//            boundMinMax(eta_, etaMin, etaMax);
        }
        else
        {
            //Just report min and max eta
            Info << "Surface elevation, eta - min: "
                 << gMin(eta_) << ", max: " << gMax(eta_) << endl;
        }

        //beach
        beach();

        //apply eta internal field to zeta boundary field for all mapped patches
        {
            // Since we're inside initEvaluate/evaluate there might be processor
            // comms underway. Change the tag we use.
            int oldTag = UPstream::msgType();
            UPstream::msgType() = oldTag + 1;

            // loop over all volume mesh boundary patches
            forAll(mesh_.boundary(), pi)
            {
                if (isA<mappedPatchBase>(mesh_.boundary()[pi].patch()))
                {
            
                    // Get the scheduling information from the mappedPatchBase
                    const mappedPatchBase& mpp =
                        refCast<const mappedPatchBase>
                        (mesh_.boundary()[pi].patch());
                        
                    const label samplePatchI = mpp.samplePolyPatch().index();

                    const fvPatchScalarField& etaPatch =
                        refCast<const fvPatchScalarField>
                        (eta_.boundaryField()[samplePatchI]);

                    // Swap to obtain full local values
                    scalarField etaIntFld(etaPatch.patchInternalField());
                    mpp.distribute(etaIntFld);
                    
                    //correct zeta
                    zeta_.boundaryField()[pi] = etaIntFld;
                }
            }

            // Restore tag
            UPstream::msgType() = oldTag;
        }
     
        solutionIndex_ = mesh_.time().timeIndex();
        
        if (sMesh_.time().outputTime() )//&& debug)
        {
                Usv.write();
                phis.write();
        }
        
        if (mesh_.time().outputTime() || mesh_.time().timeIndex() == 0 )
        {
            Info << "Writing phis to file." << endl;
            phis.write();
            
            Info << "Divergence of face fluxes." << endl;
            volScalarField divPhis("divPhis", fvc::div(phis));
            divPhis.write();

            Info << "Writing eta to file." << endl;
            eta_.write();
        }

        return true;    
    }

}

const scalarField& surfaceHeightSolver::height(label patchi) const
{
    return zeta_.boundaryField()[patchi];
}

bool surfaceHeightSolver::read()
{
    Usurface_ = vector(this->lookup("Ufs"));

    beachActive_ = this->lookupOrDefault<Switch>("beach", false);

    featureNames_ =
    (
        this->lookupOrDefault<wordList>("featurePatches", wordList(0))
    );

    beachDistance_ = this->lookupOrDefault<scalar>("beachDistance", 1.5);

    beachWidth_ = this->lookupOrDefault<scalar>("beachWidth", 0.5);

    if (this->found("etaLimit"))
    {
        etaLimit_.reset(new scalar(readScalar(this->lookup("etaLimit"))));
    }
    else if (etaLimit_.valid())
    {
        etaLimit_.clear();
    }

    initBeach();

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
