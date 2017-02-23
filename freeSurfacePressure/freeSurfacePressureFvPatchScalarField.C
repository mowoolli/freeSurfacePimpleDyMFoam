/*---------------------------------------------------------------------------*\
| Portions 2015 Copyright (C) Marc O. Woolliscroft                            |
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "freeSurfacePressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "surfaceHeightSolver.H"
#include "uniformDimensionedFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Memember Functions * * * * * * * * * * * //

void Foam::freeSurfacePressureFvPatchScalarField::initFreeSurfaceSolver()
{
    //- surface mesh 
    if (solveFreeSurface_)
    {
        if 
        (
            !this->db().time().foundObject<surfaceHeightSolver>
            ("freeSurfaceProperties")
        )
        {
            //get the other-side region
            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>(patch().patch());
            
        
            autoPtr<surfaceHeightSolver> shs
            (
                new surfaceHeightSolver
                (
                    patch().boundaryMesh().mesh(),
                    mpp.sampleRegion()
                )
            );
            shs->store(shs);
        }
    }
    else
    {
        //instantiate zeta field directly
        if 
        (
            !this->db().foundObject<volScalarField>("zeta")
        )
        {
            autoPtr<volScalarField> zeta
            (
                new volScalarField
                (
                    IOobject
                    (
                        "zeta",
                        db().time().timeName(),
                        db(),
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    patch().boundaryMesh().mesh(),
                    dimensionedScalar("zeta", dimLength, 0.0) 
                )
            );
            zeta->store(zeta);
        }
        
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freeSurfacePressureFvPatchScalarField::freeSurfacePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi"),
    rhoName_("rho"),
    p0_(p.size(), 0.0),
    solveFreeSurface_(false)
{}


Foam::freeSurfacePressureFvPatchScalarField::freeSurfacePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    p0_(p.size(), 0.0),
    solveFreeSurface_(dict.lookupOrDefault<Switch>("solveFreeSurface", true))
{
    if (dict.found("p0"))
    {
        p0_ = scalarField("p0", dict, p.size());
    }

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(p0_);
    }
    
    initFreeSurfaceSolver();
}


Foam::freeSurfacePressureFvPatchScalarField::freeSurfacePressureFvPatchScalarField
(
    const freeSurfacePressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    p0_(ptf.p0_, mapper),
    solveFreeSurface_(ptf.solveFreeSurface_)
{
    initFreeSurfaceSolver();
}


Foam::freeSurfacePressureFvPatchScalarField::freeSurfacePressureFvPatchScalarField
(
    const freeSurfacePressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    p0_(tppsf.p0_),
    solveFreeSurface_(tppsf.solveFreeSurface_)
{
    initFreeSurfaceSolver();
}


Foam::freeSurfacePressureFvPatchScalarField::freeSurfacePressureFvPatchScalarField
(
    const freeSurfacePressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    p0_(tppsf.p0_),
    solveFreeSurface_(tppsf.solveFreeSurface_)
{
    initFreeSurfaceSolver();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freeSurfacePressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    p0_.autoMap(m);
}


void Foam::freeSurfacePressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const freeSurfacePressureFvPatchScalarField& tiptf =
        refCast<const freeSurfacePressureFvPatchScalarField>(ptf);

    p0_.rmap(tiptf.p0_, addr);
}


void Foam::freeSurfacePressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const surfaceScalarField& phi = 
        db().lookupObject<surfaceScalarField>(phiName_);
		
    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    scalarField Upn = phip / patch().magSf();

    scalarField zetap(patch().size(), 0.0);

    const uniformDimensionedVectorField& g =
        db().lookupObject<uniformDimensionedVectorField>("g");

    scalarField gNorm(patch().nf() & g.value());

    if (solveFreeSurface_)
    {
        surfaceHeightSolver& fshs
        (
            const_cast<surfaceHeightSolver&>
            (
                this->db().time().lookupObject<surfaceHeightSolver>
                ("freeSurfaceProperties")
            )
        );

        fshs.correct(db().lookupObject<surfaceScalarField>(phiName_), db().lookupObject<volVectorField>(UName_));

        zetap = fshs.height(patch().index());
    }
    else 
    {
        volScalarField& zeta
        (
            const_cast<volScalarField&>
            (
                db().lookupObject<volScalarField>("zeta")
            )
        );
        
        zeta.storePrevIter();
        
        Info << "YOU ARE IN THE ELSE STATEMENT!!!" << endl;

        zeta.boundaryField()[patch().index()]
            -= 0.5*mag(Upn)*Upn/gNorm;
            
        zeta.relax();
    
        zetap = zeta.boundaryField()[patch().index()];
    }
    
    if (phi.dimensions() == dimVelocity*dimArea)
    {
        operator==
		(
                      p0_ - gNorm*zetap
		);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        FatalErrorIn
        (
            "freeSurfacePressureFvPatchScalarField::updateCoeffs()"
        )   << "dimensions of " << phiName_ << " are not correct"
            << abort(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();

}


void Foam::freeSurfacePressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    p0_.writeEntry("p0", os);
    writeEntryIfDifferent<Switch>
    (
        os, "solveFreeSurface", true, solveFreeSurface_
    );
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        freeSurfacePressureFvPatchScalarField
    );
}

// ************************************************************************* //
