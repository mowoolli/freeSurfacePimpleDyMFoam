/*---------------------------------------------------------------------------*\ 
| Modified 2010-2012 Copyright (C) Engys Ltd                                  |
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    ct_(0.0),
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "freeSurfaceHeightFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "steadyStateDdtScheme.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * *  Private member functions * * * * * * * * * * * //

void freeSurfaceHeightFvPatchScalarField::updateCoeffs(const volVectorField& ufs)
{
    label i = patch().index();
    const fvPatchVectorField& up(ufs.boundaryField()[i]);

    vectorField Uspi(up.patchInternalField());

//    Info << "1. up = " << up.patchInternalField() << endl;
//    Info << "1. Uspi = " << Uspi << endl;

//    vectorField Hcpi
//    (
//        db().lookupObject<volVectorField>("Hc")
//        .boundaryField()[i].patchInternalField()
//    );
//    Hcpi /= max(SMALL, mag(Hcpi));

    //remove plane normal component
//    Uspi -= (Uspi & Hcpi)*Hcpi;

            const scalar t = db().time().value();
//            dimensionedScalar t("t", dimTime, time().value());
            scalar tr_ = 4;
//            dimensionedScalar tr_("tr", dimTime, 1);
            vector Ufinal_(1.5, 0, 0);
//            dimensionedVector Ufinal_("Ufinal", dimVelocity, vector(2, 0, 0));
            vector swayAmp_(0, 0, 0);
//            dimensionedVector swayAmp_("swayAmp", dimLength, vector(0, 0, 0));
            scalar omega_ = 0.0;
//            dimensionedScalar omega_("omega", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0);
            vector yawAmp_(0, 0, 0);
            scalar yawPhase_ = 0;
            vector CofG_(0, 0, 0);
//            dimensionedVector CofG_("CofG", dimLength, vector(0, 0, 0));
            vector yawAxisNormal_ = yawAmp_/(mag(yawAmp_) + VSMALL);

            if (t < tr_ )
            {
                forAll(this->patch().Cf(), i)
                {
                    CofG_ -= (CofG_ & yawAxisNormal_)*yawAxisNormal_;
                    vector faceCenter_ = this->patch().Cf()[i];
                    faceCenter_ -= (faceCenter_ & yawAxisNormal_)*yawAxisNormal_;

//                    Info << "faceCenter = " << faceCenter_ << endl;

                    Uspi = (-0.5*Ufinal_*(cos((pi*t)/(tr_ + VSMALL)) - 1.0)) + (t*swayAmp_*omega_*cos(omega_*t) + swayAmp_*sin(omega_*t))/tr_ + (((t*omega_*mag(yawAmp_)*cos                       (omega_*t + yawPhase_) + (mag(yawAmp_)*sin(omega_*t + yawPhase_))/tr_)*yawAmp_/(mag(yawAmp_) + VSMALL)) ^ (faceCenter_ - CofG_));

                }
            }
            else
            {   
                forAll(this->patch().Cf(), i)
                {
                    CofG_ -= (CofG_ & yawAxisNormal_)*yawAxisNormal_;
                    vector faceCenter_ = this->patch().Cf()[i];
                    faceCenter_ -= (faceCenter_ & yawAxisNormal_)*yawAxisNormal_;
                    
//                    Info << "faceCenter = " << faceCenter_ << endl;

                    Uspi = Ufinal_ + swayAmp_*omega_*cos(omega_*t) + ((mag(yawAmp_)*cos(omega_*t + yawPhase_)*yawAmp_/(mag(yawAmp_) + VSMALL)) ^ (faceCenter_ - CofG_));
                }
            }

    //normalise
    Uspi /= max(SMALL, mag(Uspi));

//    Info << "ct_ = " << ct_ << endl;
//    Info << "cp_ = " << cp_ << endl;
//    Info << "cm_ = " << cm_ << endl;

    this->valueFraction() = -patch().nf() & Uspi;
//    Info << "valueFraction() = " << valueFraction() << endl;
    this->valueFraction() = cm_*max(0.0, 1/(1-ct_)*(this->valueFraction()-ct_));
//    Info << "valueFraction() = " << valueFraction() << endl;
    this->valueFraction() = pow(this->valueFraction(), cp_);
//    Info << "valueFraction() = " << valueFraction() << endl;
//    Info << "refValue() = " << refValue() << endl;

    fvPatchScalarField::operator=(valueFraction()*refValue() + (1.0 - valueFraction())*(this->patchInternalField() + refGrad()/this->patch().deltaCoeffs()));
//    Info << "Hello from inside updateCoeffs" << endl;
//    Info << "operator = " << (valueFraction()*refValue() + (1.0 - valueFraction())*(this->patchInternalField() + refGrad()/this->patch().deltaCoeffs())) << endl;
//    Info << "1. waterHeight_ = " << waterHeight_ << endl;
//    Info << "1. refValue() = " << refValue() << endl;
//    Info << "1. refGrad() = " << refGrad() << endl;
//    Info << "1. valueFraction() = " << valueFraction() << endl;
//    Info << "2. Uspi = " << Uspi << endl;
//    Info << "1. nf() = " << patch().nf() << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

freeSurfaceHeightFvPatchScalarField::freeSurfaceHeightFvPatchScalarField
(
	const fvPatch& p,
	const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    UsName_(""),
    waterHeight_(0.0),
    cm_(1.0),
    ct_(0.0),
    cp_(1.0)
{
    fvPatchScalarField::operator= (pTraits<scalar>::zero);
    this->refValue() = pTraits<scalar>::zero ;
    this->refGrad() = pTraits<scalar>::zero ;
    this->valueFraction() = 0.0 ;

//    Info << "Hello from 1st constructor" << endl;
}

freeSurfaceHeightFvPatchScalarField::freeSurfaceHeightFvPatchScalarField
(
	const freeSurfaceHeightFvPatchScalarField& ptf,
	const fvPatch& p,
	const DimensionedField<scalar, volMesh>& iF,
	const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    UsName_(ptf.UsName_),
    waterHeight_(ptf.waterHeight_),
    cm_(ptf.cm_),
    ct_(ptf.ct_),
    cp_(ptf.cp_)
{
//    Info << "Hello from 2nd constructor" << endl;
}

freeSurfaceHeightFvPatchScalarField::freeSurfaceHeightFvPatchScalarField
(
	const fvPatch& p,
	const DimensionedField<scalar, volMesh>& iF,
	const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    UsName_(dict.lookupOrDefault<word>("Usurface", "Usv")),
    waterHeight_(dict.lookupOrDefault<scalar>("waterHeight",0.0)),
    cm_(dict.lookupOrDefault<scalar>("CtranMul",1.0)),
    ct_(dict.lookupOrDefault<scalar>("CtranOff",0.0)),
    cp_(dict.lookupOrDefault<scalar>("CtranPow",1.0))
{
//    Info << "Hello from 3rd constructor" << endl;    
    this->refGrad() = 0.0;
    this->refValue() = waterHeight_;
    this->valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    }

    //might be in initialisation
    if (db().foundObject<volVectorField>("Usv") 
        && db().foundObject<volVectorField>("Hc"))
    {
//        Info << "Hello from const 3 before updateCoeffs" << endl;
        updateCoeffs(db().lookupObject<volVectorField>("Usv"));
//        Info << "Hello from const 3 after updateCoeffs" << endl;
//        Info << "refValue() = " << refValue() << endl;
        mixedFvPatchField<scalar>::evaluate();
    }

//    Info << "waterHeight_ = " << waterHeight_ << endl;
//    Info << "refValue() = " << refValue() << endl;
//    Info << "refGrad() = " << refGrad() << endl;
//    Info << "valueFraction() = " << valueFraction() << endl;
//    Info << "UsName_ = " << UsName_ << endl;
}

freeSurfaceHeightFvPatchScalarField::freeSurfaceHeightFvPatchScalarField
(
  const freeSurfaceHeightFvPatchScalarField& ptf
)
	:
    mixedFvPatchField<scalar>(ptf),
    UsName_(ptf.UsName_),
    waterHeight_(ptf.waterHeight_),
    cm_(ptf.cm_),
    ct_(ptf.ct_),
    cp_(ptf.cp_)
{
//    Info << "Hello from 4th constructor" << endl;
}

freeSurfaceHeightFvPatchScalarField::freeSurfaceHeightFvPatchScalarField
(
  const freeSurfaceHeightFvPatchScalarField& ptf,
  const DimensionedField<scalar, volMesh>& iF
)
	:
    mixedFvPatchField<scalar>(ptf, iF),
    UsName_(ptf.UsName_),
    waterHeight_(ptf.waterHeight_),
    cm_(ptf.cm_),
    ct_(ptf.ct_),
    cp_(ptf.cp_)
{
//    Info << "Hello from 5th constructor" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void freeSurfaceHeightFvPatchScalarField::updateCoeffs()
{
    if(this->updated()) 
    {
        return;
    }
//    Info << "Hello from updateCoeffs 1" << endl;
    updateCoeffs(db().lookupObject<volVectorField>(UsName_));
//    Info << "Hello from updateCoeffs 2" << endl;
//    mixedFvPatchField<scalar>::updateCoeffs();
//    Info << "Hello from updateCoeffs 3" << endl;
}

void freeSurfaceHeightFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os) ;
    writeEntryIfDifferent<scalar>(os, "waterHeight", 0.0, waterHeight_);
    writeEntryIfDifferent<word>(os, "Us", "Usv", UsName_);
    this->writeEntry("value", os) ;
}

makePatchTypeField(fvPatchScalarField, freeSurfaceHeightFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
