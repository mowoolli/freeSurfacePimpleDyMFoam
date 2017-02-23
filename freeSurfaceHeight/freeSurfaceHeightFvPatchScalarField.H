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

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::freeSurfaceHeightFvPatchScalarField

Description
    Foam::freeSurfaceHeightFvPatchScalarField

SourceFiles
    freeSurfaceHeightFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef freeSurfaceHeightFvPatchScalarField_H
#define freeSurfaceHeightFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam 
{

/*---------------------------------------------------------------------------*\
                     Class freeSurfaceHeightFvPatch Declaration
\*---------------------------------------------------------------------------*/

class freeSurfaceHeightFvPatchScalarField
:
	public mixedFvPatchField<scalar> 
{

protected:

    // Protected data

    //- surface velocity name
        word UsName_;

    //- water height
        scalar waterHeight_;

    //- transom coefficient a: VF *= cm*cosTheta
        scalar cm_;
        
    //- transom coefficient b: VF = max(0, cm*(cosTheta/(1-ct) - ct/(1-ct)))
        scalar ct_;

    //- transom coefficient c: VF = (VF)^b (applied after transformations ab)
        scalar cp_;

private:

    void updateCoeffs(const volVectorField&);

public:

    //- Runtime scalar information
    TypeName("freeSurfaceHeight");


    // Constructors

    //- Construct from patch and internal field
    freeSurfaceHeightFvPatchScalarField
    (
    	const fvPatch&,
    	const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    freeSurfaceHeightFvPatchScalarField
    (
    	const fvPatch&,
    	const DimensionedField<scalar, volMesh>&,
    	const dictionary&
    );

    //- Construct by mapping given freeSurfaceHeightFvPatchScalarField onto a new patch
    freeSurfaceHeightFvPatchScalarField
    (
    	const freeSurfaceHeightFvPatchScalarField&,
    	const fvPatch&,
    	const DimensionedField<scalar, volMesh>&,
    	const fvPatchFieldMapper&
    );

    //- Construct as copy
    freeSurfaceHeightFvPatchScalarField
    (
    	const freeSurfaceHeightFvPatchScalarField&
    );

    //- Construct and return a clone

    virtual tmp<fvPatchField<scalar> > clone() const 
	{
        return tmp<fvPatchField<scalar> >
        (
        	new freeSurfaceHeightFvPatchScalarField(*this)
        );
    }

    //- Construct as copy setting internal field reference
    freeSurfaceHeightFvPatchScalarField
    (
    	const freeSurfaceHeightFvPatchScalarField&,
    	const DimensionedField<scalar, volMesh>&
    );

    //- Construct and return a clone setting internal field reference

    virtual tmp<fvPatchField<scalar> > clone
    (
    	const DimensionedField<scalar, volMesh>& iF
    ) const 
	{
        return tmp<fvPatchField<scalar> >
        (
            new freeSurfaceHeightFvPatchScalarField(*this, iF)
        );
    }


    // Member functions

        // Mapping functions


    // Return defining fields


    //- Update the coefficients associated with the patch field
	    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
