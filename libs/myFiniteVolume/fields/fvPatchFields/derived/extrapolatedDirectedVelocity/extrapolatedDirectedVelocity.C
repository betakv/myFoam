/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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

#include "extrapolatedDirectedVelocity.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "totalPressureFvPatchScalarField.H"
#include "totalTemperatureFvPatchScalarField.H"
#include <stdlib.h>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extrapolatedDirectedVelocityFvPatchVectorField::
extrapolatedDirectedVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    zeroGradientFvPatchVectorField(p, iF),
    inletDir_(p.size())
{}


Foam::extrapolatedDirectedVelocityFvPatchVectorField::
extrapolatedDirectedVelocityFvPatchVectorField
(
    const extrapolatedDirectedVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchVectorField(ptf, p, iF, mapper),
    inletDir_(ptf.inletDir_, mapper)
{}


Foam::extrapolatedDirectedVelocityFvPatchVectorField::
extrapolatedDirectedVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchVectorField(p, iF),
    inletDir_("inletDirection", dict, p.size())
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}

Foam::extrapolatedDirectedVelocityFvPatchVectorField::
extrapolatedDirectedVelocityFvPatchVectorField
(
    const extrapolatedDirectedVelocityFvPatchVectorField& sfspvf
)
:
    zeroGradientFvPatchVectorField(sfspvf),
    inletDir_(sfspvf.inletDir_)
{}


Foam::extrapolatedDirectedVelocityFvPatchVectorField::
extrapolatedDirectedVelocityFvPatchVectorField
(
    const extrapolatedDirectedVelocityFvPatchVectorField& sfspvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    zeroGradientFvPatchVectorField(sfspvf, iF),
    inletDir_(sfspvf.inletDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::extrapolatedDirectedVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    zeroGradientFvPatchVectorField::autoMap(m);
    inletDir_.autoMap(m);
}


void Foam::extrapolatedDirectedVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    zeroGradientFvPatchVectorField::rmap(ptf, addr);

    const extrapolatedDirectedVelocityFvPatchVectorField& tiptf =
        refCast<const extrapolatedDirectedVelocityFvPatchVectorField>(ptf);

    inletDir_.rmap(tiptf.inletDir_, addr);
}

void Foam::extrapolatedDirectedVelocityFvPatchVectorField::evaluate
(

 )
{
  if (!updated())
    {
      updateCoeffs();
    }
  
  zeroGradientFvPatchVectorField::evaluate();

  const vectorField n(patch().nf());
  operator==( inletDir() * (n & *this) / (n & inletDir()) );    

}


void Foam::extrapolatedDirectedVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    inletDir_.writeEntry("inletDirection", os);
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        extrapolatedDirectedVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
