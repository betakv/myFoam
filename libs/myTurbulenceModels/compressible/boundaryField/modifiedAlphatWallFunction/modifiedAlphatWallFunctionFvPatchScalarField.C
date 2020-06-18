/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "modifiedAlphatWallFunctionFvPatchScalarField.H"
#include "compressibleTurbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

modifiedAlphatWallFunctionFvPatchScalarField::modifiedAlphatWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Prt_(0.85)
{}


modifiedAlphatWallFunctionFvPatchScalarField::modifiedAlphatWallFunctionFvPatchScalarField
(
    const modifiedAlphatWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Prt_(ptf.Prt_)
{}


modifiedAlphatWallFunctionFvPatchScalarField::modifiedAlphatWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Prt_(dict.lookupOrDefault<scalar>("Prt", 0.85))
{}


modifiedAlphatWallFunctionFvPatchScalarField::modifiedAlphatWallFunctionFvPatchScalarField
(
    const modifiedAlphatWallFunctionFvPatchScalarField& awfpsf
)
:
    fixedValueFvPatchScalarField(awfpsf),
    Prt_(awfpsf.Prt_)
{}


modifiedAlphatWallFunctionFvPatchScalarField::modifiedAlphatWallFunctionFvPatchScalarField
(
    const modifiedAlphatWallFunctionFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(awfpsf, iF),
    Prt_(awfpsf.Prt_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void modifiedAlphatWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    // Retrieve turbulence properties from model
    const compressibleTurbulenceModel& turbModel =
        db().lookupObject<compressibleTurbulenceModel>
        (
            IOobject::groupName
            (
                compressibleTurbulenceModel::propertiesName,
                internalField().group()
            )
        );
    const fluidThermo& thermoModel =
        db().lookupObject<fluidThermo>
        (
            basicThermo::dictName
        );    
    const scalarField& rhow = turbModel.rho().boundaryField()[patchi];
    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const tmp<scalarField> tnuw = turbModel.nu(patchi);

    const scalarField Pr = thermoModel.Cp()*thermoModel.mu()/thermoModel.kappa();
    
   scalarField viscosityRatioPr = tnutw/tnuw*Pr;
   viscosityRatioPr = max(viscosityRatioPr,SMALL);
    
    scalarField Prt = 0.5882+0.3254*viscosityRatioPr-0.09*sqr(viscosityRatioPr)*(1-exp(-3.6155 * 1/viscosityRatioPr ));
    
    Prt = 1.0/Prt;

    operator==(rhow*tnutw/Prt);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void modifiedAlphatWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("Prt", Prt_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    modifiedAlphatWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
