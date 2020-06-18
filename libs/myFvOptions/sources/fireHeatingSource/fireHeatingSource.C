/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "fireHeatingSource.H"
#include "fvMatrices.H"
#include "fvmLaplacian.H"
#include "fvcGrad.H"
#include "zeroGradientFvPatchField.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(fireHeatingSource, 0);

    addToRunTimeSelectionTable
    (
        cellSetOption,
        fireHeatingSource,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::fireHeatingSource::fireHeatingSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(sourceName, modelType, dict, mesh),
    TName_("T"),
    Q_(readScalar(coeffs_.lookup("Q"))),
    q_(coeffs_.lookup("q")),
    curTimeIndex_(-1)
{
    // Set the field name to that of the energy field from which the temperature
    // is obtained

    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    fieldNames_.setSize(1, thermo.he().name());

    applied_.setSize(fieldNames_.size(), false);

    read(dict);
    
    Info << endl;
    Info << "Fire heat source information" << endl;
    Info << "============================" << endl;
    Info << "  Zone name:    " << this->name() << endl;
    Info << "  Maximum fire heat power:              " << Q_ << endl;
    Info << "  Distribution function of heat source: " << q_ << endl;
    Info << "=======================================" << endl;
    Info << endl;
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::fireHeatingSource::~fireHeatingSource()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

void Foam::fv::fireHeatingSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    DebugInfo<< name() << ": applying source to " << eqn.psi().name() << endl;
    
    const scalarField& V = mesh_.V();
    scalar sumV = 0.0;
    
    forAll(cells_, i)
    {
        label celli = cells_[i];
        sumV += V[celli]*V[celli];
    }
    reduce(sumV, sumOp<scalar>());
    
    scalarField& heSource = eqn.source();
    
    scalar t = mesh_.time().timeIndex();
    dimensionedScalar dt = mesh_.time().deltaT();
        
    forAll(cells_, i)
    {
        label celli = cells_[i];
        
        scalar heatDistribution=0.0;
        
        for( int j = 0; j < q_.size(); j++ )
                heatDistribution += q_[j] * pow(t, j);
        heatDistribution=max(heatDistribution,0.0);         
        
        heSource[celli] += -1.0*V[celli]*V[celli]/sumV*Q_*heatDistribution;
//      heSource[celli] += -1.0*V[celli]*Q_*heatDistribution;


    }

}


bool Foam::fv::fireHeatingSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.readIfPresent("T", TName_);

        return true;
    }

    return false;
}

// ************************************************************************* //
