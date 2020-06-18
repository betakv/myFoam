/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "EAM.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> EAM<BasicTurbulenceModel>::k
(
    const tmp<volTensorField>& gradU
) const
{
    volSymmTensorField D(symm(gradU));

    volScalarField a(this->Ce_/this->delta());
    volScalarField b((2.0/3.0)*tr(D));
    volScalarField c(2*Ck_*this->delta()*(dev(D) && D));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("k", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
            sqr((-b + sqrt(sqr(b) + 4*a*c))/(2*a))
        )
    );
}

template<class BasicTurbulenceModel>
void EAM<BasicTurbulenceModel>::correctNut()
{
  correctNonlinearStress(fvc::grad(this->U_));
}

template<class BasicTurbulenceModel>
void EAM<BasicTurbulenceModel>::correctNonlinearStress(const volTensorField& gradU)
{
    const volSymmTensorField S("S", dev(symm(gradU)));
          volScalarField magS("magS", mag(S));
          magS.max(SMALL);
    
    const volTensorField     W("W", skew(gradU));
    
    const volScalarField LL
    (
           simpleFilter_(this->U_ & this->U_)
         - (simpleFilter_(this->U_) & simpleFilter_(this->U_))
    );
    
    volScalarField MM
    (
        filter_(this->delta())*filter_(this->delta()) * simpleFilter_(magS*magS)
      - simpleFilter_(this->delta()*this->delta()*magS*magS)
    );
    
    MM.max(SMALL);
    
    volScalarField c 
    (
        0.5 * simpleFilter_((LL*MM)/(MM*MM))
    );
    
    c.max(SMALL);
    
    const volScalarField c1
    (
        C1_*sqrt(C3_)/pow(2*Cs_,2.5) *pow(c,1.25)
    );
    
    volScalarField tau
    (
        C3_*1.5*Ck_ / (2*Cs_) * 1/magS *sqrt(c)
    );
    
    tau.max(SMALL);

    const volScalarField Ksgs
    (
        c*this->delta()*this->delta()*magSqr(S)
    );
    
    const volSymmTensorField Sstar(S*tau);
    const volTensorField Wstar(W*tau);
    
    const volScalarField beta1
    (
        -33/20*9/4*c1
        /(sqr(9/4*c1) + magSqr(Wstar))
    );
    
    const volScalarField beta4
    (
        -33/20*1
        /(sqr(9/4*c1) + magSqr(Wstar))
    );
    this->nut_ = tau*beta1*Ksgs;
    this->nut_.max(SMALL);
    this->nut_.correctBoundaryConditions();

    this->nonlinearStress_ = beta4*Ksgs*symm((Sstar&Wstar)-(Wstar&Sstar));
    this->nonlinearStress_.correctBoundaryConditions();
    
    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
EAM<BasicTurbulenceModel>::EAM
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    LESnonlinearEddyViscosity<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ck",
            this->coeffDict_,
            1.5
        )
    ),
    
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            3.12
        )
    ),
    
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0.91
        )
    ),
    
    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.1
        )
    ),
    
    simpleFilter_(this->mesh_),
    filterPtr_(LESfilter::New(this->mesh_, this->coeffDict())),
    filter_(filterPtr_())
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool EAM<BasicTurbulenceModel>::read()
{
    if (LESnonlinearEddyViscosity<BasicTurbulenceModel>::read())
    {
        Ck_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        Cs_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> EAM<BasicTurbulenceModel>::epsilon() const
{
    volScalarField k(this->k(fvc::grad(this->U_)));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
            this->Ce_*k*sqrt(k)/this->delta()
        )
    );
}


template<class BasicTurbulenceModel>
void EAM<BasicTurbulenceModel>::correct()
{
    LESnonlinearEddyViscosity<BasicTurbulenceModel>::correct();
    correctNonlinearStress(fvc::grad(this->U_));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
