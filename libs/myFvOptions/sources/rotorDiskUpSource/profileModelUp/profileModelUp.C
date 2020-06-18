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

#include "profileModelUp.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(profileModelUp, 0);
    defineRunTimeSelectionTable(profileModelUp, dictionary);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::profileModelUp::readFromFile() const
{
    return fName_ != fileName::null;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::profileModelUp::profileModelUp(const dictionary& dict, const word& name)
:
    dict_(dict),
    name_(name),
    fName_(dict.lookupOrDefault<fileName>("file", fileName::null))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::profileModelUp::name() const
{
    return name_;
}


Foam::autoPtr<Foam::profileModelUp> Foam::profileModelUp::New
(
    const dictionary& dict
)
{
    const word& modelName = dict.dictName();

    const word modelType(dict.get<word>("type"));

    Info<< "    - creating " << modelType << " profile " << modelName << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "profileModelUp",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<profileModelUp>(cstrIter()(dict, modelName));
}


// ************************************************************************* //
