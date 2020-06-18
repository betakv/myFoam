/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::profileModelUp

Description
    Base class for profile models

SourceFiles
    profileModelUp.C

\*---------------------------------------------------------------------------*/

#ifndef profileModelUp_H
#define profileModelUp_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class profileModelUp Declaration
\*---------------------------------------------------------------------------*/

class profileModelUp
{
protected:

    // Protected data

        //- Coefficients dictionary
        const dictionary dict_;

        //- Name of profile model
        const word name_;

        //- File name (optional)
        fileName fName_;


    // Protected Member Functions

        //- Return ture if file name is set
        bool readFromFile() const;


public:

    //- Runtime type information
    TypeName("profileModelUp");


        // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            profileModelUp,
            dictionary,
            (
                const dictionary& dict,
                const word& modelName
            ),
            (dict, modelName)
        );


    // Selectors

        //- Return a reference to the selected fvOption model
        static autoPtr<profileModelUp> New(const dictionary& dict);


    //- Constructor
    profileModelUp(const dictionary& dict, const word& modelName);


    //- Destructor
    virtual ~profileModelUp() = default;


    // Member functions

        // Access

            //- Return const access to the source name
            const word& name() const;


        // Evaluation

            //- Return the Cd and Cl for a given angle-of-attack
            virtual void Cdl
            (
                const scalar Re,
                const scalar alpha,
                scalar& Cd,
                scalar& Cl
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
