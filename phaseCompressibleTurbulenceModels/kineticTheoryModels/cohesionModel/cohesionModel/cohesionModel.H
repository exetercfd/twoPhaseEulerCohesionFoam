/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::kineticTheoryModels::cohesionModel

SourceFiles
    cohesionModel.C

\*---------------------------------------------------------------------------*/

#ifndef cohesionModel_H
#define cohesionModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                           Class cohesionModel Declaration
\*---------------------------------------------------------------------------*/

class cohesionModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        cohesionModel(const cohesionModel&);

        //- Disallow default bitwise assignment
        void operator=(const cohesionModel&);


protected:

    // Protected data

        //- Reference to higher-level dictionary for re-read
        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("cohesionModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        cohesionModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        cohesionModel(const dictionary& dict);


    // Selectors

        static autoPtr<cohesionModel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~cohesionModel();


    // Member Functions

        virtual tmp<volScalarField> frictionalPressure
        (
            const phaseModel& phase,
            const volScalarField& Theta,
            const volScalarField& da,
            const volScalarField& uTerminal
        ) const = 0;
/*
        virtual tmp<volScalarField> frictionalPressurePrime
        (
            const phaseModel& phase,
            const volScalarField& Theta,
            const volScalarField& da
        ) const = 0;
*/
        virtual tmp<volScalarField> nu
        (
            const phaseModel& phase,
            const volScalarField& pf,
            const volSymmTensorField& D,
            const volScalarField& Theta,
            const volScalarField& da
        ) const = 0;

        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
