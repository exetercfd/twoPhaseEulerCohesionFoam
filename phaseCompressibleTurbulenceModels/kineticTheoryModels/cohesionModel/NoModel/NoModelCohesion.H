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
    Foam::kineticTheoryModels::cohesionModels::NoModel

Description

SourceFiles
    NoneCohesion.C

\*---------------------------------------------------------------------------*/

#ifndef NoModel_H
#define NoModel_H

#include "cohesionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace cohesionModels
{

/*---------------------------------------------------------------------------*\
               Class NoModel Declaration
\*---------------------------------------------------------------------------*/

class NoModel
:
    public cohesionModel
{
    // Private data

        dictionary coeffDict_;

        //- Material constant for frictional normal stress
        dimensionedScalar cohesionFactor_;

        //- Material constant for frictional normal stress
        dimensionedScalar cohesionForce_;

        //- Material constant for frictional normal stress
        //dimensionedScalar p_;

        //- Angle of internal friction
        //dimensionedScalar phi_;

        //- Lower limit for (alphaMax - alpha1)
        //dimensionedScalar alphaDeltaMin_;


public:

    //- Runtime type information
    TypeName("NoModel");


    // Constructors

        //- Construct from components
        NoModel(const dictionary& dict);


    //- Destructor
    virtual ~NoModel();


    // Member functions

        virtual tmp<volScalarField> frictionalPressure
        (
            const phaseModel& phase,
            const volScalarField& Theta,
            const volScalarField& da,
            const volScalarField& uTerminal
        ) const;
/*
        virtual tmp<volScalarField> frictionalPressurePrime
        (
            const phaseModel& phase,
            const volScalarField& Theta,
            const volScalarField& da
        ) const;
*/
        virtual tmp<volScalarField> nu
        (
            const phaseModel& phase,
            const volScalarField& pf,
            const volSymmTensorField& D,
            const volScalarField& Theta,
            const volScalarField& da
         
        ) const;

        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cohesionModels
} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
