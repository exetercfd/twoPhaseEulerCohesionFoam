/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "NoModelCohesion.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace cohesionModels
{
    defineTypeNameAndDebug(NoModel, 0);

    addToRunTimeSelectionTable
    (
        cohesionModel,
        NoModel,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::cohesionModels::NoModel::
NoModel
(
    const dictionary& dict
)
:
    cohesionModel(dict),
    coeffDict_(dict.optionalSubDict(typeName + "Coeffs")),
    cohesionFactor_("cohesionFactor", dimless, coeffDict_),
    cohesionForce_("cohesionForce", dimensionSet(1, 1, -2, 0, 0), coeffDict_)
    
    /*Fr_("Fr", dimensionSet(1, -1, -2, 0, 0), coeffDict_),
    eta_("eta", dimless, coeffDict_),
    p_("p", dimless, coeffDict_),
    phi_("phi", dimless, coeffDict_),
    alphaDeltaMin_("alphaDeltaMin", dimless, coeffDict_)
    */
{
    //phi_ *= constant::mathematical::pi/180.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::cohesionModels::NoModel::
~NoModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::cohesionModels::NoModel::
frictionalPressure
(
    const phaseModel& phase,
    const volScalarField& Theta,
    const volScalarField& da,
    const volScalarField& uTerminal
) const
{
    const volScalarField& alpha = phase;

    return 0.0*(cohesionFactor_*cohesionForce_*6.0*sqrt(2.0)*sqrt(Theta)*mag(fvc::grad(alpha)))/(da*uTerminal + SMALL*dimensionedScalar("1.0", dimensionSet(0, 2, -1, 0, 0), 1.0));
    
}

/*
Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::cohesionModels::NoModel::
frictionalPressurePrime
(
    const phaseModel& phase,
    const volScalarField& Theta,
    const volScalarField& da
) const
{
    const volScalarField& alpha = phase;

    return (6.0*sqrt(2.0)*sqrt(mag(Theta))*fvc::laplacian(alpha))/(mag(fvc::grad(alpha))*da + SMALL);
}
*/

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::cohesionModels::NoModel::nu
(
    const phaseModel& phase,
    const volScalarField& pf,
    const volSymmTensorField& D,
    const volScalarField& Theta,
    const volScalarField& da
) const
{
    const volScalarField& alpha = phase;
    
    return 0.0*dimensionedScalar("1.0", dimTime, 1.0)*pf*constant::mathematical::pi/(6.0*(1.0-alpha) + SMALL);
}


bool Foam::kineticTheoryModels::cohesionModels::NoModel::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");

    cohesionFactor_.read(coeffDict_);
    cohesionForce_.read(coeffDict_);
    /*p_.read(coeffDict_);

    phi_.read(coeffDict_);
    phi_ *= constant::mathematical::pi/180.0;

    alphaDeltaMin_.read(coeffDict_);
    */
    return true;
}


// ************************************************************************* //
