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

#include "MakkawiSrivastavaCohesion.H"
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
    defineTypeNameAndDebug(MakkawiSrivastava, 0);

    addToRunTimeSelectionTable
    (
        cohesionModel,
        MakkawiSrivastava,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::cohesionModels::MakkawiSrivastava::
MakkawiSrivastava
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

Foam::kineticTheoryModels::cohesionModels::MakkawiSrivastava::
~MakkawiSrivastava()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::cohesionModels::MakkawiSrivastava::
frictionalPressure
(
    const phaseModel& phase,
    const volScalarField& Theta,
    const volScalarField& da,
    const volScalarField& uTerminal
) const
{
    const volScalarField& alpha = phase;

    return (cohesionFactor_*cohesionForce_*6.0*sqrt(2.0)*sqrt(Theta)*mag(fvc::grad(alpha)))/(da*uTerminal + SMALL*dimensionedScalar("1.0", dimensionSet(0, 2, -1, 0, 0), 1.0));
    
}

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::cohesionModels::MakkawiSrivastava::nu
(
    const phaseModel& phase,
    const volScalarField& pf,
    const volSymmTensorField& D,
    const volScalarField& Theta,
    const volScalarField& da
) const
{
    const volScalarField& alpha = phase;
    
    tmp<volScalarField> tnu
    (
        new volScalarField
        (
            IOobject
            (
                "MakkawiSrivastava:nu",
                phase.mesh().time().timeName(),
                phase.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase.mesh(),
            dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0), 0.0)
        )
    );

    volScalarField& nuf = tnu.ref();

    forAll(D, celli)
    {
        //if (alpha[celli] > alphaMinFriction.value())
        //{
            nuf[celli] =
                pf[celli]*constant::mathematical::pi
                /(
                    6.0*(1.0-alpha[celli])*sqrt( (D[celli] && D[celli]) ) //+ (Theta[celli]/sqr(da[celli]))  )
                    + SMALL
                );
        //}
    }

    const fvPatchList& patches = phase.mesh().boundary();
    //const volVectorField& U = phase.U();

    volScalarField::Boundary& nufBf = nuf.boundaryFieldRef();

    //if (debuggingKineticTheory)
    //{
        //Info << "wallFriction = " << wallFriction << endl;
        Info<< "max(nufBf) before: " << gMax(nufBf) << endl;
        Info<< "min(nufBf) before: " << gMin(nufBf) << endl;
    //}
    
    forAll(patches, patchi)
    {
        if (!patches[patchi].coupled())
        {

                nufBf[patchi] =
                    pf.boundaryField()[patchi]*constant::mathematical::pi
                    /(
                        6.0*(1.0-alpha.boundaryField()[patchi])*sqrt( (D.boundaryField()[patchi] && D.boundaryField()[patchi]) ) //+
                        //(Theta.boundaryField()[patchi]/sqr(da.boundaryField()[patchi]))  )
                        + SMALL
                    );
        }
    }
    
    //if (debuggingKineticTheory)
    //{
        Info<< "max(nufBf) after: " << gMax(nufBf) << endl;
        Info<< "min(nufBf) after: " << gMin(nufBf) << endl;
    //}
    
    // Correct coupled BCs
    nuf.correctBoundaryConditions();

    return tnu;
    
}


bool Foam::kineticTheoryModels::cohesionModels::MakkawiSrivastava::read()
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
