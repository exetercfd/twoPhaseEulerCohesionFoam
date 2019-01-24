/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
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

#include "JohnsonJacksonSchaefferSrivastavaFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{
    defineTypeNameAndDebug(JohnsonJacksonSchaefferSrivastava, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        JohnsonJacksonSchaefferSrivastava,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::
JohnsonJacksonSchaefferSrivastava::JohnsonJacksonSchaefferSrivastava
(
    const dictionary& dict
)
:
    frictionalStressModel(dict),
    coeffDict_(dict.optionalSubDict(typeName + "Coeffs")),
    Fr_("Fr", dimensionSet(1, -1, -2, 0, 0), coeffDict_),
    eta_("eta", dimless, coeffDict_),
    p_("p", dimless, coeffDict_),
    phi_("phi", dimless, coeffDict_),
    alphaDeltaMin_("alphaDeltaMin", dimless, coeffDict_)
{
    phi_ *= constant::mathematical::pi/180.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::
JohnsonJacksonSchaefferSrivastava::~JohnsonJacksonSchaefferSrivastava()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::
JohnsonJacksonSchaefferSrivastava::frictionalPressure
(
    const phaseModel& phase,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax
) const
{
    const volScalarField& alpha = phase;

    return
        Fr_*pow(max(alpha - alphaMinFriction, scalar(0)), eta_)
       /pow(max(alphaMax - alpha, alphaDeltaMin_), p_);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::
JohnsonJacksonSchaefferSrivastava::frictionalPressurePrime
(
    const phaseModel& phase,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax
) const
{
    const volScalarField& alpha = phase;

    return Fr_*
    (
        eta_*pow(max(alpha - alphaMinFriction, scalar(0)), eta_ - 1.0)
       *(alphaMax-alpha)
      + p_*pow(max(alpha - alphaMinFriction, scalar(0)), eta_)
    )/pow(max(alphaMax - alpha, alphaDeltaMin_), p_ + 1.0);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::
JohnsonJacksonSchaefferSrivastava::nu
(
    const phaseModel& phase,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const volScalarField& pf,
    const volSymmTensorField& D,
    const volScalarField& Theta,
    const volScalarField& da,
    const word& wallFriction,
    const Switch& debuggingKineticTheory
) const
{
    const volScalarField& alpha = phase;

    tmp<volScalarField> tnu
    (
        new volScalarField
        (
            IOobject
            (
                "JohnsonJacksonSchaefferSrivastava:nu",
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
        if (alpha[celli] > alphaMinFriction.value())
        {
            
            
            nuf[celli] =
                sqrt(2.0)*pf[celli]*sin(phi_.value())
                /(
                    sqrt( (D[celli] && D[celli]) + (Theta[celli]/sqr(da[celli]))  )
                    + SMALL
                );
            
            
            
            
            
            
            
            /*
             nuf[celli] =
                0.5*pf[celli]*sin(phi_.value())
               /(
                    sqrt(1.0/6.0*(sqr(D[celli].xx() - D[celli].yy())
                  + sqr(D[celli].yy() - D[celli].zz())
                  + sqr(D[celli].zz() - D[celli].xx()))
                  + sqr(D[celli].xy()) + sqr(D[celli].xz())
                  + sqr(D[celli].yz()) + (Theta[celli]/sqr(da[celli]))) + SMALL
                );
            */
            /*
            nuf[celli] =
                0.5*pf[celli]*sin(phi_.value())
               /(
                    sqrt((1.0/3.0)*sqr(tr(D[celli])) - invariantII(D[celli]) + (Theta[celli]/sqr(da[celli])))
                  + SMALL
                );
            */
        }
    }

    const fvPatchList& patches = phase.mesh().boundary();
    const volVectorField& U = phase.U();

    volScalarField::Boundary& nufBf = nuf.boundaryFieldRef();

    if (debuggingKineticTheory)
    {
        Info << "wallFriction = " << wallFriction << endl;
        Info<< "max(nufBf) before: " << gMax(nufBf) << endl;
        Info<< "min(nufBf) before: " << gMin(nufBf) << endl;
    }
    
    forAll(patches, patchi)
    {
        if (!patches[patchi].coupled())
        {
            if(wallFriction == "Coulomb")
            {
            
                nufBf[patchi] =
                    (
                        pf.boundaryField()[patchi]*sin(phi_.value())
                        /(
                            mag(U.boundaryField()[patchi].snGrad())
                            + SMALL
                        )
                    );
            }    
            else if(wallFriction == "JohnsonJacksonSchaefferSrivastava")
            {    
                
                
                nufBf[patchi] =
                    sqrt(2.0)*pf.boundaryField()[patchi]*sin(phi_.value())
                    /(
                        sqrt( (D.boundaryField()[patchi] && D.boundaryField()[patchi]) +
                        (Theta.boundaryField()[patchi]/sqr(da.boundaryField()[patchi]))  )
                        + SMALL
                    );
                
                
                
                /*
                nuf[patchi] =
                    0.5*pf.boundaryField()[patchi]*sin(phi_.value())
                    /(
                        sqrt(1.0/6.0*(sqr(D.boundaryField()[patchi].xx() - D.boundaryField()[patchi].yy())
                        + sqr(D.boundaryField()[patchi].yy() - D.boundaryField()[patchi].zz())
                        + sqr(D.boundaryField()[patchi].zz() - D.boundaryField()[patchi].xx()))
                        + sqr(D.boundaryField()[patchi].xy()) + sqr(D.boundaryField()[patchi].xz())
                        + sqr(D.boundaryField()[patchi].yz()) + (Theta.boundaryField()[patchi]/sqr(da.boundaryField()[patchi]))) + SMALL
                    );
                */
                
                /*
                nufBf[patchi] =
                    0.5*pf.boundaryField()[patchi]*sin(phi_.value())
                    /(
                        sqrt((1.0/3.0)*sqr(tr(D.boundaryField()[patchi])) - invariantII(D.boundaryField()[patchi]) + (Theta.boundaryField()[patchi]/sqr(da.boundaryField()[patchi])))
                        + SMALL
                    );
                    */
            }
            else
            {
                nufBf[patchi] =
                    (
                        pf.boundaryField()[patchi]*sin(phi_.value())
                        /(
                            mag(U.boundaryField()[patchi].snGrad())
                            + SMALL
                        )
                    );
            }
        }
    }
    
    if (debuggingKineticTheory)
    {
        Info<< "max(nufBf) after: " << gMax(nufBf) << endl;
        Info<< "min(nufBf) after: " << gMin(nufBf) << endl;
    }
    
    // Correct coupled BCs
    nuf.correctBoundaryConditions();

    return tnu;
}


bool Foam::kineticTheoryModels::frictionalStressModels::
JohnsonJacksonSchaefferSrivastava::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");

    Fr_.read(coeffDict_);
    eta_.read(coeffDict_);
    p_.read(coeffDict_);

    phi_.read(coeffDict_);
    phi_ *= constant::mathematical::pi/180.0;

    alphaDeltaMin_.read(coeffDict_);

    return true;
}


// ************************************************************************* //
