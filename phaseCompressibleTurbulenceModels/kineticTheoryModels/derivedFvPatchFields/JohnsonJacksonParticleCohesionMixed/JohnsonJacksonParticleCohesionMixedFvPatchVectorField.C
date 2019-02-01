/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "JohnsonJacksonParticleCohesionMixedFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        JohnsonJacksonParticleCohesionMixedFvPatchVectorField
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::JohnsonJacksonParticleCohesionMixedFvPatchVectorField::
JohnsonJacksonParticleCohesionMixedFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    specularityCoefficient_("specularityCoefficient", dimless, 0)
{}


Foam::JohnsonJacksonParticleCohesionMixedFvPatchVectorField::
JohnsonJacksonParticleCohesionMixedFvPatchVectorField
(
    const JohnsonJacksonParticleCohesionMixedFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    specularityCoefficient_(ptf.specularityCoefficient_)
{}


Foam::JohnsonJacksonParticleCohesionMixedFvPatchVectorField::
JohnsonJacksonParticleCohesionMixedFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    specularityCoefficient_
    (
        "specularityCoefficient",
        dimless,
        dict.lookup("specularityCoefficient")
    )
{
    if
    (
        (specularityCoefficient_.value() < 0)
     || (specularityCoefficient_.value() > 1)
    )
    {
        FatalErrorInFunction
            << "The specularity coefficient has to be between 0 and 1"
            << abort(FatalError);
    }

    fvPatchVectorField::operator=
    (
        vectorField("value", dict, p.size())
    );
}


Foam::JohnsonJacksonParticleCohesionMixedFvPatchVectorField::
JohnsonJacksonParticleCohesionMixedFvPatchVectorField
(
    const JohnsonJacksonParticleCohesionMixedFvPatchVectorField& ptf
)
:
    mixedFvPatchVectorField(ptf),
    specularityCoefficient_(ptf.specularityCoefficient_)
{}


Foam::JohnsonJacksonParticleCohesionMixedFvPatchVectorField::
JohnsonJacksonParticleCohesionMixedFvPatchVectorField
(
    const JohnsonJacksonParticleCohesionMixedFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(ptf, iF),
    specularityCoefficient_(ptf.specularityCoefficient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::JohnsonJacksonParticleCohesionMixedFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchVectorField::autoMap(m);
}


void Foam::JohnsonJacksonParticleCohesionMixedFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    mixedFvPatchVectorField::rmap(ptf, addr);
}


void Foam::JohnsonJacksonParticleCohesionMixedFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // lookup the fluid model and the phase
    const twoPhaseSystem& fluid = db().lookupObject<twoPhaseSystem>
    (
        "phaseProperties"
    );

    const phaseModel& phased
    (
        fluid.phase1().name() == internalField().group()
      ? fluid.phase1()
      : fluid.phase2()
    );

    // lookup all the fields on this patch
    const fvPatchScalarField& alpha
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            phased.volScalarField::name()
        )
    );
    
    const fvPatchScalarField& rho
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            phased.volScalarField::name()
        )
    );
    /*
    const fvPatchScalarField& d
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            phased.volScalarField::name()
        )
    );
    */
    const fvPatchScalarField& gs0
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("gs0", phased.name())
        )
    );

    const scalarField nu
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("nut", phased.name())
        )
    );

    const scalarField nuFric
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("nuFric", phased.name())
        )
    );

    
    const scalarField nuCohesion
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("nuCohesion", phased.name())
        )
    );

    /*
    const scalarField pcohesion
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("pcohesion", phased.name())
        )
    );
    */
   /*
    const fvPatchScalarField& uTerminal
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("uTerminal", phased.name())
        )
    );
    */
    dimensionedScalar cohesionFactor
    (
        "cohesionFactor",
        dimless,
        db()
       .lookupObject<IOdictionary>
        (
            IOobject::groupName("turbulenceProperties", phased.name())
        )
       .subDict("RAS")
       .subDict("kineticTheoryCoeffs")
       .lookup("cohesionFactor")
    ); 
    
    dimensionedScalar cohesionForce
    (
        "cohesionForce",
        dimless,
        db()
       .lookupObject<IOdictionary>
        (
            IOobject::groupName("turbulenceProperties", phased.name())
        )
       .subDict("RAS")
       .subDict("kineticTheoryCoeffs")
       .lookup("cohesionForce")
    ); 
    
    
    
    
    
    
    
    word ThetaName(IOobject::groupName("Theta", phased.name()));

    const fvPatchScalarField& Theta
    (
        db().foundObject<volScalarField>(ThetaName)
      ? patch().lookupPatchField<volScalarField, scalar>(ThetaName)
      : alpha
    );

    // lookup the packed volume fraction
    dimensionedScalar alphaMax
    (
        "alphaMax",
        dimless,
        db()
       .lookupObject<IOdictionary>
        (
            IOobject::groupName("turbulenceProperties", phased.name())
        )
       .subDict("RAS")
       .subDict("kineticTheoryCoeffs")
       .lookup("alphaMax")
    );

    // calculate the slip value fraction
    scalarField c
    (
        constant::mathematical::pi
       *alpha
       *gs0
       *specularityCoefficient_.value()
       *sqrt(3.0*Theta)
       /max(6.0*(nu - nuFric)*alphaMax.value(), SMALL)
    );

    tmp<vectorField> n = patch().nf();
     // Get the tangential component from the internalField (zero-gradient)
    vectorField Ut(patchInternalField());
    Ut = n()*(Ut & n());

    //Info << "ut" << Ut;
    
    

    this->refValue() = (-1.0*Ut*nuCohesion)/max(mag(Ut)*c*rho*nu, SMALL);
    
    //-1.0*n*cohesionFactor.value()*6.0*sqrt(2.0)*cohesionForce.value()*alphaMax.value()//*mag(alpha.snGrad())
           ///max( uTerminal*d*rho*specularityCoefficient_.value()*alpha*gs0*sqrt(3.0)*(scalar(1) - alpha), SMALL );

    this->refGrad() = vector::zero;
    
    this->valueFraction() = c/(c + patch().deltaCoeffs());

    mixedFvPatchVectorField::updateCoeffs();
}


void Foam::JohnsonJacksonParticleCohesionMixedFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("specularityCoefficient")
        << specularityCoefficient_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// ************************************************************************* //
