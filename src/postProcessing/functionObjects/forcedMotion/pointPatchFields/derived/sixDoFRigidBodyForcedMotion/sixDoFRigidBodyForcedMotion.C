/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICLUAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "sixDoFRigidBodyForcedMotion.H"
#include "mathematicalConstants.H"
#include "fvMesh.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyForcedMotion::sixDoFRigidBodyForcedMotion()
:
    motionState_(),
    initialCentreOfMass_(vector::zero),
    initialQ_(I),
    translationAmplitude_(vector::zero),
    translationFrequency_(vector::zero),
    rotationAmplitude_(0),
    rotationFrequency_(0)
{}


Foam::sixDoFRigidBodyForcedMotion::sixDoFRigidBodyForcedMotion
(
    const point& centreOfMass,
    const tensor& Q,
    const point& initialCentreOfMass,
    const tensor& initialQ,
    const vector translationAmplitude,
    const vector translationFrequency,
    const scalar rotationAmplitude,
    const scalar rotationFrequency
)
:
    motionState_
    (
        centreOfMass,
        Q
    ),
    initialCentreOfMass_(initialCentreOfMass),
    initialQ_(initialQ),
    translationAmplitude_(translationAmplitude),
    translationFrequency_(translationFrequency),
    rotationAmplitude_(rotationAmplitude),
    rotationFrequency_(rotationFrequency)
{}


Foam::sixDoFRigidBodyForcedMotion::sixDoFRigidBodyForcedMotion
(
    const dictionary& dict
)
:
    motionState_(dict),
    initialCentreOfMass_
    (
        dict.lookupOrDefault("initialCentreOfMass", centreOfMass())
    ),
    initialQ_
    (
        dict.lookupOrDefault("initialOrientation", Q())
    ),
    translationAmplitude_
    (
        dict.lookupOrDefault("translationAmplitude", translationAmplitude())
    ),
    translationFrequency_
    (
        dict.lookupOrDefault("translationFrequency", translationFrequency())
    ),
    rotationAmplitude_
    (
        dict.lookupOrDefault("rotationAmplitude", rotationAmplitude())
    ),
    rotationFrequency_
    (
        dict.lookupOrDefault("rotationFrequency", rotationFrequency())
    )
{}


Foam::sixDoFRigidBodyForcedMotion::sixDoFRigidBodyForcedMotion
(
    const sixDoFRigidBodyForcedMotion& sDoFRBFM
)
:
    motionState_(sDoFRBFM.motionState()),
    initialCentreOfMass_(sDoFRBFM.initialCentreOfMass()),
    initialQ_(sDoFRBFM.initialQ()),
    translationAmplitude_(sDoFRBFM.translationAmplitude()),
    translationFrequency_(sDoFRBFM.translationFrequency()),
    rotationAmplitude_(sDoFRBFM.rotationAmplitude()),
    rotationFrequency_(sDoFRBFM.rotationFrequency())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyForcedMotion::~sixDoFRigidBodyForcedMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyForcedMotion::updatePosition
(
    scalar t,
    scalar deltaT
)
{
    // Daniel Wei, 12/10/2010
    // I changed the updatePosition function, so that
    // in forced vibration, the position is the function of time,
    // not just deltaT, this is the difference between 
    // sixDoFRigidBodyForcedMotion and sixDoFRigidBodyMotion

    if (Pstream::master())
    {
        // Daniel Wei, 12/10/2010
        // t        :   the current time;
        // deltaT   :   the time increasement;
        // These two variables are got when the sixDoFBodyForcedDisplacement
        //   function is called;
        
        //if (t == 0) { t = 50;}
        scalar oldTime = t - deltaT;    
        scalar curTime = t;
        

        scalar pi = Foam::mathematicalConstant::pi;

        /*// An input template
        scalar rotationAmplitude_(0.707);
        scalar rotationFrequency_(0.0);
        vector translationFrequency_(0.0,0.25,0.0);
        vector translationAmplitude_(0.0,5.0,0.0);
        */

        // For heaving motion
        vector translationVector(vector::zero);
        translationVector.x() = translationAmplitude_.x()
            *(
                Foam::sin( 2 * pi * translationFrequency_.x() * curTime )
               -Foam::sin( 2 * pi * translationFrequency_.x() * oldTime )  
            );
        translationVector.y() = translationAmplitude_.y()
            *(
                Foam::sin( 2 * pi * translationFrequency_.y() * curTime )
               -Foam::sin( 2 * pi * translationFrequency_.y() * oldTime ) 
            );
        translationVector.z() = 0.0;
        
        // For rotating motion
        scalar alphaOld = rotationAmplitude_
            *Foam::sin( 2 * pi * rotationFrequency_ * oldTime );
        scalar alphaCur = rotationAmplitude_
            *Foam::sin( 2 * pi * rotationFrequency_ * curTime );
        tensor RzOld
            (
                Foam::cos(alphaOld), Foam::sin(alphaOld), 0,
                -Foam::sin(alphaOld), Foam::cos(alphaOld), 0,
                0, 0, 1
            );
        tensor RzCur
            (
                Foam::cos(alphaCur), Foam::sin(alphaCur), 0,
                -Foam::sin(alphaCur), Foam::cos(alphaCur), 0,
                0, 0, 1
            );
        tensor rotationTensor = RzCur - RzOld;        
        
        centreOfMass() += translationVector;
        
        //Q() += rotationTensor;
        Q() = RzCur;
        
        // For debug's sake!
        Info<< "translationVector: " << translationVector << nl
            << "rotationTensor   : " << rotationTensor << nl
            << "centreOfMass     : " << centreOfMass() << nl
            << "Orientation      : " << Q() << nl
            << endl;
    }

    Pstream::scatter(motionState_);
}

// ************************************************************************* //
