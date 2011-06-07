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
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "sixDoFRigidBodyForcedMotionState.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyForcedMotionState::sixDoFRigidBodyForcedMotionState()
:
    centreOfMass_(vector::zero),
    Q_(I)
{}


Foam::sixDoFRigidBodyForcedMotionState::sixDoFRigidBodyForcedMotionState
(
    const point& centreOfMass,
    const tensor& Q
)
:
    centreOfMass_(centreOfMass),
    Q_(Q)
{}


Foam::sixDoFRigidBodyForcedMotionState::sixDoFRigidBodyForcedMotionState
(
    const dictionary& dict
)
:
    centreOfMass_(dict.lookup("centreOfMass")),
    Q_(dict.lookupOrDefault("orientation", tensor(I)))
{}


Foam::sixDoFRigidBodyForcedMotionState::sixDoFRigidBodyForcedMotionState
(
    const sixDoFRigidBodyForcedMotionState& sDoFRBFMS
)
:
    centreOfMass_(sDoFRBFMS.centreOfMass()),
    Q_(sDoFRBFMS.Q())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyForcedMotionState::~sixDoFRigidBodyForcedMotionState()
{}


// ************************************************************************* //
