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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline const Foam::sixDoFRigidBodyForcedMotionState&
Foam::sixDoFRigidBodyForcedMotion::motionState() const
{
    return motionState_;
}


inline const Foam::point&
Foam::sixDoFRigidBodyForcedMotion::initialCentreOfMass() const
{
    return initialCentreOfMass_;
}


inline const Foam::tensor&
Foam::sixDoFRigidBodyForcedMotion::initialQ() const
{
    return initialQ_;
}


inline const Foam::tensor& Foam::sixDoFRigidBodyForcedMotion::Q() const
{
    return motionState_.Q();
}


inline const Foam::vector&
Foam::sixDoFRigidBodyForcedMotion::translationAmplitude() const
{
    return translationAmplitude_;
}


inline const Foam::vector&
Foam::sixDoFRigidBodyForcedMotion::translationFrequency() const
{
    return translationFrequency_;
}


inline const Foam::scalar&
Foam::sixDoFRigidBodyForcedMotion::rotationAmplitude() const
{
    return rotationAmplitude_;
}


inline const Foam::scalar&
Foam::sixDoFRigidBodyForcedMotion::rotationFrequency() const
{
    return rotationFrequency_;
}


inline Foam::point& Foam::sixDoFRigidBodyForcedMotion::initialCentreOfMass()
{
    return initialCentreOfMass_;
}


inline Foam::tensor& Foam::sixDoFRigidBodyForcedMotion::initialQ()
{
    return initialQ_;
}


inline Foam::tensor& Foam::sixDoFRigidBodyForcedMotion::Q()
{
    return motionState_.Q();
}


inline Foam::vector&
Foam::sixDoFRigidBodyForcedMotion::translationAmplitude()
{
    return translationAmplitude_;
}


inline Foam::vector&
Foam::sixDoFRigidBodyForcedMotion::translationFrequency()
{
    return translationFrequency_;
}


inline Foam::scalar&
Foam::sixDoFRigidBodyForcedMotion::rotationAmplitude()
{
    return rotationAmplitude_;
}


inline Foam::scalar&
Foam::sixDoFRigidBodyForcedMotion::rotationFrequency()
{
    return rotationFrequency_;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::tmp<Foam::pointField>
Foam::sixDoFRigidBodyForcedMotion::currentPosition(const pointField& pInitial) const
{
    return
    (
        centreOfMass()
      + (Q() & (pInitial - initialCentreOfMass_)))
//      + (Q() & initialQ_.T() & (pInitial - initialCentreOfMass_)))
    ;
}


inline Foam::point Foam::sixDoFRigidBodyForcedMotion::currentPosition
(
    const point& pInitial
) const
{
    return
    (
        centreOfMass()
      + (Q() & (pInitial - initialCentreOfMass_))
//      + (Q() & initialQ_.T() & (pInitial - initialCentreOfMass_))
    );
}


inline const Foam::tensor&
Foam::sixDoFRigidBodyForcedMotion::orientation() const
{
    return Q();
}


inline Foam::tensor&
Foam::sixDoFRigidBodyForcedMotion::orientation()
{
    return Q();
}


inline const Foam::point& Foam::sixDoFRigidBodyForcedMotion::centreOfMass() const
{
    return motionState_.centreOfMass();
}


inline Foam::point& Foam::sixDoFRigidBodyForcedMotion::centreOfMass()
{
    return motionState_.centreOfMass();
}

// ************************************************************************* //
