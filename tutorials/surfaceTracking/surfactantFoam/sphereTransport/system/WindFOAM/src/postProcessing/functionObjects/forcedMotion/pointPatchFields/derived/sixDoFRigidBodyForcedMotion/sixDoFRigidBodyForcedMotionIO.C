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

#include "sixDoFRigidBodyForcedMotion.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyForcedMotion::write(Ostream& os) const
{
    motionState_.write(os);

    os.writeKeyword("initialCentreOfMass")
        << initialCentreOfMass_ << token::END_STATEMENT << nl;
    os.writeKeyword("initialOrientation")
        << initialQ_ << token::END_STATEMENT << nl;
    os.writeKeyword("translationAmplitude")
        << translationAmplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("translationFrequency")
        << translationFrequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("rotationAmplitude")
        << rotationAmplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("rotationFrequency")
        << rotationFrequency_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, sixDoFRigidBodyForcedMotion& sDoFRBFM)
{
    is  >> sDoFRBFM.motionState_
        >> sDoFRBFM.initialCentreOfMass_
        >> sDoFRBFM.initialQ_
        >> sDoFRBFM.translationAmplitude_
        >> sDoFRBFM.translationFrequency_
        >> sDoFRBFM.rotationAmplitude_
        >> sDoFRBFM.rotationFrequency_;

    // Check state of Istream
    is.check
    (
        "Foam::Istream& Foam::operator>>"
        "(Foam::Istream&, Foam::sixDoFRigidBodyForcedMotion&)"
    );

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const sixDoFRigidBodyForcedMotion& sDoFRBFM
)
{
    os  << sDoFRBFM.motionState()
        << token::SPACE << sDoFRBFM.initialCentreOfMass()
        << token::SPACE << sDoFRBFM.initialQ()
        << token::SPACE << sDoFRBFM.translationAmplitude()
        << token::SPACE << sDoFRBFM.translationFrequency()
        << token::SPACE << sDoFRBFM.rotationAmplitude()
        << token::SPACE << sDoFRBFM.rotationFrequency();

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<(Foam::Ostream&, "
        "const Foam::sixDoFRigidBodyForcedMotion&)"
    );

    return os;
}


// ************************************************************************* //
