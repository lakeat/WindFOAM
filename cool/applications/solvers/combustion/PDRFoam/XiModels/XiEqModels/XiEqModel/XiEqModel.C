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

#include "XiEqModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(XiEqModel, 0);
    defineRunTimeSelectionTable(XiEqModel, dictionary);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiEqModel::XiEqModel
(
    const dictionary& XiEqProperties,
    const hhuCombustionThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su
)
:
    XiEqModelCoeffs_
    (
        XiEqProperties.subDict
        (
            word(XiEqProperties.lookup("XiEqModel")) + "Coeffs"
        )
    ),
    thermo_(thermo),
    turbulence_(turbulence),
    Su_(Su)
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::XiEqModel::~XiEqModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::XiEqModel::read(const dictionary& XiEqProperties)
{
    XiEqModelCoeffs_ = XiEqProperties.subDict(type() + "Coeffs");

    return true;
}


// ************************************************************************* //
