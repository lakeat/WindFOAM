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

Class
    sixDOFMSDqODE

Description
    6-DOF solver using quaternions

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "sixDOFMSDqODE.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Runtime type information
// Not possible because of I/O error: incorrect type, expecting dictionary
// HJ, 11/Feb/2008
// namespace Foam
// {
//     defineTypeNameAndDebug(sixDOFMSDqODE, 0);
// }


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sixDOFMSDqODE::setCoeffs()
{
    // Set ODE coefficients from position and rotation

    // Linear displacement relative to spring equilibrium
    const vector& Xval = Xrel_.value();
    coeffs_[0] = Xval.x();
    coeffs_[1] = Xval.y();
    coeffs_[2] = Xval.z();

    // Linear velocity
    const vector& Uval = U_.value();
    coeffs_[3] = Uval.x();
    coeffs_[4] = Uval.y();
    coeffs_[5] = Uval.z();

    // Rotational angle relative to spring equilibrium
    const vector& thetaVal = Trel_.value();
    coeffs_[6] = thetaVal.x();
    coeffs_[7] = thetaVal.y();
    coeffs_[8] = thetaVal.z();

    // Rotational velocity
    const vector& omegaVal = omega_.value();
    coeffs_[9] = omegaVal.x();
    coeffs_[10]= omegaVal.y();
    coeffs_[11]= omegaVal.z();
}


Foam::dimensionedVector Foam::sixDOFMSDqODE::A
(
    const dimensionedVector& xR,
    const dimensionedVector& uR
) const
{
    return
    (
       - (linSpringCoeffs_ & xR)    // translation spring
       - (linDampingCoeffs_ & uR)   // translation damping
       + force()
    )/mass_;
}


Foam::dimensionedVector Foam::sixDOFMSDqODE::omegaDot
(
    const dimensionedVector& thetaR,
    const dimensionedVector& omegaR
) const
{
    dimensionedScalar massTmp = mass_;
    massTmp.value() = 1;
    dimensionedVector tmpVector =
    (
       - (rotSpringCoeffs_ & thetaR)    // rotation spring
       - (rotDampingCoeffs_ & omegaR)   // rotation damping
       + moment()
    )/massTmp;

    tmpVector.value().x() = tmpVector.value().x()/momentOfInertia_.value().xx();
    tmpVector.value().y() = tmpVector.value().y()/momentOfInertia_.value().yy();
    tmpVector.value().z() = tmpVector.value().z()/momentOfInertia_.value().zz();

    return tmpVector;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::sixDOFMSDqODE::sixDOFMSDqODE(const IOobject& io)
:
    IOdictionary(io),
    mass_(lookup("mass")),
    momentOfInertia_(lookup("momentOfInertia")),
    Xequilibrium_(lookup("translationalEquilibriumPosition")),
    Tequilibrium_(lookup("rotationalEquilibriumPosition")),
    linSpringCoeffs_(lookup("translationalSpring")),
    linDampingCoeffs_(lookup("translationalDamping")),
    rotSpringCoeffs_(lookup("rotationalSpring")),
    rotDampingCoeffs_(lookup("rotationalDamping")),
    Xrel_(lookup("Xrel")),
    Trel_(lookup("Trel")),
    Xdiff_
    (
        "Xdiff",
        dimensionSet(0, 1, 0, 0, 0),
        vector(0,0,0)
    ),
    Tdiff_
    (
        "Tdiff",
        dimensionSet(0, 0, 0, 0, 0),
        vector(0,0,0)
    ),
    U_(lookup("U")),
    omega_(lookup("omega")),
    force_(lookup("force")),
    moment_(lookup("moment")),
    coeffs_(12, 0.0)
{
    Info<< "Constructing the sixDOFMSDbodies ...." << endl;
    setCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDOFMSDqODE::~sixDOFMSDqODE()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sixDOFMSDqODE::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    // Set the derivatives for translation displacement
    dydx[0] = y[3];
    dydx[1] = y[4];
    dydx[2] = y[5];

    dimensionedVector curX
	(
		"curX",
		dimLength,
		vector(y[0], y[1], y[2])
	);
    dimensionedVector curU
	(
		"curU",
		dimVelocity,
		vector(y[3], y[4], y[5])
	);

    const vector linAccel = A(curX, curU).value();

    dydx[3] = linAccel.x();
    dydx[4] = linAccel.y();
    dydx[5] = linAccel.z();

    // Set the derivatives for rotational displacement
	dydx[6] = y[9];
    dydx[7] = y[10];
    dydx[8] = y[11];

    dimensionedVector curTheta
    (
        "curTheta",
        dimensionSet(0, 0, 0, 0, 0),
        vector(y[6], y[7], y[8])
    );
    dimensionedVector curOmega
    (
        "curOmega",
        dimensionSet(0, 0, -1, 0, 0),
        vector(y[9], y[10], y[11])
    );

    const vector rotAccel = omegaDot(curTheta, curOmega).value();

    dydx[9] = rotAccel.x();
    dydx[10]= rotAccel.y();
    dydx[11]= rotAccel.z();
}


void Foam::sixDOFMSDqODE::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdx,
    scalarSquareMatrix& dfdy
) const
{
    Info << "jacobian(...)" << endl;
    notImplemented("sixDOFMSDqODE::jacobian(...) const");
}


void Foam::sixDOFMSDqODE::update(const scalar delta)
{
    // Update position
    vector Xold = Xrel_.value();

    Info<< "Xrel old = " << Xrel_ << endl;

    vector& Xval = Xrel_.value();

    Xval.x() = coeffs_[0];
    Xval.y() = coeffs_[1];
    Xval.z() = coeffs_[2];

    Info<< "Xrel new = " << Xrel_ << endl;

    Xdiff_.value() = Xval - Xold;

    vector& Uval = U_.value();

    Uval.x() = coeffs_[3];
    Uval.y() = coeffs_[4];
    Uval.z() = coeffs_[5];

    // Update omega
    vector Told = Trel_.value();

    Info<< "Trel old = " << Trel_ << endl;

    vector& Tval = Trel_.value();

    Tval.x() = coeffs_[6];
    Tval.y() = coeffs_[7];
    Tval.z() = coeffs_[8];

    Info<< "Trel new = " << Trel_ << endl;

    Tdiff_.value() = Tval - Told;

    vector& omegaVal = omega_.value();

    omegaVal.x() = coeffs_[9];
    omegaVal.y() = coeffs_[10];
    omegaVal.z() = coeffs_[11];
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

bool Foam::sixDOFMSDqODE::writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const sixDOFMSDqODE& sds)
{
    os.writeKeyword("mass") << tab << sds.mass_ << token::END_STATEMENT << nl;
    os.writeKeyword("momentOfInertia") << tab << sds.momentOfInertia_
        << token::END_STATEMENT << nl << nl;

    os.writeKeyword("translationalEquilibriumPosition") << tab << sds.Xequilibrium_
        << token::END_STATEMENT << nl;
    os.writeKeyword("rotationalEquilibriumPosition") << tab << sds.Tequilibrium_
        << token::END_STATEMENT << nl;
    os.writeKeyword("translationalSpring") << tab << sds.linSpringCoeffs_
        << token::END_STATEMENT << nl;
    os.writeKeyword("translationalDamping") << tab << sds.linDampingCoeffs_
        << token::END_STATEMENT << nl << nl;
    os.writeKeyword("rotationalSpring") << tab << sds.rotSpringCoeffs_
        << token::END_STATEMENT << nl;
    os.writeKeyword("rotationalDamping") << tab << sds.rotDampingCoeffs_
        << token::END_STATEMENT << nl << nl;

    os.writeKeyword("Xrel") << tab << sds.Xrel() << token::END_STATEMENT << nl;
    os.writeKeyword("Trel") << tab << sds.Trel() << token::END_STATEMENT << nl;
    os.writeKeyword("U") << tab << sds.U() << token::END_STATEMENT << nl;
    os.writeKeyword("omega") << tab << sds.omega()
        << token::END_STATEMENT << nl << nl;

    os.writeKeyword("force") << tab << sds.force()
        << token::END_STATEMENT << nl;
    os.writeKeyword("moment") << tab << sds.moment()
        << token::END_STATEMENT << endl;

    return os;
}


// ************************************************************************* //

