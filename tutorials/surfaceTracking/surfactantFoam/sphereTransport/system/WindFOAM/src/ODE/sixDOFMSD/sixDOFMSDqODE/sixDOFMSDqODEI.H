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

Description
    6-DOF solver using quaternions

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.
    Updates by Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "finiteRotation.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dimensionedScalar& Foam::sixDOFMSDqODE::mass() const
{
    return mass_;
}


const Foam::dimensionedDiagTensor& Foam::sixDOFMSDqODE::momentOfInertia() const
{
    return momentOfInertia_;
}


const Foam::dimensionedVector& Foam::sixDOFMSDqODE::Xequilibrium() const
{
    return Xequilibrium_;
}


const Foam::dimensionedVector& Foam::sixDOFMSDqODE::Tequilibrium() const
{
    return Tequilibrium_;
}


const Foam::dimensionedVector& Foam::sixDOFMSDqODE::Xrel() const
{
    return Xrel_;
}


const Foam::dimensionedVector& Foam::sixDOFMSDqODE::Trel() const
{
    return Trel_;
}


const Foam::dimensionedVector& Foam::sixDOFMSDqODE::U() const
{
    return U_;
}


const Foam::dimensionedVector& Foam::sixDOFMSDqODE::omega() const
{
    return omega_;
}


Foam::dimensionedVector Foam::sixDOFMSDqODE::A() const
{
    return A(Xrel_, U_);
}


Foam::dimensionedVector Foam::sixDOFMSDqODE::omegaDot() const
{
    return omegaDot(Trel_, omega_);
}


Foam::dimensionedVector Foam::sixDOFMSDqODE::X() const
{
    return Xequilibrium_ + Xrel_;
}


Foam::dimensionedVector Foam::sixDOFMSDqODE::T() const
{
    return Tequilibrium_ + Trel_;
}


Foam::dimensionedVector Foam::sixDOFMSDqODE::Xdiff() const
{
    return Xdiff_;
}


Foam::dimensionedVector Foam::sixDOFMSDqODE::Tdiff() const
{
    return Tdiff_;
}


Foam::dimensionedVector& Foam::sixDOFMSDqODE::U()
{
    return U_;
}


Foam::dimensionedVector& Foam::sixDOFMSDqODE::omega()
{
    return omega_;
}


const Foam::dimensionedVector& Foam::sixDOFMSDqODE::force() const
{
    return force_;
}


Foam::dimensionedVector& Foam::sixDOFMSDqODE::force()
{
    return force_;
}


const Foam::dimensionedVector& Foam::sixDOFMSDqODE::moment() const
{
    return moment_;
}


Foam::dimensionedVector& Foam::sixDOFMSDqODE::moment()
{
    return moment_;
}


// ************************************************************************* //

