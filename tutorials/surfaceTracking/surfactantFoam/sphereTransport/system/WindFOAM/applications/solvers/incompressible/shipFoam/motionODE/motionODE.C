/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

Application
Description
\*---------------------------------------------------------------------------*/
#include "motionODE.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

motionODE::motionODE()
    {}

label motionODE::nEqns() const
    {
        return 6;
    }

void motionODE::derivatives
    (
        const scalar x,
        const scalarField& y,
        scalarField& dydx
    ) const
    {
        dydx[0] = y[1];
	dydx[1] = (y[5]-y[3]*y[1]-y[4]*y[0]) / y[2];
//	m*y" + b*y' + k*y = F
//	y=y[0]
//	substitution:  y' = y[1]
//	this yields: y" = (F - b*y[1] - k*y) / m
//	constants: y[2] = m; y[3] = b; y[4] = k; y[5] = F
	dydx[2] = 0;			
	dydx[3] = 0;			
	dydx[4] = 0;			
	dydx[5] = 0;			
    }

void motionODE::jacobian
    (
        const scalar x,
        const scalarField& y,
        scalarField& dfdx,
        scalarSquareMatrix& dfdy
    ) const
    {
    }
