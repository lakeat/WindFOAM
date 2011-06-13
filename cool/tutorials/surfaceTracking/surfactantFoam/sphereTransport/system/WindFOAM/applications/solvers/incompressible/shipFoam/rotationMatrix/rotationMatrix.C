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

#include "rotationMatrix.H"

//constructor
rotationMatrix::rotationMatrix(vector bodyRotation)
{
	scalar radX(bodyRotation[0]);
	scalar radY(bodyRotation[1]);
	scalar radZ(bodyRotation[2]);
	
	Rx=tensor
	(
	1, 0, 0,
	0, Foam::cos(radX), -Foam::sin(radX), 
	0, Foam::sin(radX), Foam::cos(radX)
	);

	Ry=tensor
	(
	Foam::cos(radY), 0, Foam::sin(radY), 
	0, 1, 0, 
	-Foam::sin(radY), 0, Foam::cos(radY)
	);

	Rz=tensor
	(
	Foam::cos(radZ), -Foam::sin(radZ), 0, 
	Foam::sin(radZ), Foam::cos(radZ), 0, 
	0, 0, 1
	);
}
//empty constructor
rotationMatrix::rotationMatrix()
{}


//member function
tensor rotationMatrix::rotXYZ()		//returns rotation matrix for 3 eulerian angles
	{
	return (Rz & Ry & Rx);
	}

