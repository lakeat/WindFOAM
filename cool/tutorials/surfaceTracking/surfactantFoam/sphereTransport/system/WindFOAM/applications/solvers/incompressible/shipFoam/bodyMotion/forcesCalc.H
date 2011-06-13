/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software, you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY, without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM, if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
calculates pressure forces and moments on motion patch
move previous pressureForces back in pForce
\*---------------------------------------------------------------------------*/

const volScalarField& p_ = mesh_.lookupObject<volScalarField>("p");

scalar sumWF(0);
for (int ii=2; ii>0; ii--)
{	
Fstore[ii] = Fstore[ii-1];
Mstore[ii] = Mstore[ii-1];
}

for (int ii=0;ii<3;ii++)
{
sumWF += weightFactor[ii];
}

Fstore[0] =
gSum
(
	p_.boundaryField()[patchI]*
	(global2body & mesh_.Sf().boundaryField()[patchI])
);
//Explanation: 
//gSum = sum but when running in parallel, uses all processors
// p value at every patch of the boundaryField multiplied with:
//mesh.Sf()... = patch area vector (projected areas in three directions)
//inner product global2body.rotXYZ(): rotation matrix which takes account for local coordinate system
//in cases the body has a different orientation than global coordinate system. Yet to be checked!

//update position of centre of gravity
CoG += mesh_.time().deltaT().value() * U_cog;  

//Moment around centre of gravity
Mstore[0] =
gSum
(
	(global2body & (mesh_.Cf().boundaryField()[patchI] - CoG))
	^ (
	p_.boundaryField()[patchI]*
	(global2body & mesh_.Sf().boundaryField()[patchI])
	)
);

for (int jj=0; jj<3; jj++)
{
	aveForce[jj] = 0;
	aveMoment[jj] = 0;
	for (int ii=0; ii<3; ii++)
	{
		aveForce[jj] += weightFactor[ii] * Fstore[ii][jj] / sumWF;
		aveMoment[jj] += weightFactor[ii] * Mstore[ii][jj] / sumWF;
	}
}


Fnett = aveForce + (global2body & (mass * g_.value()));

Info <<"forcesCalc on body: " << motionPatch_ << nl;
Info <<"actual pressure force : " << Fstore[0] << nl;
Info <<"Calculated moment around CoG : " << Mstore[0] << nl << nl;
