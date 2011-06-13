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
Calculate spring terms by moving the motion patch with 0.01 m and 
recalculating p field. Next: calculate new pressure force on the motion patch.
Difference with original p-force divided by 0.01 is spring term
Fp0 - Fp1 = k * dx

\*---------------------------------------------------------------------------*/

//scalar dX(0.01);
//vector dRot(0.02, 0.02, 0.02);		//angle to rotate in radians, appr. 1 degree

//dimensionedVector dist("dist", dimensionSet(0,1,0,0,0,0,0),vector(dX, dX, dX));
//dimensionedVector dimensionedCoG("dimensionedCoG", dimensionSet(0,1,0,0,0,0,0), CoG);

//rotationMatrix XYZ(dRot);		//creates object "rotXYZ" of class "rotationMatrix"
dimensionedVector dist("dist", dimensionSet(0,1,0,0,0,0,0),vector(0.1, 0.1, 0.1));

volScalarField gh_1("gh", g_ & (mesh_.C()+dist));
const volScalarField& rho_ = mesh_.lookupObject<volScalarField>("rho");
 
volScalarField p_1 = rho_*gh_1;		//question: needs rho to be updated too? Denk het niet

vector pForce_1 = 
gSum
(
	p_1.boundaryField()[patchI]*
	(global2body & mesh_.Sf().boundaryField()[patchI])
);

//Explanation: 
//gSum = sum but when running in parallel, uses all processors
// p value at every patch of the boundaryField multiplied with:
//mesh.Sf()... = patch area vector (projected areas in three directions)
//inner product global2body.rotXYZ(): rotation matrix which takes account for local coordinate system
//in cases the body has a different orientation than global coordinate system. Yet to be checked!

vector dFp = Fstore[0] - pForce_1;	//plus of min?

k_hydroStatTrans =  dFp / scalar(0.1);

//tensor rotations = rotationMatrix(dRot).rotXYZ();
dimensionedVector dimensionedCoG("dimensionedCoG", dimensionSet(0,1,0,0,0,0,0), CoG);

volScalarField gh_2("gh", g_ & (rotations & (mesh_.C()-dimensionedCoG)));
//const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
 
volScalarField p_2 = rho_*gh_2;		

vector pMoment_1 = 
gSum
(
	(global2body & (mesh_.Cf().boundaryField()[patchI] - CoG))
	^ (
	p_2.boundaryField()[patchI]*
	(global2body & mesh_.Sf().boundaryField()[patchI])
	)
);

vector dMp = Mstore[0] - pMoment_1;	//Mstore[0] = vector of current timestep

k_hydroStatRot =  dMp / scalar(0.2);	//Moment per radian rotation

Info << "new spring terms for translation = " << k_hydroStatTrans << nl;
Info << "new spring terms for rotation = " << k_hydroStatRot << nl;


