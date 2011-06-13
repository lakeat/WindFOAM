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
Setting of initial ship conditions and variables. 
\*---------------------------------------------------------------------------*/
motionPatch_ = motionPatch;
writeCounter = 0;

word fileName__ = motionPatch_ + "_motion.txt";

motionFile = autoPtr<writeMotionFile>
(
	new writeMotionFile (fileName__)
);

odeSolver = autoPtr<ODESolver>(ODESolver::New(ODESolverName, ode));

y = scalarField(6,0);

U_cog=vector::zero;		//velocity of cog in body coordinate system
//U_cog_old(vector::zero), 	//velocity of cog previous time step in body coordinate system
U_cog_global=vector::zero;	//velocity of cog in global coordinate system

Omega_cog=vector::zero;		//rotation speed of cog in body coordinate system
//Omega_cog_old(vector::zero),	//idem for previous time step
Omega_cog_global=vector::zero;	//rotation speed of cog in global coordinate system

translation=vector::zero;	//total translation per time step
rotation=vector::zero;
totalDisplacement=vector::zero;
totalRotation=vector::zero;

aveForce=vector::zero;
aveMoment=vector::zero;
Fnett=vector::zero;
Fstore = pointField(3,vector::zero);
Mstore = pointField(3,vector::zero);

vectX=vector(1,0,0);
vectY=vector(0,1,0);
vectZ=vector(0,0,1);

k_hydroStatTrans=vector::zero;		//spring coefficients which can be updated every time step
k_hydroStatRot=vector::zero;

bodysubDict = dictionary
(
	shipDict.subDict("bodies").subDict(motionPatch_)
);

patchI=label(mesh_.boundaryMesh().findPatchID(motionPatch_));  //patchI: label number of Hull patch

aMax=scalar(readScalar(bodysubDict.lookup("aMax")));
aMin=scalar(-1*aMax);

CoG=vector(bodysubDict.lookup("CoG"));		//intial centre of gravity

// Initialise rotation matrices in case body coordinate system is not aligned with global coordinate system
// This can be useful when e.g. ship is sailing at an angle relative to incoming waves
bodyRotation=vector(bodysubDict.lookup("bodyRotation"));

//create object "global2body" of class "rotationMatrix" for rotating from global to local coord. sys
global2body = rotationMatrix(bodyRotation*pi_/180).rotXYZ();

//create object "body2global" of class "rotationMatrix" for rotating from local back to global coord. sys
body2global = rotationMatrix(-1*bodyRotation*pi_/180).rotXYZ();


//Translation coefficients
calcTranslationDOF=bodysubDict.lookup("calcTranslationDOF");

U_cog_old=bodysubDict.lookup("initialVelocity");

constantForce=bodysubDict.lookup("constantForce");

mass=readScalar(bodysubDict.lookup("mass"));
d=bodysubDict.lookup("linearDamping");
k=bodysubDict.lookup("linearSpring");

//Rotation coeffients
calcRotationDOF=bodysubDict.lookup("calcRotationDOF");

J=bodysubDict.lookup("momentOfInertia");
d_rot=bodysubDict.lookup("linearDamping_rot");
k_rot=bodysubDict.lookup("linearSpring_rot");

Omega_cog_old=bodysubDict.lookup("initialRotationSpeed");
constantMoment=bodysubDict.lookup("constantMoment");

check();

if (Pstream::parRun())
	{
	dirName_ = mesh_.time().path()/".."/"CoGtrack";
	}
else
	{
	dirName_ = mesh_.time().path()/"CoGtrack";
	}

mkDir(dirName_);
fileName_ = dirName_/motionPatch_ + "_" + mesh_.time().timeName() + ".txt";

IFstream CoGtrackStream(fileName_);

if (CoGtrackStream.good())
	{
	restart = false;

	List<vector> CoGtrackStream_ 
	(
		CoGtrackStream
	);
	totalDisplacement += CoGtrackStream_[0];
	CoG += totalDisplacement;
	totalRotation = CoGtrackStream_[1];
	U_cog_old = CoGtrackStream_[2];
	Omega_cog_old = CoGtrackStream_[3];
	Omega_cog_old *= pi_/180;				//convert from degrees/s to radians/s
	Fstore[0] = CoGtrackStream_[4];
	Fstore[1] = CoGtrackStream_[5];
	Mstore[0] = CoGtrackStream_[6];
	Mstore[1] = CoGtrackStream_[7];

	}
else
	{
	restart = true;
	}

// ************************************************************************* //
