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
 Description
    Solves the force balance equation to provide the velocity of motion
    for the next time step.
*   Uses a generalised setup which solves the basic equation:
     m*x'' + c*x' + k*x = Fres
Only pressure forces and gravity forces are taken into account.
\*---------------------------------------------------------------------------*/

//Calculate body forces

Info<< "forceBalance for body "<< motionPatch_ << endl;

totalDisplacement +=mesh_.time().deltaT().value() * U_cog_old;	//uses calculated velocity which is calculated in previous timeStep
totalRotation +=rotation*180/pi_;

//Loop to calculate 3 DOF translation velocities
y[2] = mass;	//mass

if (springCoeffUpdateIntervalCounter >= springCoeffUpdateInterval)
{
	# include "calcSpringCoeff.H"
	springCoeffUpdateIntervalCounter = 0;
}
else
{
	springCoeffUpdateIntervalCounter++;
}

for (int ii=0; ii<3; ii++)
{	
	if (calcTranslationDOF[ii]==int(1))
	{
	y[0] = 0;				//reset initial displacement for every timestep assumes static equilibrium
	y[1] = U_cog_old[ii]; 
	y[3] = d[ii];				//velocity dampi_ng term b
	y[4] = mag(totalDisplacement[ii])*k[ii] + k_hydroStatTrans[ii];	//spring displacement term k
	y[5] = Fnett[ii] + constantForce[ii];	//Force term F
	# include "motionCalc.H"
	U_cog[ii] = y[1];			//y[1] is velocity
	translation[ii] = y[0];			//y[0] is displacement term
	
	//limit motion based on maximum acceleration aMax
	U_cog[ii] = min(aMax * mesh_.time().deltaT().value() + U_cog_old[ii], U_cog[ii]);
	U_cog[ii] = max(aMin * mesh_.time().deltaT().value() + U_cog_old[ii], U_cog[ii]);
	if (y[1] != U_cog[ii])
		{
		Info << "Limiting CoG motion on maximum accelleration from " << y[1] << " to " << U_cog[ii] << nl;
		}
	}
}
Info << "U_cog [m/s]= " << U_cog << nl;
Info << "total displacement [m] = " << totalDisplacement << nl;
U_cog_old = U_cog;
U_cog_global = body2global & U_cog;	//rotate local velocities to global velocities

//Loop to calculate 3 DOF rotation velocities
for (int ii=0; ii<3; ii++)
{	
	if (calcRotationDOF[ii]==int(1))
	{
	y[0] = 0;		//reset initial displacement for every timestep assumes static equilibrium
	y[1] = Omega_cog_old[ii]; 
	y[2] = J[ii];					//moment of inertia
	y[3] = d_rot[ii];				//velocity dampi_ng term b
	y[4] = mag(totalRotation[ii]*(pi_/180))*k_rot[ii] + k_hydroStatRot[ii];		//spring displacement term k
	y[5] = aveMoment[ii] + constantMoment[ii];	//Force term F

	# include "motionCalc.H"
	Omega_cog[ii] = y[1];				//y[1] is rotation velocity
	rotation[ii] = y[0];				//y[0] is rotation term
	
	}
}

Info << "Omega_cog [degrees/s]= " << Omega_cog*180/pi_ << nl;
Info << "total rotation [degrees]= " << totalRotation << nl;

Omega_cog_old = Omega_cog;
Omega_cog_global = body2global & Omega_cog;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
pointCentres = mesh_.boundaryMesh()[patchI].localPoints() - CoG;

//Explanation for 1st line (Xdot)
//pointCentres ^ vectX: cross product between two vectors yield vector, in this example
//: (0,z,-y)
//Next: this vector & Omega_cog: inner product, yields scalar(Field)
//Finally: scalar(Field) * vectX: scalar multiplication yields vector(Field), :(x,0,0)
Xdot = (((pointCentres ^ vectX) & Omega_cog_global) * vectX);
Ydot = (((pointCentres ^ vectY) & Omega_cog_global) * vectY);
Zdot = (((pointCentres ^ vectZ) & Omega_cog_global) * vectZ);

//Do not forget the double ==!
//Map the calculated velocities on the motion patch
pointMotionU.boundaryField()[patchI] == 	
(
Xdot + Ydot + Zdot +
U_cog_global
);

if (mesh_.time().write())
{
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
	
	OFstream writeCoGtrack
		(
		fileName_
		);
	writeCoGtrack << "8" << nl << "(" << nl;
	writeCoGtrack << totalDisplacement << nl << totalRotation << nl << U_cog_old << nl << Omega_cog_old*180/pi_ << nl;
	writeCoGtrack << Fstore[0] << nl << Fstore[1] << nl << Mstore[0] << nl << Mstore[1] << nl;
	writeCoGtrack << ")" << nl;
}

if (writeCoG)
	{
		writeCounter++;
		
	if (writeCounter >= writeInterval)
		{
		motionFile->write(mesh_.time().timeName(), totalDisplacement, totalRotation);
		writeCounter = 0;
		}
	}

Info <<"end forceBalance" << nl<<nl;
