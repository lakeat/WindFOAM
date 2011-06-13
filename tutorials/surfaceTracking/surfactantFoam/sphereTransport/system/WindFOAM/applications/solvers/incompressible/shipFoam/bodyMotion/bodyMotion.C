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
\*---------------------------------------------------------------------------*/

#include "bodyMotion.H"

namespace Foam
{
//constructor
bodyMotion::bodyMotion
(
	const fvMesh& mesh,
	const wordList motionPatches,
	const uniformDimensionedVectorField& g
)
:
	pi_(mathematicalConstant::pi),
	mesh_(mesh),
	shipDict		//create object shipDict of class IOdictionary
	(
	IOobject
		(
			"shipDict",
			mesh_.time().constant(),
			mesh_,
			IOobject::MUST_READ,	
			IOobject::NO_WRITE
		)
	),
	weightFactor(vector (shipDict.lookup("weightFactor"))),

	//this is the field which obtains the velocities to move the motionpatch
	pointMotionU 
	(
	const_cast<pointVectorField&>
	(
		mesh_.objectRegistry::lookupObject<pointVectorField>
		(
			"pointMotionU"
		)
	)),

	ODESolverName(shipDict.subDict("ODECoeffs").lookup("ODESolver")),
	eps(readScalar(shipDict.subDict("ODECoeffs").lookup("eps"))),
	hEst(readScalar(shipDict.subDict("ODECoeffs").lookup("hEst"))),
	
	g_(g.value()),
	motionPatches_(motionPatches),
	writeCoG(Switch(shipDict.lookup("writeCoG"))),
	writeInterval(readInt(shipDict.lookup("writeInterval"))),
	springCoeffUpdateInterval(int(5)),
	springCoeffUpdateIntervalCounter(int(5))

{
    if (shipDict.found("springCoeffUpdateInterval"))
    {
        springCoeffUpdateInterval = readInt(shipDict.lookup("springCoeffUpdateInterval"));
	springCoeffUpdateIntervalCounter = springCoeffUpdateInterval;
	Info << "Spring coefficients are updated every " << springCoeffUpdateInterval << " timesteps. " << nl;
    }
    else
    {
	Info << "Spring coefficients are updated every 5 timesteps. " << nl;
    }
	
	vector dRot(0.02, 0.02, 0.02);		//angle to rotate in radians, appr. 1 degree
	rotations = rotationMatrix(dRot).rotXYZ();		//creates object "rotXYZ" of class "rotationMatrix"
}

void bodyMotion::initialize(word motionPatch)
{
#include "initialize.H"
}

scalar bodyMotion::startUpdate()
{

	if (restart == false)
		{
		return 0;		//startUpdate delay not to be used after a restart => 0
		}
	else
		{
		return(readScalar(shipDict.lookup("startUpdate")));
		}
}

void bodyMotion::forcesCalc()
{
#include "forcesCalc.H"
}

void Foam::bodyMotion::forceBalance()
{
#include "forceBalance.H"
}

void bodyMotion::check()
{
bool fail = true;
forAll(mesh_.boundaryMesh().names(), nameI)
	{
		if (nameI == motionPatch_)
		{
		fail = false;
		}
	}

if (fail==true)
	{
	FatalErrorIn("bodyMotion::check() const")
		<< "motionPatch named " << motionPatch_
		<< " not found.  Available patches: "
		<< mesh_.boundaryMesh().names()
		<< abort(FatalError);
	}
}

} //end namespace Foam
