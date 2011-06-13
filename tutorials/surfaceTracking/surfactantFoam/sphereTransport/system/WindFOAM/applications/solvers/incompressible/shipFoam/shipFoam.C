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
    shipFoam

writen by: M. Couwenberg, febr 2009

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.
shipFoam is based on interDyMFoam, extended with the possibility to move boundary patches based in hydrostatic forces/moments and additional external forces/moments.
At every time step the following ordinary differntial equation is calculated:
F = mass (/inertia) * x" + damping * x' + spring * x

with:
F = external force (or moment) + gravity force + pressure force. Although turbulence forces are also calculated and printed in the log file, these are not added.
mass / inertia: speaks for itself
x" = accelleration of DOF x (which can be 3 translations and 3 rotations)
damping: speaks for itself
x' = velocity of DOF x. This velocity is mapped on the field pointMotionU on the respective motionPatch.
spring: external spring as supplied in the shipDict + hydrostatic spring term.
Hydrostatic spring is calculated by using the actual pressure force and recalculate this pressure force with a slightly translated motionPatch, which is similar to pushing the motionpatch a small amount in a certain direction and measuring the difference in forces.
x = translation / rotation DOF. Depending on user input, any of the 6 DOF's (or all) can be calculated 

Apart from the motion classes, a small functionality is added compared with interDyMFoam:
one can enter a relaxationFactor for pd. This is important to obtain stable pressure results when the bodies are moved.

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"
#include "probes.H"

#include "bodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "readPISOControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"
    #include "correctPhi.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

#include "countMotionPatches.H"
//create object and construct
bodyMotion body1(mesh, motionPatches, g);
//create object and copy construct
bodyMotion body2(body1);
bodyMotion body3(body1);
bodyMotion body4(body1);

//bodyMotion*() bodies = new bodyMotion[numMotionPatches]();
//forAll (motionPatches, mPI)
//{
//bodies = new bodyMotion[mPI] (mesh, motionPatches[mPI]);
//}
body1.initialize(motionPatches[0]);
if(numMotionPatches > 1)
{Info <<"Initialize body 2 ..." << nl;
body2.initialize(motionPatches[1]);}
if(numMotionPatches > 2)
{Info <<"Initialize body 3 ..." << nl;
body3.initialize(motionPatches[2]);}
if(numMotionPatches > 3)
{Info <<"Initialize body 4 ..." << nl;
body4.initialize(motionPatches[3]);}


scalar counter(0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;


    while (runTime.run())
    {
        #include "readControls.H"
        #include "CourantNo.H"


        // Make the fluxes absolute
        fvc::makeAbsolute(phi, U);

        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

        // Do any mesh changes
        mesh.update();

        if (mesh.changing())
        {
            Info<< "Execution time for mesh.update() = "
                << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                << " s" << endl;
        }

//        volScalarField gh("gh", g & mesh.C());
//        surfaceScalarField ghf("ghf", g & mesh.Cf());

        if (mesh.changing() && correctPhi)
        {
            #include "correctPhi.H"
        }

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);

        if (mesh.changing() && checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
        }

        twoPhaseProperties.correct();

        #include "alphaEqnSubCycle.H"

        #include "UEqn.H"

        // --- PISO loop
        for (int corr=0; corr<nCorr; corr++)
        {
            #include "pEqn.H"
        }

        turbulence->correct();

	body1.forcesCalc();
	if(numMotionPatches > 1)
	{body2.forcesCalc();}
	if(numMotionPatches > 2)
	{body3.forcesCalc();}
	if(numMotionPatches > 3)
	{body4.forcesCalc();}

	if (counter >= body1.startUpdate())
	{
		body1.forceBalance();
		if(numMotionPatches > 1)
		{body2.forceBalance();}
		if(numMotionPatches > 2)
		{body3.forceBalance();}
		if(numMotionPatches > 3)
		{body4.forceBalance();}
	}
	else
	{
	counter++;
	}
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
