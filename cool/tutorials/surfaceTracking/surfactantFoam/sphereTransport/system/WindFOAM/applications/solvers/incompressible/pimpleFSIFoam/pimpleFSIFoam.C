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

Application
    pimpleFSIFoam.C

Description
    This is a Fluid Structure Interaction solver, this could be done
    theoretically in a library level as well as in a solver level.
    However, considering the complexity, I just make it in a solver application.

    The Solver codes is grabbed from some different code, the basic structure
    is based on the pimpleDyMFoam, so to make it easy to suite for the moving
    mesh case.
    The mesh is moving on the basis of a RBF solver, which is simple enough to
    handle the problem we meet.
    There online is a icoFsiFoam, which I think is helpful in its structure,
    however, in our case , the body is rigid, no pointPatch to pointPatch
    intersection is needed, no body-shape solver is needed.
    And since the body motion is a basically 6-DoF motion, (in our bridge case,
    it is actually just a 2-DoF motion, a heaving motion together with a
    rotating motion). So the motion is solved by a sixDOFMSDbodies solver that is
    already built in the OpenFOAM-1.6-ext. However, in order to make it useful
    for our case, some minor modification to the sixDOFMSDqODE and sixDOFMSDbodies
    classes, to allow the change of the forces from external.

Updates
    06/07/2011
    ----------
    -Currently, the dead force, like gravity is not considered yet.
    -The rotation motion has not worked very well. This needs to be further
    revisited.
    -Should it be forces=fTmp, or forces+=fTmp???

Authorship
    Daniel Wei, NatHaz Lab.  All rights reserved.
    06/07/2011

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "dynamicFvMesh.H"
#include "RBFMotionSolver.H"
#include "ODESolver.H"
#include "sixDOFMSDbodies.H"
#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "readPIMPLEControls.H"
#   include "initContinuityErrs.H"
#   include "createFields.H"
#   include "readTimeControls.H"
#   include "readFSIProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    fileName bodyMotionDir;
    if (Pstream::master())
    {
        if (Pstream::parRun())
        {
            bodyMotionDir = runTime.path()/"..";
        }
        else
        {
            bodyMotionDir = runTime.path();
        }
    }
    OFstream of(bodyMotionDir/"bodyMotion.dat");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

//    while (runTime.run())
//    {
    for (runTime++; !runTime.end(); runTime++)
    {
#       include "readControls.H"
#       include "CourantNo.H"

        // Make the fluxes absolute
        fvc::makeAbsolute(phi, U);

#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        bool meshChanged = mesh.update();

        if (correctPhi && (mesh.moving() || meshChanged))
        {
#           include "correctPhi.H"
        }

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);

        if (mesh.moving() && checkMeshCourantNo)
        {
#           include "meshCourantNo.H"
        }

#       include "setPressure.H"
#       include "setMotion.H"
#       include "solveFluid.H"

        of << runTime.value() << tab
           << structure()[0].X().value().y() << tab
           << structure()[0].T().value().z() << endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

