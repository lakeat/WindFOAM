{
    // Setting pressure on solid patch
    Info << "Setting force and moment on the moving body" << endl;

//    // Read fluid patch data
//    word movingPatchName(FSIProperties.lookup("movingPatch"));

//    label movingPatchID =
//        mesh.boundaryMesh().findPatchID(movingPatchName);

    forces bodyForces
        (
            "bodyForces",
            mesh,
            FSIProperties,
            0
        );
    forcesMoments fm = bodyForces.calcForcesMoment();

    totForce = fm.first().first() + fm.first().second();
    totMoment = fm.second().first() + fm.second().second();

    Info << "Total forces  = " << totForce << endl;
    Info << "Total moments = " << totMoment << endl;

    //-For our bridge deck cases, it is just a 2-DOF system
    totForce.x()  = 0;
    totForce.z()  = 0;
    totMoment.x() = 0;
    totMoment.y() = 0;


    dimensionedScalar bodyWeight
    (
        FSIProperties.lookup("bodyWeight")
    );

    totForce.y()  -= bodyWeight.value();


/*
    if (movingPatchID < 0)
    {
        FatalErrorIn("void RBFMotionSolver::makeControlIDs()")
            << "Patch " << movingPatches_[patchI] << " not found.  "
            << "valid patch names: " << mesh().boundaryMesh().names()
            << abort(FatalError);
    }
*/

/*
    //-Just use the totalPressureForce is not enough, so I changed it.
    //-Daniel Wei
    //-06/06/2011
    vector totalPressureForce =
        sum
        (
            p.boundaryField()[movingPatchID]*
            mesh.Sf().boundaryField()[movingPatchID]
        );
    totalPressureForce *= rhoFluid.value();
    totalPressureForce.x()=0;
    totalPressureForce.z()=0;

    Info << "Total pressure force = " << totalPressureForce << endl;
*/


    //-Just use the totalPressureForce is not enough, so I changed it.
    //-Daniel Wei
    //-06/06/2011
}

