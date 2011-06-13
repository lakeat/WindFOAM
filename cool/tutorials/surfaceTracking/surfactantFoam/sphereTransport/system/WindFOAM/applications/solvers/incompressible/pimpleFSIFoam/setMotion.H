{
    const polyMesh& mesh =
        runTime.lookupObject<polyMesh>(polyMesh::defaultRegion);

    // Grab RBF motion solver
    RBFMotionSolver& ms =
        const_cast<RBFMotionSolver&>
        (
            mesh.lookupObject<RBFMotionSolver>("dynamicMeshDict")
        );

 	pointField statPoints = ms.movingPoints();

//-Note, here is the point, in ideal, the new position of the solid body
//-is obviously a function of external forces and time. However, in current
//-OpenFOAM implementation of sixDOFMSDbodies's solver, it is just the function
//-of time (time here means both the current time, old time, and time step).
//-Here the sixDOFMSDbodies class has its external forces included, actually,
//-it is not sixDOFMSDbodies class has the force included, it is the
//-sixDOFMSDqODE class has the force included. See its H file.
//-So, in order to continue the solid body motion solving process,
//-I hereby update the forces on sixDOFMSDbodies first, before continuing.
//-Daniel Wei
//-06/03/2011

//-Dead force needs to be added, i.e. body weight.
//-Daniel Wei
//-06/10/2011


    dimensionedVector dimensionedTotForce
    (
        "totalForce",
        dimensionSet(1, 1, -2, 0, 0),
        totForce
    );
    dimensionedVector dimensionedTotMoment
    (
        "totalMoment",
        dimensionSet(1, 2, -2, 0, 0),
        totMoment
    );

    structure.getOde()[0].editForce(dimensionedTotForce);
    structure.getOde()[0].editMoment(dimensionedTotMoment);
    structure.solve();

    // Output the ODE solution results to the screen
    Info<< "Current Xrel in time instant = "
        << structure()[0].Xrel().value() << nl
        << "Current Trel in time instant = "
        << structure()[0].Trel().value() << nl
        << "Current X in time instant = "
        << structure()[0].X().value() << nl
        << "Current T in time instant = "
        << structure()[0].T().value() << nl
        << "Current velocity in time instant = "
        << structure()[0].U().value() << nl
        << "Current omega in time instant = "
        << structure()[0].omega().value()  << nl
        << "Current forces in time instant = "
        << structure()[0].force().value() << nl
        << "Current moment in time instant = "
        << structure()[0].moment().value() << nl
        << endl;

    // Motion is a vectorField of all moving boundary points
    vectorField motion(ms.movingPoints().size(), vector::zero);

    vector bodyTranslationVector = structure()[0].Xdiff().value();
    vector bodyRotationVector    = structure()[0].Tdiff().value();

    // Force the 6-DOF to 2-DOF system
    // In a flow around a bridge deck, the deck can only supposed to move in
    // the vertical (y direction) and torsional (around z axis)
    bodyTranslationVector.x() = 0;
    bodyTranslationVector.z() = 0;
    bodyRotationVector.x() = 0;
    bodyRotationVector.y() = 0;

    // Building the rotation matrix
    scalar bRx = bodyRotationVector.x();
    scalar bRy = bodyRotationVector.y();
    scalar bRz = bodyRotationVector.z();

    tensor Rx
        (
            1,              0,              0,
            0, Foam::cos(bRx),-Foam::sin(bRx),
            0, Foam::sin(bRx), Foam::cos(bRx)
        );
    tensor Ry
        (
            Foam::cos(bRy), 0, Foam::sin(bRy),
            0,              1,              0,
           -Foam::sin(bRy), 0, Foam::cos(bRy)
        );
    tensor Rz
        (
            Foam::cos(bRz),-Foam::sin(bRz), 0,
            Foam::sin(bRz), Foam::cos(bRz), 0,
            0,              0,              1
        );
    tensor RRcur = Rx & Ry & Rz;
//    tensor RRold
//        (
//            1,              0,              0,
//            0,              1,              0,
//            0,              0,              1
//        );
//    tensor RR = RRcur - RRold;

    //Info << "tensor bodyRotation = " << RR << endl;
    vectorField bodyRotation = (RRcur & statPoints) - statPoints;


    motion = bodyTranslationVector + bodyRotation;
    //motion = bodyTranslationVector;
    //motion = bodyRotation;


	ms.setMotion(motion);
	//movePoints(ms.newPoints());
}

