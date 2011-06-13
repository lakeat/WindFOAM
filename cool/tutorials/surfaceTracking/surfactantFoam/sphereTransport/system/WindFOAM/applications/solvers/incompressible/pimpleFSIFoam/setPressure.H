{
    Info << "Calculating forces and moments on the moving body" << endl;

    // Tuple which holds the pressure (.first()) and viscous (.second) forces
    typedef Tuple2<vector, vector> pressureViscous;

    // Tuple which holds the forces (.first()) and moment (.second)
    // pressure/viscous forces Tuples.
    typedef Tuple2<pressureViscous, pressureViscous> forcesMoments;

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

//    Info << "Total forces  = " << totForce  << nl
//         << "Total moments = " << totMoment << endl;

    //-For bridge deck simulation, it is just a 2-DOF system,
    // so all the other four DOF is constrained.
    totForce.x()  = 0;
    totForce.z()  = 0;
    totMoment.x() = 0;
    totMoment.y() = 0;


    //-Body weight, currently this is input automatically.
    dimensionedScalar bodyWeight
    (
        "bodyWeight",
        dimensionSet(1, 1, -2, 0, 0),
        9.8*structure()[0].mass().value()
    );

    totForce.y() -= bodyWeight.value();


    //-Just use the totalPressureForce is not enough, so I changed it.
    //-Daniel Wei
    //-06/06/2011
//    vector totalPressureForce =
//        sum
//        (
//            p.boundaryField()[movingPatchID]*
//            mesh.Sf().boundaryField()[movingPatchID]
//        );
//    totalPressureForce *= rhoFluid.value();
//    totalPressureForce.x()=0;
//    totalPressureForce.z()=0;
//    Info << "Total pressure force = " << totalPressureForce << endl;

}

