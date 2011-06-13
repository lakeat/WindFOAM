    Info << "Reading FSI properties" << endl;
    IOdictionary FSIProperties
    (
        IOobject
        (
            "FSIProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    sixDOFMSDbodies structure(runTime);



//    // Annoying warnings removal
//    transonic=0;
//    nOuterCorr=0;
//    nCorr=0;
//    nNonOrthCorr=0;
//    transonic=0;
//    momentumPredictor=0;
//    maxCo=0;
//    maxDeltaT=0;

