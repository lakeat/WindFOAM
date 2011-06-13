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

    //-Note, in the following code, since we just have one moving body here,
    // so structure()[0].memberFunction is used to call "sixDOFMSDqODE" class.

