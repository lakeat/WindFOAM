//Info << "Reading the ship parameters dictionary" << endl,
IOdictionary shipDict
(
IOobject
	(
		"shipDict",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ,	
		IOobject::NO_WRITE
	)
);

wordList motionPatches = shipDict.subDict("bodies").toc();

int numMotionPatches = motionPatches.size();

Info << "motionPatches : " << motionPatches << nl;

