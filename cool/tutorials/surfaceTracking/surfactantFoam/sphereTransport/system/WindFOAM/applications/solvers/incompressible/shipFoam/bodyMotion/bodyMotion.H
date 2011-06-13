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

#include "IFstream.H"
#include "word.H"
#include "fvCFD.H"
#include "fvMesh.H"
//#include "turbulenceModel.H"
#include "pointFields.H"
#include "motionODE.H"
#include "rotationMatrix.H"
#include "writeMotionFile.H"

namespace Foam
{

class bodyMotion
{
scalar pi_;
//- Reference to mesh
const fvMesh& mesh_;
IOdictionary shipDict;
vector weightFactor;
pointVectorField& pointMotionU;

word ODESolverName;
scalar eps;
scalar hEst;

dimensionedVector g_;
wordList motionPatches_;
int writeCounter;
bool writeCoG;
int writeInterval;
int springCoeffUpdateInterval, springCoeffUpdateIntervalCounter;
word motionPatch_;
fileName dirName_, fileName_;
scalarField y;
motionODE ode;
autoPtr<ODESolver> odeSolver;
autoPtr<writeMotionFile> motionFile;

tensor global2body, body2global; 
vector U_cog, U_cog_old, U_cog_global;		
vector Omega_cog, Omega_cog_old, Omega_cog_global;
vector translation, rotation, totalDisplacement, totalRotation;
vector aveForce, aveMoment, Fnett;
pointField Fstore, Mstore;
vector vectX, vectY, vectZ;
	
dictionary bodysubDict;

label patchI;  //patchI: label number of Hull patch
scalar aMax, aMin;
vector CoG;
vector bodyRotation;
vector calcTranslationDOF, calcRotationDOF;
vector constantForce;
scalar mass;
vector d, k;
vector J, d_rot, k_rot;
vector constantMoment;
	
pointField pointCentres;
vectorField Xdot;
vectorField Ydot;
vectorField Zdot;
List<vector> CoGtrackStream_ ;
bool restart;
vector k_hydroStatTrans, k_hydroStatRot;
vector dRot;
tensor rotations;		
// ************************************************************************* //
public:

//constructor
bodyMotion
(    
	const fvMesh& mesh,
	const wordList motionPatches,
	const uniformDimensionedVectorField& g
);

//member functions
void initialize(word motionPatch);
scalar startUpdate();
void forcesCalc();
void forceBalance();
void check();
}; 

} //end namespace Foam
