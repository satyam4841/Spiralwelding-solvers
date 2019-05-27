/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    icoFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

    \heading Solver details
    The solver uses the PISO algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}}
          + \div \left( \vec{U} \vec{U} \right)
          - \div \left(\nu \grad \vec{U} \right)
          = - \grad p
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"
#include "reader.h"
#include <vector>
#include <string>
#include <fstream>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //#include "addCheckCaseOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //==================================Constants====================================
    volVectorField C = mesh.C();
    double deltaTime=runTime.deltaT().value();
    double totalTime=runTime.endTime().value();
    double totalSteps=totalTime/deltaTime;
    int flagCell=0;
    Info<< "Time of simulation = " << totalTime <<"s"<<nl << endl;
    Info<< "Total timesteps: = " << totalSteps << nl << endl;
        
//====================================================================================      

        int nline=0;
        std::string line;        
        std::ifstream myfile;
        myfile.open("/home/satyam/SpiralData/spiraldata2.txt");
        
        if(myfile.is_open()){
            while(!myfile.eof()){
                getline(myfile,line);
                nline++;
            }
            myfile.close();
        }

        cout<<"number of coordinates :"<<nline<<"\n";    
        std::vector<Cartesian> spiral(nline); 
        int t=0;
        std::fstream inputfile;
        inputfile.open("/home/satyam/SpiralData/spiraldata2.txt");
        
        while(! inputfile.eof()){
            inputfile>>spiral[t].x>>spiral[t].y;
            spiral[t].z=0.05;
            t++;
        }
        inputfile.close();

        std::vector<int> cellii(nline);
        for (int t=0;t<nline;t++){
            vector position(spiral[t].x, spiral[t].y, 0.0047);
            cellii[t] = mesh.findCell(position);
            //Info<< "spiralcelli = " << cellii[t]<<"\t For point : " << t << nl << endl;
        }
            
        std::vector<int> spiralcells;
        for (int i=1;i<cellii.size();i++){
            if ((cellii[i-1]!=cellii[i]) || i==1){
                    spiralcells.push_back(cellii[i]);
            }
        }
        
        Info<< "Total input:" << cellii.size() <<"\t" << "Unique: " << spiralcells.size() << nl << endl;
        for (int i=0;i<spiralcells.size();i++){
            cout<<spiralcells[i]<<"\t";   
        }
        

//====================================================================================      

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        //Info<< "Time = " << runTime.timeName().value() << nl << endl;

        #include "CourantNo.H"

        // Momentum predictor

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
            Info<< "continuity ok."<< nl << endl;
        }
//===================================================================

        //double currentTime=runTime.time().value();
        double currentTime=runTime.time().value();
        int currentStep=currentTime/deltaTime;
        int spiCell=spiralcells.size();
        int duration=totalSteps/spiCell;
        double ro = r.value();
        double Q = q.value(); 

        Info<< "Duration at each cell:"<< duration << endl;
        Info<< "Gaussian radius:"<< ro << endl;
        if(flagCell==spiCell-10){
                g[spiralcells[flagCell]]=Q;
        }
        
        g[spiralcells[flagCell]]=Q;
        
        if(currentStep%duration==0){
            
            g[spiralcells[flagCell]]=Q;  
            
            forAll (C, cellI){  
                    vector cellCentroid = C[cellI];
                    vector X = C[spiralcells[flagCell]];
                    const double x = (cellCentroid.x() - X.x())*(cellCentroid.x() - X.x());
                    const double y = (cellCentroid.y() - X.y())*(cellCentroid.y() - X.y());
                    const double z = (cellCentroid.z() - X.z())*(cellCentroid.z() - X.z());
                    g[cellI]=(Q/(3.14*ro*ro))*Foam::exp((x+y+z)*(1/(-ro*ro)));
            }
            
            flagCell++;
            
            forAll (C, cellI) //To make last heat source in last cell zero.
            {
                g[spiralcells[flagCell]]=0;
            }
        
        }

        Info<< "Current cell:"<< spiralcells[flagCell] << endl;
//====================================================================

        fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT, T)
             ==
                g
            );
        
            TEqn.solve();

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
