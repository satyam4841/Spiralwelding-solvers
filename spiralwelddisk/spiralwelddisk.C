/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    laplacianFoam

Group
    grpBasicSolvers

Description
    Laplace equation solver for a scalar quantity.

    \heading Solver details
    The solver is applicable to, e.g. for thermal diffusion in a solid.  The
    equation is given by:

    \f[
        \ddt{T}  = \div \left( D_T \grad T \right)
    \f]

    Where:
    \vartable
        T     | Scalar field which is solved for, e.g. temperature
        D_T   | Diffusion coefficient
    \endvartable

    \heading Required fields
    \plaintable
        T     | Scalar field which is solved for, e.g. temperature
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "Cartesian.h"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addCheckCaseOptions.H"
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    volVectorField C = mesh.C();
    double deltaTime=runTime.deltaT().value();
    double totalTime=runTime.endTime().value();
    double totalSteps=totalTime/deltaTime;
    int flagCell=0;
    Info<< "Time of simulation = " << totalTime <<"s"<<nl << endl;
    Info<< "Total timesteps: = " << totalSteps << nl << endl;
        
    int nline=0;
        std::string line;        
        std::ifstream myfile;
        myfile.open("/home/satyam/SpiralData/spiraldata30May.txt");
        
        if(myfile.is_open())
        {
            while(!myfile.eof())
            {
                getline(myfile,line);
                nline++;
            }
            myfile.close();
        }

        cout<<"number of coordinates :"<<nline<<"\n";    
        std::vector<Cartesian> spiral(nline); 
        int t=0;
        std::fstream inputfile;
        inputfile.open("/home/satyam/SpiralData/spiraldata30May.txt");
        
        while(! inputfile.eof())
        {
            inputfile>>spiral[t].x>>spiral[t].y;
            spiral[t].z=0.05;
            t++;
        }
        inputfile.close();

        std::vector<int> cellii(nline);
        for (int t=0;t<nline;t++)
        {
            vector position(spiral[t].x, spiral[t].y, -0.0001);
            cellii[t] = mesh.findCell(position);
            //Info<< "spiralcelli = " << cellii[t]<<"\t For point : " << t << nl << endl;
        }
            
        std::vector<int> spiralcells;
        for (int i=1;i<cellii.size();i++)
        {
            if ((cellii[i-1]!=cellii[i]) || i==1)
            {
                    spiralcells.push_back(cellii[i]);
            }
        }
        
        Info<< "Total input:" << cellii.size() <<"\t" << "Unique: " << spiralcells.size() << nl << endl;
        
        /*for (int i=0;i<spiralcells.size();i++)
        {
            cout<<spiralcells[i]<<"\t";   
        }*/
            
        int spiCell=spiralcells.size();
        int duration=totalSteps/spiCell;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;
    
    Info<< "Duration at each cell:"<< duration << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        //double currentTime=runTime.time().value();
        int currentStep=runTime.time().value()/deltaTime;
        //int spiCell=spiralcells.size();
        //int duration=totalSteps/spiCell;
        double ro = r.value();
        double Q = q.value(); 

        Info<< "Gaussian radius:"<< ro << endl;
        if(flagCell==spiCell-10)
        {
                g[spiralcells[flagCell]]=Q;
        }
        
        g[spiralcells[flagCell]]=Q;
        
        if(currentStep%duration==0)
        {
            
            g[spiralcells[flagCell]]=Q;  
            
            forAll (C, cellI)
            {  
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

       // Info<< "Current cell:"<< spiralcells[flagCell] << endl;

// =========================================================================================  //
        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT, T)
             ==
                g
            );

            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }

        runTime.write();
        //#include "write.H"  // Uncomment to obtain gradT 

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
