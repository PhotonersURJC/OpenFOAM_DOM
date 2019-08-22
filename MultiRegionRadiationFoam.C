/*---------------------------------------------------------------------------*\
             Discrete Ordinate Method Radiation Model for OpenFOAM.


Code corresponding to the article entitled:


"Validation of a Novel Comprehensive Open Source Discrete Ordinate Method

Code for Radiation Transport Simulation"


by:


José Moreno, Cintia Casado, Javier Marugán

Department of Chemical and Environmental Technology,

ESCET, Universidad Rey Juan Carlos,

C/Tulipán s/n, 28933 Móstoles (Madrid), Spain

Tel. +34 91 664 7466; E-mail: javier.marugan@urjc.es

 

submitted for publication to the journal:


Computer Methods in Applied Mechanics and Engineering


on the 5th of March, 2018.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "newRadiationModel.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "fvm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    regionProperties rp(runTime);
    #include "createFluidMeshes.H"
	#include "createSolidMeshes.H"
    bool Converged = false;
	solverPerformance::debug=0;
    while (runTime.loop())
    {
        Converged = true;
        Info<< "Time = " << runTime.timeName() << nl << endl;
        forAll(fluidRegions, i)
        {
            Info<< "\nSolving for fluid region "
                << fluidRegions[i].name() << endl;
			
	    radiation::newRadiationModel& rad = fluidRadiation[i];
            if (!rad.correct()) {Converged=false;}
        }
        forAll(solidRegions, i)
        {
            Info<< "\nSolving for solid region "
                << solidRegions[i].name() << endl;
	    radiation::newRadiationModel& rad = solidRadiation[i];
            if (!rad.correct()) {Converged=false;}
        }
		if(Converged)
			runTime.writeAndEnd();
        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
