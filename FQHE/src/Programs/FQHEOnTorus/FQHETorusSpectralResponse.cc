#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"
#include "Tools/FQHEFiles/FQHETorusPseudopotentialTools.h" //added by ba340
#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "MainTask/FQHEOnTorusMainTask.h" //added by ba340

#include "Hamiltonian/ParticleOnTorusGenericHamiltonian.h" //added by ba340

#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "Operator/ParticleOnSphereDensityOperator.h"

#include "Operator/ParticleOnTorusDensityOperator.h" //added by ba340

#include "LanczosAlgorithm/LanczosManager.h" //added by ba340
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h" //added by ba340
#include "LanczosAlgorithm/BasicLanczosAlgorithm.h" // added by ba340

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusShort.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h" //added by ba340
#include "Architecture/ArchitectureOperation/MainTaskOperation.h" //added by ba340

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <cstring>
#include <cstdlib>
#include <iomanip>
#include <cassert>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;
using std::setw; //added by ba340

ParticleOnTorus* GetHilbertSpace(bool Statistics,int NbrParticles, int NbrFluxQuanta, int Momentum)
{
	ParticleOnTorus* Space;
	if ( Statistics == true )
			Space = new FermionOnTorus ( NbrParticles, NbrFluxQuanta, Momentum );
		else
		{

	#ifdef  __64_BITS__
			if ( ( NbrFluxQuanta + NbrParticles - 1 ) < 63 )
	#else
			if ( ( NbrFluxQuanta + NbrParticles - 1 ) < 31 )
	#endif
			{
				Space = new BosonOnTorusShort ( NbrParticles, NbrFluxQuanta, Momentum );
			}
			else
			{
				Space = new BosonOnTorus ( NbrParticles, NbrFluxQuanta, Momentum );
			}
			cout << " Hilbert space dimension = " << Space->GetHilbertSpaceDimension() << endl;
		}
    return Space;
}

int main ( int argc, char** argv )
{
    cout.precision ( 14 );

    // some running options and help
    OptionManager Manager ( "FQHETorusSpectralResponse" , "0.01" );
    OptionGroup* SystemGroup = new OptionGroup ( "system options" );
    OptionGroup* PlotOptionGroup = new OptionGroup ( "plot options" );
    OptionGroup* PrecalculationGroup = new OptionGroup ( "precalculation options" );
    OptionGroup* MiscGroup = new OptionGroup ( "misc options" );

    ArchitectureManager Architecture;
    LanczosManager Lanczos(false);

    Manager += SystemGroup;
    Manager += PlotOptionGroup;
    Architecture.AddOptionGroup ( &Manager );
    Lanczos.AddOptionGroup(&Manager);
    Manager += PrecalculationGroup;
    Manager += MiscGroup;

    ( *SystemGroup ) += new SingleStringOption ( '\0', "state", "name of the vector file describing the state whose density has to be plotted" );
    (*SystemGroup) += new SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
    (*SystemGroup) += new SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
    
    (*SystemGroup) += new SingleIntegerOption ('\n', "sr-save-interval", "number of Lanczos iterations after which the spectral response is printed, -1 for final only (default = -1)",-1);
    (*SystemGroup) += new SingleDoubleOption ('\n', "sr-omega-min", "spectral response omega min (default = 0.0)",0.0);
    (*SystemGroup) += new SingleDoubleOption ('\n', "sr-omega-max", "spectral response omega max (default = 10.0)",10.0);
    (*SystemGroup) += new SingleDoubleOption ('\n', "sr-epsilon", "spectral response epsilon (default = 1E-2)",1E-2);
    (*SystemGroup) += new SingleDoubleOption ('\n', "sr-omega-interval", "spectral response omega step size (default = 1E-2)",1E-2);
    (*SystemGroup) += new SingleDoubleOption ('\n', "sr-spectral-resolution", "spectral response omega step size (default = 1E-2)",1E-2);
    
    ( *PlotOptionGroup ) += new SingleStringOption ( '\n', "output", "output file ame (default output name replace the .vec extension of the input file with .rho or .rhorho)", 0 );
    
    (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 1024);

    ( *MiscGroup ) += new BooleanOption ( 'h', "help", "display this help" );
    
    if ( Manager.ProceedOptions ( argv, argc, cout ) == false )
    {
        cout << "see man page for option syntax or type FCICorrelation -h" << endl;
        return -1;
    }
    if ( Manager.GetBoolean ( "help" ) == true )
    {
        Manager.DisplayHelp ( cout );
        return 0;
    }
    if ( Manager.GetString ( "state" ) == 0 )
    {
        cout << "FQHETorusSpectralResponse requires an input state" << endl;
        return -1;
    }
    if ( IsFile ( Manager.GetString ( "state" ) ) == false )
    {
        cout << "can't find vector file " << Manager.GetString ( "state" ) << endl;
        return -1;
    }
     if ( Manager.GetDouble ( "sr-omega-min" ) > Manager.GetDouble ( "sr-omega-max" ) )  
    {
        cout << "incorrect range" << endl;
        return -1;
    }


    int NbrParticles = 0;
    int NbrFluxQuanta = 0;
    int Momentum = 0;
    double Ratio = 0;
    bool Statistics = false;
    
    long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
    if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
    
    if (FQHEOnTorusFindSystemInfoFromVectorFileName_SpectralResponse(Manager.GetString("state"), NbrParticles, NbrFluxQuanta, Momentum, Ratio, Statistics)==false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("state") << endl;
      return -1;
    }
   
    cout << setw ( 20 ) << std::left << "NbrParticles" << setw ( 20 ) << std::left << NbrParticles << endl;
    cout << setw ( 20 ) << std::left << "NbrFluxQuanta" << setw ( 20 ) << std::left << NbrFluxQuanta << endl;
    cout << setw ( 20 ) << std::left << "Momentum" << setw ( 20 ) << std::left << Momentum << endl;
    cout << setw ( 20 ) << std::left << "Ratio" << setw ( 20 ) << std::left << Ratio << endl;
    cout << setw ( 20 ) << std::left << "Statistics" << setw ( 20 ) << std::left << Statistics << endl;
    
    double* PseudoPotentials;
	  int NbrPseudoPotentials = 0;
	  if (Manager.GetString("interaction-file") == 0)
		{
		  cout << "an interaction file has to be provided" << endl;
		  return -1;
		}
	  else
		{
		  if (FQHETorusGetPseudopotentials(Manager.GetString("interaction-file"), NbrPseudoPotentials, PseudoPotentials) == false)
		return -1;
		}

	  char* OutputNamePrefix = new char [1024];
	  sprintf (OutputNamePrefix, "fermions_torus_spec_resp_kysym_%s_n_%d_2s_%d_ratio_%f", Manager.GetString("interaction-name"), NbrParticles, NbrFluxQuanta, Ratio);
    
    RealVector* RealState = new RealVector();
    
    if ( RealState->ReadVector ( Manager.GetString ( "state" ) ) == false )
    {
        cout << "can't open vector file " << Manager.GetString ( "state" ) << endl;
        return -1;
    }
    ComplexVector *ComplexState = new ComplexVector(*RealState);
    
    char OutputFileName[1024];
    sprintf(OutputFileName,"%s.dat", OutputNamePrefix);
    char EigenvectorName[1024];
    
    ParticleOnTorus* Space = GetHilbertSpace(Statistics, NbrParticles, NbrFluxQuanta, Momentum);
    //check dimension of space matches (cassert)
    Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
    
    bool FirstRun=true;
    // qx and qy are the Fourier modes of the density operator
    for (int qy=0;qy<NbrFluxQuanta;++qy)
      {

	ParticleOnTorus* TargetSpace = GetHilbertSpace(Statistics, NbrParticles, NbrFluxQuanta, (Momentum+qy)%NbrFluxQuanta);
	Space->SetTargetSpace(TargetSpace);
	ComplexVector* TargetVector = new ComplexVector(TargetSpace->GetHilbertSpaceDimension(),true);
	ComplexVector* TmpTargetVector = new ComplexVector(TargetSpace->GetHilbertSpaceDimension());
	for (int qx=0;qx<NbrFluxQuanta;++qx)
	{
	  //ky labels momentum eigenstates on the torus
	  for (int ky=0;ky<NbrFluxQuanta;++ky)
	  {
	    ParticleOnTorusDensityOperator Operator (Space,(ky+qy)%NbrFluxQuanta,ky,qx,Ratio);
	    VectorOperatorMultiplyOperation Operation(&Operator,ComplexState,TmpTargetVector);
	    Operation.ApplyOperation(Architecture.GetArchitecture());
	    (*TargetVector) += (*TmpTargetVector);
	  }
	}
	delete TmpTargetVector; //remember to delete these pointers
	sprintf(EigenvectorName,"%s_qy_%d", OutputNamePrefix, qy);
	
	//create hamiltonian
	
	AbstractQHEHamiltonian* Hamiltonian = new ParticleOnTorusGenericHamiltonian (TargetSpace, NbrParticles, NbrFluxQuanta, Ratio, NbrPseudoPotentials, PseudoPotentials, Architecture.GetArchitecture(), /*1024*/ Memory);
	double Shift = -10.0;	
	Hamiltonian->ShiftHamiltonian(Shift);
	
	//main task
	cout <<  "Manager at " <<  &Manager <<  endl;
	FQHEOnTorusMainTask Task(&Manager, Space, &Lanczos, Hamiltonian, Momentum, Shift, OutputFileName, FirstRun, EigenvectorName,  qy,  TargetVector);
	MainTaskOperation TaskOperation (&Task);
	TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	
	if (FirstRun==true)
	    FirstRun = false;
	    
	delete Hamiltonian;
	delete TargetSpace;

      }
      
      delete RealState;
      delete Space;
      
   

    return 0;
}
