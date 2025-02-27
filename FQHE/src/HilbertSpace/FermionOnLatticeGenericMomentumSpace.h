////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//  class of bosons on linearly indexed system for a system size such that    //
//               NbrStates  < 63 or 31 (64 bits or 32bits systems)            //
//                                                                            //
//                  states indexed by y-momentum and x-position               //
//                    (Landau gauge, A(r) = \alpha x \vec e_y)                //
//                                                                            //
//                        last modification : 11/02/2008                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef FERMIONONLATTICEGENERICMOMENTUMSPACE_H
#define FERMIONONLATTICEGENERICMOMENTUMSPACE_H

#include "config.h"
#include "MathTools/Complex.h"
#include "HilbertSpace/ParticleOnLattice.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "GeneralTools/GarbageFlag.h"

class FermionOnLatticeGenericMomentumSpace :  public ParticleOnLattice
{

 protected:

  // number of Fermions
  int NbrFermions;

  // number of States
  int NbrStates;

  // number of K-points
  int NbrKPoints;

  // length in x-direction
  int Nx;
  // length in y-direction
  int Ny;
  // total momenta in both directions
  int Kx;
  int Ky;
  // number of bands
  int NbrBands;


  // array describing each state
  unsigned long* StateDescription;
  // array giving the highest bit reached for a bososn in a given state
  int* StateHighestBit;

  // temporary state used when applying ProdA operator
  unsigned long ProdATemporaryState;
  // Lz maximum value associated to temporary state used when applying ProdA operator
  int ProdAHighestBit;  

  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table in a given HighestBit sector
  unsigned long LookUpTableMemorySize;
  // shift used in each lzmax sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used lzmax value of the state an the second 
  int** LookUpTable;

  // a table containing ranging from 0 to 2^MaximumSignLookUp - 1
  double* SignLookUpTable;
  // a table containing the mask on the bits to keep for each shift that is requested by sign evaluation
  unsigned long* SignLookUpTableMask;
  // number to evalute size of SignLookUpTable
  int MaximumSignLookUp;


  // pointer to the target space when an index is require after applying basic operation
  FermionOnLatticeGenericMomentumSpace* TargetSpace;


 public:

  enum 
    {
      NoSymmetry = 0x0,
    };

  // default constructor
  FermionOnLatticeGenericMomentumSpace();

  // basic constructor -> yields a square lattice in Landau gauge
  // 
  // nbrFermions = number of bosons
  // nx = number of momenta in simulation cell in x-direction
  // ny = number of momenta in simulation cell in y-direction
  // kx = total momentum in x-direction
  // ky = total momentum in y-direction
  // nBands = number of single particle bands
  // memory = memory that can be allocated for precalculations
  // verbose = flag indicating if any output is wanted
  FermionOnLatticeGenericMomentumSpace (int nbrFermions, int nx, int ny, int kx, int ky, int nBands, unsigned long memory = 10000000, bool verbose = true);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  FermionOnLatticeGenericMomentumSpace(const FermionOnLatticeGenericMomentumSpace& bosons);

  // virtual destructor
  //
  virtual ~FermionOnLatticeGenericMomentumSpace ();

  // assignment (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeGenericMomentumSpace& operator = (const FermionOnLatticeGenericMomentumSpace& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic();

  // get the quantization axis 
  //
  // return value = particle statistic
  virtual char GetLandauGaugeAxis();

  // get the number of sites
  //
  // return value = number of sites
  virtual int GetNbrSites();

  // get the number of sublattices
  //
  // return value = number of sublattices
  virtual int GetNbrSublattices(){return this->NbrBands;}


  // get information about any additional symmetry of the Hilbert space
  //
  // return value = symmetry id
  virtual int GetHilbertSpaceAdditionalSymmetry();

  // check whether HilbertSpace impaements ordering of operators
  //
  virtual bool HaveOrder ();
  
  // check whether a given operator \prod c^\dagger_m \prod c_n increases or decreases the index of a state
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // return value = 1, if created state is of higher value, 0 if equal, and -1 if lesser value
  virtual int CheckOrder (int* m, int* n, int nbrIndices);

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  virtual List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  virtual AbstractQuantumNumber* GetQuantumNumber (int index);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  virtual AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
						 SubspaceSpaceConverter& converter);

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  virtual void SetTargetSpace(ParticleOnLattice* targetSpace);

  // return Hilbert space dimension of the target space
  //
  // return value = Hilbert space dimension
  virtual int GetTargetHilbertSpaceDimension();

  // it is possible to change the flux through the simulation cell, keep previous periodic boundary conditions
  // Attention: this does require the Hamiltonian to be recalculated!!
  // nbrFluxQuanta = number of quanta of flux piercing the simulation cell
  virtual void SetNbrFluxQuanta(int nbrFluxQuanta);

  // change flux through cell and periodic boundary conditions
  // Attention: this does require the Hamiltonian to be recalculated!!
  // nbrFluxQuanta = number of quanta of flux piercing the simulation cell
  // solenoidX = new solenoid flux through torus in x-direction
  // solenoidY = new solenoid flux through torus in y-direction
  virtual void SetNbrFluxQuanta(int nbrFluxQuanta, double solenoidX, double solenoidY);

  // request solenoid fluxes
  // solenoidX = new solenoid flux through torus in x-direction
  // solenoidY = new solenoid flux through torus in y-direction
  //
  void GetSolenoidFluxes(double &solenoidX, double &solenoidY);

  // obtain the current setting of the flux piercing this lattice
  virtual int GetNbrFluxQuanta();

  // apply creation operator to a word, using the conventions
  // for state-coding and quantum numbers of this space
  // state = word to be acted upon
  // q = quantum number of boson to be added
  // coefficient = reference on the double where the multiplicative factor has to be stored
  virtual unsigned long Ad (unsigned long state, int q, double& coefficient);

  // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value = index of the destination state 
  virtual int AdAdAA (int index, int m1, int m2, int n1, int n2, double &coefficient);

  // apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AA (int index, int n1, int n2);

  // apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on a multiplicative, trivially one here (not changed)
  // return value = index of the destination state 
  virtual int AdAd (int m1, int m2, double& coefficient);

  // apply a^+_m a_m operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = the destination state, if non-zero exectation, otherwise: HilbertSpaceDimension
  virtual double AdA (int index, int m);

  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the multiplicative factor (result always 1.0 here)
  // return value = index of the destination state 
  virtual int AdA (int index, int m, int n, double& coefficient);

 
  // apply \sum q U_q a^+_q a_q ( a^+_q a_q - 1 )
  // index = index of the state on which the operator has to be applied
  // nbrInteraction = number of q-values in sum, if equals NbrStates, ordered sequence 0,...,NbrStates-1 assumed
  // qValues = array of quantum numbers where an interaction is present
  // interactionPerQ = coefficient U_q of the interaction
  //
  virtual double AdAdAADiagonal(int index, int nbrInteraction, double *interactionPerQ, int *qValues);

  // code set of quantum numbers posx, posy into a single integer
  // posx = position along x-direction
  // posy = position along y-direction
  // sublattice = sublattice index
  // translationPhase = returns phase occurred from translating the
  //                    site to the fundamental region [0,Lx-1] x [0,Ly-1]
  virtual int EncodeQuantumNumber(int posx, int posy, int sublattice, Complex &Phase);

  // code set of quantum numbers posx, posy into a single integer
  // posx = position along x-direction
  // posy = position along y-direction
  // sublattice = sublattice index
  // translationPhase = returns phase occurred from translating the
  //                    site to the fundamental region [0,Lx-1] x [0,Ly-1]
  virtual int FastEncodeQuantumNumber(int posx, int posy, int sublattice);

  // decode a single encoded quantum number q to the set of quantum numbers posx, posy, sublattice
  // posx = position along x-direction
  // posy = position along y-direction
  // sublattice = sublattice index
  virtual void DecodeQuantumNumber(int q, int &posx, int &posy, int &sublattice);

  // obtain a list of quantum numbers in state
  // index = index of many-body state to be considered
  // quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles
  // normalization = indicating the multiplicity of the state for bosonic spaces
  virtual void ListQuantumNumbers(int index, int *quantumNumbers, double &normalization);

  // obtain a list of quantum numbers in state
  // index = index of many-body state to be considered
  // quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles
  virtual void ListQuantumNumbers(int index, int *quantumNumbers);

  // translate a state by a multiple of the lattice vectors
  // shiftX = length of translation in x-direction
  // shiftY = length of translation in y-direction
  // translationPhase = returns phase inccurred by translation
  // return value = index of translated state
  virtual int TranslateState(int index, int shiftX, int shiftY, Complex &translationPhase);

  // find whether there is a translation vector from state i to state f
  // i = index of initial state
  // f = index of final state
  // shiftX = length of translation in x-direction
  // shiftY = length of translation in y-direction
  // return value = final state can be reached by translation
  virtual bool IsTranslation(int i, int f, int &shiftX, int &shiftY);


  // evaluate wave function in real space using a given basis
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis);

  // evaluate wave function in real space using a given basis, using time coherence
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // nextCoordinates = index of the coordinate that will be changed during the next time iteration
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
							 AbstractFunctionBasis& basis, int nextCoordinates);

  // evaluate wave function in real space using a given basis and only for agiven range of components
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
					int firstComponent, int nbrComponent);                                
  
  // evaluate wave function in real space using a given basis and only for a given range of components, using time coherence
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // nextCoordinates = index of the coordinate that will be changed during the next time iteration
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
							 AbstractFunctionBasis& basis, 
							 int nextCoordinates, int firstComponent, int nbrComponent);

  // initialize evaluation of wave function in real space using a given basis and only for a given range of components and
  //
  // timeCoherence = true if time coherence has to be used
  virtual void InitializeWaveFunctionEvaluation (bool timeCoherence = false);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // carefully test whether state is in Hilbert-space and find corresponding state index
  //
  // stateDescription = unsigned integer describing the state
  // highestBit = maximum nonzero bit reached by a particle in the state (can be given negative, if not known)
  // return value = corresponding index, or dimension of space, if not found
  virtual int CarefulFindStateIndex(unsigned long stateDescription, int highestBit);

 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int lzmax);

  
  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // currentQ = current largest combined quantum number
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int nbrFermions, int currentQ, int currentTotalKx, int currentTotalKy);


  // generate states with kx and ky symmetries
  //
  // nbrFermions = number of fermions
  // currentQ = current consolidated quantum number
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // pos = position where to start filling states
  // return value = Hilbert space dimension
  
  long GenerateStates(int nbrFermions, int currentQ, int currentTotalKx, int currentTotalKy, long pos);
    

  // generate look-up table associated to current Hilbert space
  // 
  // memory = memory size that can be allocated for the look-up table
  
  void GenerateLookUpTable(unsigned long memory);

  
};


// decode a single encoded quantum number q to the set of quantum numbers kx, ky
// kx = momentum along x-direction
// ky = momentum along y-direction
inline void FermionOnLatticeGenericMomentumSpace::DecodeQuantumNumber(int q, int &kx, int &ky, int &band)
{
  band=q/(this->NbrKPoints);
  q=q%this->NbrKPoints;
  ky=q/this->Nx;
  kx=q%this->Nx;
}


// code set of quantum numbers posx, posy into a single integer
// kx = momentum along x-direction
// ky = momentum along y-direction
// band = band index
// 
inline int FermionOnLatticeGenericMomentumSpace::EncodeQuantumNumber(int kx, int ky, int band, Complex &translationPhase)
{
  //cout << "Encoding " << kx<<", "<<ky<<": ";
  while (kx<0)
    kx+=this->Nx;
  while (kx>=this->Nx)
    kx-=this->Nx;
  while (ky<0)
    ky+=this->Ny;
  while (ky>=this->Ny)
    ky-=this->Ny;
  int rst = kx + this->Nx*ky;
  rst+=band*this->NbrKPoints;
  translationPhase=0.0;
  return rst;
}

// code set of quantum numbers posx, posy into a single integer - same as above, but shortened interface
// code set of quantum numbers posx, posy into a single integer
// kx = momentum along x-direction
// ky = momentum along y-direction
// band = band index
inline int FermionOnLatticeGenericMomentumSpace::FastEncodeQuantumNumber(int kx, int ky, int band)
{
  int rst = kx + this->Nx*ky;
  rst+=band*this->NbrKPoints;
  return rst;
}

#endif


