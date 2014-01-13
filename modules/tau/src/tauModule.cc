/**************************************************************************
 * BASF2 (Belle Analysis Framework 2)                                     *
 * Copyright(C) 2013 - Belle II Collaboration                             *
 *                                                                        *
 * Author: The Belle II Collaboration                                     *
 * Contributors: rsobie                                                   *
 *                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/
#include <iomanip>

#include </home/rsobie/release/mypkg/modules/tau/include/tauModule.h>
#include </home/rsobie/release/mypkg/modules/tau/include/PdtLund.h>

#include <framework/datastore/StoreObjPtr.h>
#include <framework/datastore/StoreArray.h>  
#include <framework/dataobjects/EventMetaData.h>

#include <generators/dataobjects/MCParticle.h>  

#include <tracking/dataobjects/TrackFitResult.h>

#include <TVector3.h>

#include </home/rsobie/release/mypkg/modules/tau/src/TauDecayMode.h>

using namespace Belle2;

//-----------------------------------------------------------------
//                 Register the Module
//-----------------------------------------------------------------
REG_MODULE(tau)

//-----------------------------------------------------------------
//                 Implementation
//-----------------------------------------------------------------

tauModule::tauModule() : Module()
{
  // Set module properties
  setDescription("testing out my new tau");

  // Parameter definitions

}

tauModule::~tauModule()
{
}

void tauModule::initialize()
{
   std::cout << "entered tau initialize" << std::endl;
}

void tauModule::beginRun()
{
}

void tauModule::event()
{

  StoreObjPtr<EventMetaData> eventmetadata;
  if(!eventmetadata) {
    B2INFO("an object called '" << eventmetadata.getName() << "' does not exist in the data store.");
    } else {
  //object exists, you can now access its data
   //  B2INFO("we're currently in event " << eventmetadata->getEvent() << "!");

   //std::cout << "evnt info: " << eventmetadata->getEvent() <<" "
//	<< eventmetadata->getRun() <<" "
//	<< eventmetadata->getExperiment() << std::endl;
  }



  // ---------------------------------------------------
  // charged tracks
  // ---------------------------------------------------
  //printTracks();

  // ---------------------------------------------------
  // Tau MC generator information
  // ---------------------------------------------------
  // (not needed anymore?) StoreArray<MCParticle> mcList;
  TauDecayMode myMode;
  int tauPlus  = myMode.getTauMCID( (int)1);
  int tauMinus = myMode.getTauMCID( (int)-1);

  int nProngsPlus  = myMode.getNumberProngs( tauPlus);
  int nProngsMinus = myMode.getNumberProngs( tauMinus);
  // std::cout << "nprongs = " << nProngsPlus <<" "<< nProngsMinus << std::endl;

  //myMode.printTauMCInfo( tauPlus );
  //myMode.printTauMCInfo( tauMinus );

  int tmp1 = myMode.getDecayMode ( tauPlus );
  if(tmp1 >= 0 ) myMode.printFinalState( tauPlus );

  int tmp2 = myMode.getDecayMode ( tauMinus );
  if(tmp2 >= 0 ) myMode.printFinalState( tauMinus );

  std::cout << "tau decay modes : " << tmp1 <<" "<< tmp2 << std::endl;


  // ---------------------------------------------------
  // find MC tau leptons
  // ---------------------------------------------------
  //printDecay();

}

void tauModule::endRun()
{
}

void tauModule::terminate()
{
}

void tauModule::printTracks() 
{
  TVector3 tkmom;
  StoreArray<TrackFitResult> trackList;

  // ---------------------------------------------------
  // track information
  // ---------------------------------------------------
  std::cout << std::endl;
  std::cout << "*** Charged track properties *** " << std::endl;
  std::cout << "Number of tracks = " << trackList.getEntries() << std::endl;
  std::cout << " N"
            << "     P"
            << "  cosT"
            << "   phi"
            << "  Q"
            << std::endl;

   for (int i = 0; i < trackList.getEntries(); i++) {
        tkmom = trackList[i]->getMomentum();
        std::cout << std::setiosflags(std::ios::fixed);
        std::cout << std::setw(2) << i 
                  << std::setw(6) << std::setprecision(2) << tkmom.Mag() 
                  << std::setw(6) << std::setprecision(2) << tkmom.CosTheta() 
                  << std::setw(6) << std::setprecision(2) << tkmom.Phi()
                  << std::setw(3) << std::setprecision(0) << trackList[i]->getCharge()
                  << std::endl;
   }
}

// -------------------------------------------------------------------------
// printDecay - print tau MC information
// December 2013
// -------------------------------------------------------------------------
void tauModule::printDecay()
{
  // ---------------------------------------------
  // find MC tau leptons
  // ---------------------------------------------
  short TauMinus = -1, TauPlus = -1;
  StoreArray<MCParticle> mcList;
  for (int i = 0; i < mcList.getEntries(); i++) {
   if(TauMinus == -1) { if(mcList[i]->getPDG() == PdtLund::tau_minus)  TauMinus = i; }
   if(TauPlus  == -1) { if(mcList[i]->getPDG() == PdtLund::tau_plus )  TauPlus  = i; }
  }

  // ---------------------------------------------
  // print out tau properties
  // ---------------------------------------------
  std::cout << std::endl;
  std::cout << "*** Tau lepton properties ***" << std::endl;
  std::cout << "TauQ"
            << "     P"
            << "  cosT"
            << "   phi"
            << std::endl;

  TVector3 taumv = mcList[TauMinus]->getMomentum();
  for(int i=0; i<2; i++) {
     int myTau = TauPlus;
     if( i==1 ) myTau = TauMinus;
     TVector3 tauV = mcList[myTau]->getMomentum();
     if(mcList[myTau]->getCharge() ==  1 ) std::cout << "Tau+";
     if(mcList[myTau]->getCharge() == -1 ) std::cout << "Tau-";
     std::cout << std::setw(6) << std::setprecision(2) << tauV.Mag()
               << std::setw(6) << std::setprecision(2) << tauV.CosTheta()
               << std::setw(6) << std::setprecision(2) << tauV.Phi()
               << std::endl;
  }

  // ---------------------------------------------
  // identify tau decay mode
  // ---------------------------------------------
  //    const std::vector<Belle2::MCParticle*> getDaughters() const;

  std::vector<MCParticle*> tauDaugVector = mcList[TauPlus]->getDaughters();

  std::cout << "n daug = " << tauDaugVector.size() << std::endl;

  std::cout << "Tau daughters : ";
  for(unsigned int i=0; i<tauDaugVector.size(); i++) {
     std::cout << tauDaugVector[i]->getPDG() <<" ";
  }
  std::cout << std::endl;
            

  // ---------------------------------------------
  // number of prongs
  // ---------------------------------------------
  unsigned int nProng = 0;
  for(unsigned int i=0; i<tauDaugVector.size(); i++) {
    if( tauDaugVector[i]->getCharge() != 0 ) nProng++;
  }
  std::cout << "N-prongs (tau+ only) = " << nProng << std::endl;


  // ---------------------------------------------
  // count the particle types
  // ---------------------------------------------
  int _nnue(0), _nnum(0), _nnut(0);
  int _nele(0), _nmu(0), _npi(0), _npi0(0), _nrho(0), _na1(0);
  int _nK(0), _nKL(0), _nKS(0);
  int _nw(0);
  for(unsigned int i=0; i<tauDaugVector.size(); i++) {
     if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::nu_e    ) ) _nnue++;
     if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::nu_mu   ) ) _nnum++;
     if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::nu_tau  ) ) _nnut++;

     if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::e_plus  ) ) _nele++;
     if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::mu_plus ) ) _nmu++;
     if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::pi_plus ) ) _npi++;
     if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::pi0     ) ) _npi0++;
     if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::rho_plus) ) _nrho++;
     if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::a_1_plus) ) _na1++;

     if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::K_plus  ) ) _nK++;
     if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::K_L0    ) ) _nKL++;
     if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::K_S0    ) ) _nKS++;

     if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::W_plus  ) ) _nw++;
   }
 
  // double check sum of all counters with number of expected particles

   std::cout << "Number of decay particles: nu(e,m,t)=" 
             << _nnue <<" "<< _nnum <<" "<< _nnut << "  W = " << _nw 
             << std::endl
             << "(e,mu,pi,pi0,rho,a1,K,KL,KS) : " 
             << _nele <<" "<< _nmu <<" "<< _npi <<" "<< _npi0 <<" "<< _nrho <<" "<< _na1 <<" "
             << _nK <<" "<< _nKL <<" "<< _nKS 
             << std::endl;

}


