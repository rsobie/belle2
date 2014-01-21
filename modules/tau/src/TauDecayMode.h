// ----------------------------------------------------------------
//
// Tau MC properties
//
// December 2013
// ----------------------------------------------------------------
// Modes: (drop neutrinos)

// nFirst = 1 and not-W decay:
// 1 		e- 
// 2		mu-
// 3		pi-
// 4 		K-
// 5            rho- (assuming all decays to pi+ pi-)
// 11           a1- -> pi- pi- pi+
// 12           a1- -> pi- pi0 pi0 
// 13           a1 other decays
// 51           K*- -> pi0 K-
// 52           K*- -> pi- KL
// 53           K*- -> pi- KS

// ----------------------------------------------------------------


//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class TauDecayMode {

public:

  TauDecayMode();
  virtual ~TauDecayMode();
  int  getTauMCID       ( int );
  void printTauMCInfo   ( int );
  void printFinalState  ( int );
  void countMesons      ( int, int );
  int countResonanceParticles ( int, int );
  int  getDecayMode     ( int );
  int  getNumberProngs  ( int );
 
  int countParticles(   std::vector<Belle2::MCParticle*>  );

private:
  int _tauPtr;
  int _nProng;

  bool _wDecay;
  int  _ptr;


  bool _Resonance, _Strange;
  int _nnue, _nnum, _nnut;
  int _ng, _nele, _nmu,  _npi0, _nrho, _na1, _nw;
  int _nomega, _neta, _nf1;
  int _nK0, _nKL, _nKS, _nKstar0, _nKstarp;
  int _npi, _npip, _npim;
  int _nK, _nKp, _nKm;

  int _nrespi, _nrespi0, _nresgam, _nreseta, _nresetap, _nresomega;
  int _nresK, _nresKS, _nresKL, _nresK0;

  int _pdgLabel[100], _pdgCount[100];

  int _mode;

};


// constructor
TauDecayMode::TauDecayMode() {
  _tauPtr = 0;
  _nProng = 0;
}

// destructor
TauDecayMode::~TauDecayMode() {}

// ----------------------------------------------------------------
// get tau mcID pointer
// ----------------------------------------------------------------
int TauDecayMode::getTauMCID( int tauCharge ) {
 
  // ---------------------------------------------
  // find MC tau leptons
  // ---------------------------------------------
  _tauPtr = 0;
  Belle2::StoreArray<Belle2::MCParticle> mcList;
  for (int i = 0; i < mcList.getEntries(); i++) {
     if( tauCharge == 1 ) {
        if(mcList[i]->getPDG() == PdtLund::tau_plus )  { _tauPtr = i; break;}
     } else {
        if(mcList[i]->getPDG() == PdtLund::tau_minus ) { _tauPtr = i; break;}
     }
  }

  return _tauPtr;
}

// ----------------------------------------------------------------
// ----------------------------------------------------------------
void TauDecayMode::printTauMCInfo( int tauPtr ) 
{
  Belle2::StoreArray<Belle2::MCParticle> mcList;

  TVector3 tauVec  = mcList[tauPtr]->getMomentum();
  std::cout << std::endl << "Tau MC Info (P,cosTheta,Phi,Q): ";
  std::cout << std::setw(6) << std::setprecision(2) << tauVec.Mag()
            << std::setw(6) << std::setprecision(2) << tauVec.CosTheta()
            << std::setw(6) << std::setprecision(2) << tauVec.Phi()
            << std::setw(3) << std::setprecision(0) << mcList[tauPtr]->getCharge()
            << std::setw(3) << std::setprecision(0) << _nProng
            << std::endl;
}

// ----------------------------------------------------------------
// get number of charged tracks (prongs) in this decay 
// ----------------------------------------------------------------
int TauDecayMode::getNumberProngs( int tauPtr ) 
{
  _nProng = 0;
  Belle2::StoreArray<Belle2::MCParticle> mcList;
  std::vector<Belle2::MCParticle*> tauDaugVector = mcList[tauPtr]->getDaughters();
  for(unsigned int i=0; i<tauDaugVector.size(); i++) {
    if( tauDaugVector[i]->getCharge() != 0 ) _nProng++;
  }

// this does not work for resonances (W and a1)

  return _nProng;
}


// ----------------------------------------------------------------
// print final state 
// ----------------------------------------------------------------
void TauDecayMode::printFinalState( int tauPtr )  
{
  // ----------------------------------------
  // check particle pointers
  // ----------------------------------------
  if( tauPtr == 0 )        return;
  if( _wDecay && _ptr==0 ) return;

  // ----------------------------------------
  // get list of daughers from tau or W-boson
  // ----------------------------------------
  Belle2::StoreArray<Belle2::MCParticle> mcList;
  std::vector<Belle2::MCParticle*> tauDV = mcList[tauPtr]->getDaughters();
  std::vector<Belle2::MCParticle*> tauDaugVector;
  if( _wDecay) {
      tauDaugVector = tauDV[ _ptr ]->getDaughters();
  } else {
      tauDaugVector = mcList[tauPtr]->getDaughters();
  }

  // ----------------------------------------
  // tau 4-vector 
  // ----------------------------------------
   if( mcList[tauPtr]->getCharge() == 1 ) {
      std::cout << "Tau+ (p,CT,Phi) = "; 
   } else {
      std::cout << "Tau- (p,CT,Phi) = ";
   }
   TVector3 tauVec  = mcList[tauPtr]->getMomentum();
   std::cout << std::setw(6) << std::setprecision(2) << tauVec.Mag()
             << std::setw(6) << std::setprecision(2) << tauVec.CosTheta()
             << std::setw(6) << std::setprecision(2) << tauVec.Phi();

   std::cout << "  Daughters : "; 
   for(unsigned int i=0; i<tauDaugVector.size(); i++) {
       std::cout << std::setw(8) << tauDaugVector[i]->getPDG();
   }
   std::cout << std::endl;

  // ----------------------------------------
  // particles in primary final state
  // ----------------------------------------
   std::cout << "g N(emt) L(em) pi(*+-) K(*+-) M(pz rh om et f1 a1) K0(0LS *0-) " << std::endl;;
   std::cout << std::setw(3) << std::setprecision(0);
   std::cout << std::setw(1) << _ng
             << std::setw(4) << _nnue   << _nnum << _nnut
             << std::setw(5) << _nele   << _nmu
             << std::setw(6) << _npi    << _npip << _npim
             << std::setw(5) << _nK     << _nKp << _nKm
             << std::setw(6) << _npi0 << std::setw(3) << _nrho << std::setw(3) << _nomega 
             << std::setw(3) << _neta << std::setw(3) << _nf1  << std::setw(3) << _na1
             << std::setw(6) << _nK0  << _nKL  << _nKS 
             << std::setw(3) << _nKstar0 << _nKstarp
             << std::endl;

  // -----------------------------------------------
  // search for resonances (omega, eta, f1, a1, K0s)
  // -----------------------------------------------
  if( _neta!=1 && _nomega!=0 ) return;

  int rPtr = 0;
  for(unsigned int i=0; i<tauDaugVector.size(); i++) {
     if( tauDaugVector[i]->getPDG() == PdtLund::eta ||
         tauDaugVector[i]->getPDG() == PdtLund::omega ||
         tauDaugVector[i]->getPDG() == PdtLund::f_1 ||
         abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::a_1_plus))  { rPtr = i; break; }
  }
  if( rPtr == 0 ) return;
  std::vector<Belle2::MCParticle*> daugVector = tauDaugVector[rPtr]->getDaughters();

  std::cout << "===>RJS test: " << tauDaugVector[rPtr]->getPDG() << "  Daughters : ";
   for(unsigned int i=0; i<daugVector.size(); i++) {
       std::cout << std::setw(6) << daugVector[i]->getPDG();
   }
   std::cout << std::endl;

}

// ----------------------------------------------------------------
// count mesons from a1 and W 
// ----------------------------------------------------------------
void TauDecayMode::countMesons( int tauPtr, int ptr ) 
{
  Belle2::StoreArray<Belle2::MCParticle> mcList;
  std::vector<Belle2::MCParticle*> tauDaugVector = mcList[tauPtr]->getDaughters();

  std::vector<Belle2::MCParticle*> daugVector = tauDaugVector[ptr]->getDaughters();

  std::cout << "===>RJS countMeson test: " << tauDaugVector[ptr]->getPDG() <<" "
            << daugVector.size() << std::endl;


  // temporary change from normal counters
  int xpi=0, xpi0=0, xK=0, xg=0;

  for(unsigned int i=0; i<daugVector.size(); i++) {
    std::cout << "===>RJS daug = " << daugVector[i]->getPDG() << std::endl;

     if(      abs(daugVector[i]->getPDG()) == abs(PdtLund::pi_plus ) ) { xpi++;  }
     else if( abs(daugVector[i]->getPDG()) == abs(PdtLund::pi0     ) ) { xpi0++; }
     else if( abs(daugVector[i]->getPDG()) == abs(PdtLund::K_plus  ) ) { xK++;   }
     else if( abs(daugVector[i]->getPDG()) == abs(PdtLund::gamma   ) ) { xg++;   }
     else { std::cout << "===>RJS: Did not find the particle " << daugVector[i]->getPDG() << std::endl; }
  }

  std::cout << std::setw(3) << std::setprecision(0);
  std::cout << std::setw(3) << "===>RJS: testing a1 (pi,pi0,K) = " << xpi << xpi0 << xK << std::endl;

}

// ----------------------------------------------------------------
// count resonance particles
// ----------------------------------------------------------------
int TauDecayMode::countResonanceParticles( int tauPtr, int ptr ) 
{
   int nResParticles = 0;

   _nrespi=0;  _nrespi0=0; _nresgam=0; 
   _nresK=0; _nresKS=0; _nresKL=0; _nresK0=0;
   _nreseta=0; _nresetap=0; _nresomega=0;
    
   Belle2::StoreArray<Belle2::MCParticle> mcList;
   std::vector<Belle2::MCParticle*> tauDaugVector = mcList[tauPtr]->getDaughters();
 
   std::vector<Belle2::MCParticle*> rDaugVector = tauDaugVector[ptr]->getDaughters();
   for(unsigned int i=0; i<rDaugVector.size(); i++) {
       if(      abs(rDaugVector[i]->getPDG()) == abs(PdtLund::pi_plus   ) ) _nrespi++;
       else if( abs(rDaugVector[i]->getPDG()) == abs(PdtLund::pi0       ) ) _nrespi0++;
       else if( abs(rDaugVector[i]->getPDG()) == abs(PdtLund::gamma     ) ) _nresgam++;
       else if( abs(rDaugVector[i]->getPDG()) == abs(PdtLund::eta       ) ) _nreseta++;
       else if( abs(rDaugVector[i]->getPDG()) == abs(PdtLund::eta_prime ) ) _nresetap++;
       else if( abs(rDaugVector[i]->getPDG()) == abs(PdtLund::omega     ) ) _nresomega++;
       else if( abs(rDaugVector[i]->getPDG()) == abs(PdtLund::K_plus    ) ) _nresK++;
       else if( abs(rDaugVector[i]->getPDG()) == abs(PdtLund::K0        ) ) _nresK0++;
       else if( abs(rDaugVector[i]->getPDG()) == abs(PdtLund::K_S0      ) ) _nresKS++;
       else if( abs(rDaugVector[i]->getPDG()) == abs(PdtLund::K_L0      ) ) _nresKL++;
       else { std::cout << "===>RJS: Unidentified resonance decay particle " 
                        << rDaugVector[i]->getPDG() << std::endl; }
       if( rDaugVector[i]->getPDG() != PdtLund::gamma ) nResParticles++;
     }

       std::cout << std::setw(3) << std::setprecision(0);
       std::cout << "Resonance: " << tauDaugVector[ptr]->getPDG()
             << " N (g) (pi pi0 eta etap omega) (K K0 KS KL) = "
             << std::setw(3) << nResParticles
             << std::setw(3) << _nresgam
             << std::setw(3) << _nrespi << _nrespi0 << _nreseta << _nresetap << _nomega
             << std::setw(3) << _nresK  << _nresK0 << _nresKS << _nresKL
             << std::endl;

   return nResParticles;
}


// ----------------------------------------------------------------
// get decay mode 
// ----------------------------------------------------------------
int TauDecayMode::getDecayMode( int tauPtr ) 
{
  _mode = 0;
  Belle2::StoreArray<Belle2::MCParticle> mcList;

  // ----------------------------------------------------------------
  // identify decays via a W-boson
  // ----------------------------------------------------------------
  _ptr = 0;
  _wDecay = false;
  std::vector<Belle2::MCParticle*> tauDV = mcList[tauPtr]->getDaughters();
  for(unsigned int i=0; i<tauDV.size(); i++) {
     if( abs(tauDV[i]->getPDG()) == abs(PdtLund::W_plus) )  { 
         _ptr = i;
         _wDecay = true;
         break;
     }     
  }

  std::vector<Belle2::MCParticle*> tauDaugVector;
  if( _wDecay) {
      tauDaugVector = tauDV[ _ptr ]->getDaughters();
  } else {
      tauDaugVector = mcList[tauPtr]->getDaughters();
  }


  // ----------------------------------------------------------------
  // count number of first decay particles (excluding neutrinos/gammas)
  // ----------------------------------------------------------------
  int nFirst = 0;
  for(unsigned int i=0; i<tauDaugVector.size(); i++) {
    if( abs(tauDaugVector[i]->getPDG()) != abs(PdtLund::gamma)  &&
        abs(tauDaugVector[i]->getPDG()) != abs(PdtLund::nu_e)   &&
        abs(tauDaugVector[i]->getPDG()) != abs(PdtLund::nu_mu)  &&
        abs(tauDaugVector[i]->getPDG()) != abs(PdtLund::nu_tau)  )  nFirst++;
  }

  // ----------------------------------------------------------------
  // count all particle-types in first decay
  // ----------------------------------------------------------------
  int resPtr = countParticles( tauDaugVector );


//  int tauDecayMode(0);
  int wPtr(0),  kstarPtr(0);
  int nP;

  // ----------------------------------------------------------------
  // decays with a single primary decay particle
  // ----------------------------------------------------------------
  if( nFirst == 1 && !_wDecay )  {
    if(       _nele == 1 ) { _mode=1; }  // e
     else if( _nmu  == 1 ) { _mode=2; }  // mu
     else if( _npi  == 1 ) { _mode=3; }  // pi
     else if( _nK   == 1 ) { _mode=4; }  // K
     else if( _nrho == 1 ) { _mode=5; }  // rho

     else if (_na1  == 1 ) {             // a1
        nP = countResonanceParticles( tauPtr, resPtr );
        if(        nP==3 && _nrespi==3 && _nrespi0==0) { _mode = 11; } // a1- -> pi- pi- pi+
         else if(  nP==3 && _nrespi==1 && _nrespi0==2) { _mode = 12; } // a1- -> pi- pi0 pi0
         else                                          { _mode = 13; } // other a1 decays
     } // end-a1

     else if (_nKstarp  == 1 ) {        // K*-
        nP = countResonanceParticles( tauPtr, resPtr );
        if(        nP==2 && _nrespi0==1 && _nresK==1)  { _mode = 51; } // K*- -> pi0 K-
         else if(  nP==2 && _nrespi==1  && _nresKL==1) { _mode = 52; } // K*- -> pi- KL
         else if(  nP==2 && _nrespi==1  && _nresKS==1) { _mode = 53; } // K*- -> pi- KS
         else                                          { _mode = 59; } // other K*- decays
     }  // end-K*-

     else { _mode = 99; }

  } // end-NF1 non-W


  if( _wDecay ) std::cout << "===>RJS: W-decay nFirst = " << nFirst 
                         << " (pi pi0) = " << _npi << _npi0 << std::endl;

  // ---------------------------------------------------------------------
  // W decays
  // ---------------------------------------------------------------------
  if( nFirst == 2 && _wDecay ) {             
     if(      _npi==1 && _npi0==1 && _ng>=1 )  { _mode = 1200; } // W -> pi- pi0 (>= 1gamma)
     else if( _npi==1 && _neta==1 )            { _mode = 1210; } // W -> pi- eta   (** none observed **)
     else if( _npi==1 && _nomega==1 )          { _mode = 1220; } // W -> pi- omega (** none observed **)
     else if( _nK==1  && _nKL==1    )          { _mode = 1230; } // W -> K- KL
     else if( _nK==1  && _nKS==1    )          { _mode = 1240; } // W -> K- KS
     else                                      { _mode = 1299; } // other 2-body W-decays
  }

  // Notes: 
  // (i)   W -> Kpipi (K has same charge as tau/W)
  // (ii)  pi- KL KL and pi- KS KS are odd decays
  // (iii) no K-K-K+ decays observed in the MC
  if( nFirst == 3 && _wDecay ) {             
    if(       _npi==1 && _npi0==1 && _neta==1   )     { _mode = 1300; } // W -> pi- pi0 eta
     else if( _npi==1 && _npi0==1 && _neta==1   )     { _mode = 1310; } // W -> pi- pi0 omega
     else if( _npi==1 && _npi0==1 && _nKL==1    )     { _mode = 1320; } // W -> pi- pi0 KL
     else if( _npi==1 && _npi0==1 && _nKS==1    )     { _mode = 1321; } // W -> pi- pi0 KS
     else if( _npi==1 && _nKL==1  && _nKS==1    )     { _mode = 1330; } // W -> pi- KL KS
     else if( _npi==1 && _nKL==2                )     { _mode = 1331; } // W -> pi- KL KL
     else if( _npi==1 && _nKS==2                )     { _mode = 1332; } // W -> pi- KS KS
     else if( _npi==1 && _nK==2                 )     { _mode = 1335; } // W -> pi- K- K+

     else if( _nK==1  && _npi==2                )     { _mode = 1340; } // W -> K-  pi+ pi-
     else if( _nK==1  && _npi0==2               )     { _mode = 1341; } // W -> K-  pi0 pi0
     else if( _nK==1  && _npi0==1 && _nKL==1    )     { _mode = 1350; } // W -> K-  pi0 KL
     else if( _nK==1  && _npi0==1 && _nKS==1    )     { _mode = 1351; } // W -> K-  pi0 KS

     else if( _nK==3                            )     { _mode = 1360; } // W -> K-  K- K+ (** none observed **)
     else                                             { _mode = 1399; } // other 3-body W-decays
  }

  if( nFirst == 4 && _wDecay ) {             
     if(      _npi==1 && _npi0==3 )    { _mode = 1400; } // W -> pi 3pi0
     else if( _npi==3 && _npi0==1 )    { _mode = 1410; } // W -> 3pi pi0
     else                              { _mode = 1499; } // other 4-body W-decays
  }

  if( nFirst == 5 && _wDecay ) {             
     if(      _npi==5  )               { _mode = 1500; } // W -> 5pi
     else if( _npi==3 && _npi0==2 )    { _mode = 1510; } // W -> 3pi 2pi0 
     else if( _npi==1 && _npi0==4 )    { _mode = 1520; } // W ->  pi 4pi0
     else                              { _mode = 1599; } // other 5-body W-decays
  }

  if( nFirst == 6 && _wDecay ) {
     if( _npi==6 && _npi0==1  )        { _mode = 1600; } // W -> 5pi pi0
     else                              { _mode = 1699; } // other 6-body W-decays
  }

  return _mode;


    if( nFirst == 1 && _wDecay ) {              // W -> quark-quark 
       countResonanceParticles( tauPtr, wPtr );
       int Strange =   (_nresK >=1   || _nresKS >=1 || _nresKL >=1);
       int Resonance = (_nreseta >=1 || _nresetap >=1);

       // pi- npi0
       if( _nrespi==1 && !Strange && !Resonance ) {
          if(       _nrespi0==2 )                 { _mode = 21; } // W -> pi 2pi0
           else if( _nrespi0==3 )                 { _mode = 22; } // W -> pi 3pi0
           else if( _nrespi0>=4 )                 { _mode = 23; } // W -> pi 4pi0
           else if( _nrespi0==1 && _nresgam>=1  ) { _mode = 24; } // W -> pi- pi0 gamma
           else                                   { _mode =-25; } // other W -> pi decays
       }

       // 3pi- npi0
       if( _nrespi==3 && !Strange && !Resonance ) {
          if(       _nrespi0==1 )                 { _mode = 31; } // W -> 3pi pi0
           else if( _nrespi0==2 )                 { _mode = 32; } // W -> 3pi 2pi0
           else if( _nrespi0>=3 )                 { _mode = 33; } // W -> 3pi 3pi0
           else                                   { _mode =-35; } // other W -> 3pi decays
       }

       // pi- eta
       if( _nrespi==1 &&_nreseta==1 && !Strange && Resonance ) {
          if(       _nrespi0==0 )                 { _mode = 110; } // W -> eta pi- 
           else if( _nrespi0==1 )                 { _mode = 120; } // W -> eta pi- pi0
           else if( _nrespi0>=2 )                 { _mode = 130; } // W -> eta pi- 2pi0
           else                                   { _mode =-140; } // other W -> pi- eta decays
       }

       // 3pi- eta
       if( _nrespi==3 &&_nreseta==1 && !Strange && Resonance ) {
          if(       _nrespi0==0 )                 { _mode = 150; } // W -> eta 3pi-
           else if( _nrespi0==1 )                 { _mode = 160; } // W -> eta 3pi- pi0
           else if( _nrespi0>=2 )                 { _mode = 170; } // W -> eta 3pi- 2pi0
           else                                   { _mode =-180; } // other W -> 3pi- eta decays
       }

       // pi- omega 
       if( _nrespi==1 &&_nresomega==1 && !Strange && Resonance ) {
          if(       _nrespi0==0 )                 { _mode = 210; } // W -> omega pi-
           else if( _nrespi0==1 )                 { _mode = 220; } // W -> omega pi- pi0
           else if( _nrespi0>=2 )                 { _mode = 230; } // W -> omega pi- 2pi0
           else                                   { _mode =-240; } // other W -> pi- omega decays
       }

       // 3pi- omega
       if( _nrespi==3 &&_nresomega==1 && !Strange && Resonance ) {
          if(       _nrespi0==0 )                 { _mode = 250; } // W -> omega 3pi-
           else if( _nrespi0==1 )                 { _mode = 260; } // W -> omega 3pi- pi0
           else if( _nrespi0>=2 )                 { _mode = 270; } // W -> omega 3pi- 2pi0
           else                                   { _mode =-280; } // other W -> 3pi- omega decays
       }

       // pi- eta'
       if( _nrespi==3 &&_nresetap==1 && !Strange && Resonance ) {
         if(        _nrespi0==0 )                 { _mode = 310; } // W -> eta' pi-
           else if( _nrespi0>=1 )                 { _mode = 320; } // W -> eta' pi- pi0
           else                                   { _mode =-330; } // other W -> eta' pi-  decays
        }

       // pi pi K
       if( (_nrespi+_nrespi0)==2 && Strange && !Resonance ) {
          if(      _nrespi==1 && _nrespi0==1 && _nresK==0 && _nresKS==1 && _nresKL==0) {_mode=410;} // pi- pi0 KS
          else if( _nrespi==1 && _nrespi0==1 && _nresK==0 && _nresKS==0 && _nresKL==1) {_mode=420;} // pi- pi0 KL
          else if( _nrespi==2 && _nrespi0==0 && _nresK==1 && _nresKS==0 && _nresKL==0) {_mode=430;} // pi- pi+ K-
          else if( _nrespi==0 && _nrespi0==2 && _nresK==1 && _nresKS==0 && _nresKL==0) {_mode=440;} // pi0 pi0 K-
          else {_mode=-450;} // *** CAREFUL WITH THIS ONE - MAY PICK UP THE WRONG DECAYS ***
          // PUT THE (pi K K) and (K K) decays in this subsection
       }

       // pi K K (mode 520, 540 or 560 observed so far)
       if( !Resonance ) {
         if(        _nrespi==1 && _nrespi0==0 && _nresK==2 && _nresKS==0 && _nresKL==0 )  {_mode=510;} // K- K+ pi- 
           else if( _nrespi==1 && _nrespi0==1 && _nresK==2 && _nresKS==0 && _nresKL==0 )  {_mode=520;} // K- K+ pi- pi0

         if(        _nrespi==0 && _nrespi0==1 && _nresK==1 && _nresKS==1 && _nresKL==0 )  {_mode=530;} // K- KS pi0
          else if ( _nrespi==0 && _nrespi0==2 && _nresK==1 && _nresKS==1 && _nresKL==0 )  {_mode=540;} // K- KS pi0 pi0
          else if ( _nrespi==0 && _nrespi0==1 && _nresK==1 && _nresKS==0 && _nresKL==1 )  {_mode=550;} // K- KL pi0 
          else if ( _nrespi==0 && _nrespi0==2 && _nresK==1 && _nresKS==0 && _nresKL==1 )  {_mode=560;} // K- KL pi0 pi0

         if(        _nrespi==0 && _nrespi0==0 && _nresK==1 && _nresKS==1 && _nresKL==0 )  {_mode=610;} // K- KS
           else if( _nrespi==0 && _nrespi0==0 && _nresK==1 && _nresKS==0 && _nresKL==1 )  {_mode=620;} // K- KL

         if(        _nrespi==1 && _nrespi0==0 && _nresK==0 && _nresKS==1 && _nresKL==1 )  {_mode=710;} // pi- KS KL
          else if ( _nrespi==1 && _nrespi0==0 && _nresK==0 && _nresKS==2 && _nresKL==0 )  {_mode=720;} // pi- KS KS
          else if ( _nrespi==1 && _nrespi0==0 && _nresK==0 && _nresKS==0 && _nresKL==2 )  {_mode=730;} // pi- KL KL
       }

       if( !Resonance && !Strange ) {
        if(        _nrespi==5 && _nrespi0==0 ) {_mode=810;} // 3pi- 2pi+
         else if(  _nrespi==5 && _nrespi0==1 ) {_mode=820;} // 3pi- 2pi+ pi0
       }


     } // end-W
   
   // K*
   else if (_nKstarp  == 1 ) {             // kstar
        countResonanceParticles( tauPtr, kstarPtr ); 
        if(      _nrespi==1 && _nrespi0==0 && _nresKS==1 && _nresKL==0 && _nresK==0) {_mode=910;} // pi- KS
        else if( _nrespi==1 && _nrespi0==0 && _nresKS==0 && _nresKL==1 && _nresK==0) {_mode=920;} // pi- KL
        else if( _nrespi==0 && _nrespi0==1 && _nresKS==0 && _nresKL==0 && _nresK==1) {_mode=930;} // pi0 K-
        else {_mode=-940;} // other kstar decays
   }


   if( _mode ==0 ) _mode = 999;
  

  std::cout << "===>RJS: temporary exit with mode = " << _mode << std::endl;
  return _mode;

  // ---------------------------------------------------
  // 2-body decays (excluding neutrinos and photons)
  // ---------------------------------------------------
  if( nFirst == 2 ) {
   if(       _npi==1 && _nf1==1      ) { _mode=1010; }   // pi- f1
    else if( _nK==1 && _neta==1      ) { _mode=1020; }   // K- eta
    else if( _nK==1 && _nomega==1    ) { _mode=1030; }   // K- omega
    else if( _nKstarp==1 && _neta==1 ) { _mode=1040; }   // K*- eta
    else { _mode=-1050; }   // other 2-body decays

   if(_mode==0 )  _mode = 2000;
  }


  // ---------------------------------------------------
  // 3-body decays (excluding neutrinos and photons)
  // ---------------------------------------------------
  if( nFirst == 3 ) {
   if(  _npi==1 && _npi0==1 && _nomega==1 ) { _mode=1110; }   // pi- pi0 omega
   if(  _npi==1 && _npi0==2 )               { _mode=1120; }   // pi- 2pi0 

   if( _mode==0 ) _mode = 3000;
  }

  // ---------------------------------------------------
  // 4-body decays (excluding neutrinos and photons)
  // ---------------------------------------------------
  if( nFirst == 4 ) {
   if(  _npi==3 && _neta==1 )             { _mode=1210; }   // pi- pi- pi+ eta
   if(  _npi==1 && _npi0==2 && _neta==1 ) { _mode=1220; }   // pi- 2pi0  eta
   if(  _npi==1 && _npi0==2 && _nomega==1){ _mode=1230; }   // pi- 2pi0  omega
   if(  _npi==1 && _npi0==1 && _nK==2   ) { _mode=1240; }   // pi- K+ K- pi0  eta (** check charge state)
   if(  _nK==1  && _npi0==3 )             { _mode=1250; }   // pi- K+ K- pi0  eta (** check charge state)

   if(  _npi==3 && _nK0==1  )             { _mode=1260; }   // pi- pi- pi+ K0
   if(  _npi==2 && _npi0==1 && _nK==1 )   { _mode=1270; }   // pi- pi+ pi0 K-     (** check charge state)
   if(  _npi==1 && _npi0==1 && _nK0==2 )  { _mode=1280; }   // pi- pi0 K0 K0
   if(  _npi==1 && _npi0==2 && _nK0==1 )  { _mode=1290; }   // pi- pi0 pi0 K0

   if( _mode==0 ) _mode = 4000;
  }


  // ---------------------------------------------------
  // 5-body decays (excluding neutrinos and photons)
  // ---------------------------------------------------
  if( nFirst == 5 ) {
   if(  _npi==1 && _npi0==4 )            { _mode=1310; }   // pi- 4pi0

   if( _mode==0 ) _mode = 5000;
  }


  if( _mode==0) _mode=9999;

  return _mode;

// ***************************************************************************************


  // K-(321), K0(311) and K*-(892)(323) are primary final states

  if( _nKL!=0 || _nKS!=0 || _nKstar0!=0 )
    std::cout << "===>RJS: Strange Decay in primary chain XXX" << std::endl;

  // check for multiple resonances in final state
  if( _nK0>1 || _nKstarp>1 || _na1>1 || _nw>1 ) 
    std::cout << "===>RJS: multple particles in final state (K0,K*-,a1,W)= "
              << _nK0 << _nKstarp << _na1 << _nw << std::endl;

  // what are the eta(221), omega(223) and f1(20223) decays? (K-eta)
  // (pi- pi0 eta), (pi- 2pi0 eta),(3pi- eta), (K- eta), (K*- eta)
  // (pi- pi0 omega), (K- omega)
  // (pi- f1)

  if( _nomega==1 || _neta==1 || _nf1==1 ) 
   std::cout << "===>RJS: what are these decays? " <<_nomega << _neta << _nf1 << std::endl;

  // test new function
  if( _na1==1) countMesons( tauPtr, resPtr );


  // -------------------------------------------
  // identify primary decays modes
  // -------------------------------------------
  _mode = 0;
  bool electron    = _nnue==1 && _nele==1 && _nnum==0 && _nmu==0;
  bool muon        = _nnue==0 && _nele==0 && _nnum==1 && _nmu==1;
  bool noLepton    = _nnue==0 && _nele==0 && _nnum==0 && _nmu==0;
  bool noMeson     = _npi==0  && _npi0==0 && _nK==0   && _nrho==0;
  bool noResonance = _na1==0  && _neta==0 && _nomega==0 && _nf1==0 && _nK0==0 && _nKstarp==0 && _nw==0;


  if( noResonance ) {
    if(       noMeson && electron  )                                   { _mode=10; } // e
     else if( noMeson && muon )                                        { _mode=20; } // mu
     else if( noLepton && _npi==1 && _npi0==0 && _nK==0 && _nrho==0 )  { _mode=30; } // pi
     else if( noLepton && _npi==0 && _npi0==0 && _nK==1 && _nrho==0 )  { _mode=40; } // K
     else if( noLepton && _npi==0 && _npi0==0 && _nK==0 && _nrho==1 )  { _mode=50; } // rho
     else {std::cout << "===>RJS: unknown noResonant decay mode" << std::endl;}
     return _mode;
  }



// ----------------------------------------------------------------
// Resonance decays
// (not handling 2KS decays)
// ----------------------------------------------------------------
  _Strange = false;
  _nrespi=0;  _nrespi0=0; _nresgam=0; _nresK=0; _nresKS=0; _nresKL=0; _nreseta=0; _nresetap=0;
  if( _nw==1 || _na1==1 || _nKS==1 || _nKstar0==1 || _nKstarp==1) {
     std::vector<Belle2::MCParticle*> rDaugVector = tauDaugVector[resPtr]->getDaughters();  
     for(unsigned int i=0; i<rDaugVector.size(); i++) {
        std::cout << i << " pid = " << rDaugVector[i]->getPDG() << std::endl;
        if( abs(rDaugVector[i]->getPDG()) == abs(PdtLund::pi_plus   ) ) _nrespi++;
        if( abs(rDaugVector[i]->getPDG()) == abs(PdtLund::pi0       ) ) _nrespi0++;
        if( abs(rDaugVector[i]->getPDG()) == abs(PdtLund::gamma     ) ) _nresgam++;
        if( abs(rDaugVector[i]->getPDG()) == abs(PdtLund::eta       ) ) _nreseta++;
        if( abs(rDaugVector[i]->getPDG()) == abs(PdtLund::eta_prime ) ) _nresetap++;
        if( abs(rDaugVector[i]->getPDG()) == abs(PdtLund::K_plus    ) ) _nresK++;
        if( abs(rDaugVector[i]->getPDG()) == abs(PdtLund::K_S0      ) ) _nresKS++;
        if( abs(rDaugVector[i]->getPDG()) == abs(PdtLund::K_L0      ) ) _nresKL++;
     }
   if( _nresK!=0 || _nresKS!=0 || _nresKL!=0) _Strange = true;
       std::cout << std::setw(3) << std::setprecision(0);
       std::cout << "Resonance (g) (pi pi0 eta etap) (K KS KL) = " 
             << std::setw(3) << _nresgam 
             << std::setw(3) << _nrespi << _nrespi0 << _nreseta << _nresetap
             << std::setw(3) << _nresK  << _nresKS << _nresKL
             << " Strange = " << _Strange
             << std::endl;

    // -------------------------------------------
    // W -> non-strange decays
    // -------------------------------------------
    if( _nw == 1 && !_Strange ) {

      if( _nrespi==1 ) {
          if(      _nrespi0==2 )                 { return (int) 21; } // W -> pi 2pi0
          else if( _nrespi0==3 )                 { return (int) 22; } // W -> pi 3pi0
          else if( _nrespi0>=4 )                 { return (int) 23; } // W -> pi 4pi0
          else if( _nrespi0==1 && _nresgam>=1  ) { return (int) 24; } // W -> pi- pi0 gamma
          else if( _nrespi0==1 && _nreseta==1  ) { return (int) 25; } // W -> eta pi- pi0
          else if( _nrespi0==1 && _nresetap==1 ) { return (int) 26; } // W -> eta' pi- 
          else { return (int)-27; }                   // other W -> pi decays
       }  else if( _nrespi==3 ) {
          if(      _nrespi0==1 ) { return (int) 31; } // W -> 3pi pi0
          else if( _nrespi0==2 ) { return (int) 32; } // W -> 3pi 2pi0
          else if( _nrespi0>=3 ) { return (int) 33; } // W -> 3pi 3pi0
          else { return (int)-34; }                   // other W -> 3pi decays
       }  else if( _nrespi==5 ) {
          if(      _nrespi0==0 ) { return (int) 35; } // W -> 5pi
          else if( _nrespi0==1 ) { return (int) 36; } // W -> 5pi pi0
          else if( _nrespi0>=2 ) { return (int) 37; } // W -> 5pi 2pi0
          else { return (int)-38; }                   // other W -> 5pi decays
       }  else { return (int)-39; }                   // unknown W decays

    } // end of non-strange W decays

    // -------------------------------------------
    // W -> strange decays (excl K* mesons)
    // -------------------------------------------
    if( _nw == 1 && _Strange ) {

      // *** this will include all charge combinations eg (K- pi- pi+)  ***
      if( _nresK==1 && _nresKS==0 && _nresKL==0 && _nrespi==2 ) {
          if(       _nrespi0==0 )  { return (int)71; }  // W -> K+ pi- pi-
           else if( _nrespi0==1 )  { return (int)72; }  // W -> K+ pi- pi- pi0
           else { return (int)-73; } // unknown decay
       }

      // *** this will include all charge combinations eg (K- K- pi+) ***
       if( _nresK==2 && _nresKS==0 && _nresKL==0 && _nrespi==1 ) {
          if(       _nrespi0==0 )  { return (int)74; }  // W -> K+ K- pi-
           else if( _nrespi0==1 )  { return (int)75; }  // W -> K+ K- pi- pi0
           else { return (int)-76; } // unknown decay
       }

     if( _nresK==3 && _nresKS==0 && _nresKL==0 && _nrespi==0 ) {
          if(       _nrespi0==0 )  { return (int)77; }  // W -> K+ K- K-
           else if( _nrespi0==1 )  { return (int)78; }  // W -> K+ K- K- pi0
           else { return (int)-79; } // unknown decay
       }

      if( _nresK==1 && _nresKS==1 && _nresKL==0 && _nrespi==0 ) {
          if(       _nrespi0==0 )  { return (int)81; }  // W -> K+ KS
           else if( _nrespi0==1 )  { return (int)82; }  // W -> K+ KS pi0 
           else { return (int)-83; } // unknown decay
       }

       if( _nresK==1 && _nresKL==1 && _nresKS==0 && _nrespi==0 ) {
          if(       _nrespi0==0 )  { return (int)84; }  // W -> K+ KL 
           else if( _nrespi0==1 )  { return (int)85; }  // W -> K+ KL pi0
           else { return (int)-86; } // unknown decay
       }

       if( _nresK==0 && _nresKL==1 && _nresKS==1 && _nrespi==1 ) {
          if(       _nrespi0==0 )  { return (int)87; }  // W -> pi- KS KL (non-resonant)
           else if( _nrespi0==1 )  { return (int)88; }  // W -> pi- KS KL pi0 (non-resonant)
           else { return (int)-89; } // unknown decay
       }

       if( _nresK==1 && _nresKS==1 && _nresKL==0 && _nrespi==0 ) {
          if(       _nrespi0==0 )  { return (int)91; }  // W -> K- KS
           else if( _nrespi0==1 )  { return (int)92; }  // W -> K- KS pi0
           else { return (int)-93; } // unknown decay
       }

      if( _nresK==1 && _nresKS==0 && _nresKL==1 && _nrespi==0 ) {
          if(       _nrespi0==0 )  { return (int)94; }  // W -> K- KL
           else if( _nrespi0==1 )  { return (int)95; }  // W -> K- KL pi0
           else { return (int)-96; } // unknown decay
       }

 
    } // end of strange W decays

    // -------------------------------------------
    // a1 decays
    // -------------------------------------------
    if( _na1==1) {
        if(      _nrespi==3 && _nrespi0==0 ) { return (int) 15; } // a1- -> 2pi- pi+
        else if( _nrespi==1 && _nrespi0==2 ) { return (int) 16; } // a1- -> pi- 2pi0
        else { return (int)-17; }  // other a1 decays
    }

   // -------------------------------------------
   // KS decays
   // -------------------------------------------
   if( _nKS==1) {
        if(      _nrespi==2 && _nrespi0==0 ) { return (int) 31; } // KS -> pi- pi+
        else if( _nrespi==0 && _nrespi0==2 ) { return (int) 32; } // KS -> 2pi0
        else { return (int)-33; }  // other KS decays
    }
   if( _nKS==2) return (int)-34; // 2KS decay mode
   

   // -------------------------------------------
   // Kstar decays
   // -------------------------------------------
   if( _nKstar0==1) {
        if(      _nresK==1  && _nrespi==1 && _nrespi0==0 ) { return (int) 51; } // Kstar0 -> K- pi-
        else if( _nresKS==1 && _nrespi==0 && _nrespi0==1 ) { return (int) 52; } // Kstar0 -> KS0 pi0
        else if( _nresKL==1 && _nrespi==0 && _nrespi0==1 ) { return (int) 53; } // Kstar0 -> KL0 pi0
        else { return (int)-54; }  // other Kstar0 decays
    }

   if( _nKstarp==1) {
        if(      _nresK ==1 && _nrespi==0 && _nrespi0==1 ) { return (int) 61; } // Kstar+ -> K- pi0
        else if( _nresKS==1 && _nrespi==1 && _nrespi0==0 ) { return (int) 62; } // Kstar+ -> KS0 pi-
        else if( _nresKL==1 && _nrespi==1 && _nrespi0==0 ) { return (int) 63; } // Kstar+ -> KL0 pi-
        else { return (int)-64; }  // other Kstar0 decays
    }

  } //end of resonance section

   
// ----------------------------------------------------------------
// ----------------------------------------------------------------
  if( _nnue==1 && _nele==1)  return (int)1;   // e 
  if( _nnum==1 && _nmu==1 )  return (int)2;   // mu
  if( _npi==1)               return (int)3;   // pi-
  if( _nrho==1)              return (int)4;   // rho-

  // -------------------------------------------
  // charged kaon decays
  // -------------------------------------------
  if( _nK==1 && _npi==0 && _npi0==0) return (int)40;  // K-


  return (int)99;
}

// *******************************************************************
// Count the types of daughter particles
// *******************************************************************
int TauDecayMode::countParticles(   std::vector<Belle2::MCParticle*> tauDaugVector ) 
{
   // ---------------------------------------------
   // get the list of daughter particles
   // ---------------------------------------------
//   Belle2::StoreArray<Belle2::MCParticle> mcList;
//   std::vector<Belle2::MCParticle*> tauDaugVector = mcList[tauPtr]->getDaughters();

   // ---------------------------------------------
   // initialize counters
   // ---------------------------------------------
   int ptr = 0;

   // photons
   _ng = 0;
   // leptons
   _nnue=0;   _nnum=0;  _nnut=0; _nele=0;  _nmu=0;
   // primary mesons
   _npi=0;    _npip=0;  _npim=0;
   _npi0=0;   _nrho=0;  _na1=0;  
   _nomega=0; _neta=0;  _nf1=0;
   // strange mesons
   _nK=0;       _nKp=0;     _nKm=0;
   _nK0=0;      _nKL=0;     _nKS=0;   
   _nKstar0=0;  _nKstarp=0;

  for(unsigned int i=0; i<tauDaugVector.size(); i++) {
     if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::gamma   ) ) {  _ng++; }

     else if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::nu_e    ) ) { _nnue++; }
     else if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::nu_mu   ) ) { _nnum++; }
     else if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::nu_tau  ) ) { _nnut++; }

     else if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::e_plus  ) ) { _nele++; }
     else if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::mu_plus ) ) { _nmu++;  }
     else if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::pi0     ) ) { _npi0++; }
     else if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::rho_plus) ) { _nrho++; }
     else if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::omega   ) ) { _nomega++;   }
     else if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::eta     ) ) { _neta++;   }
     else if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::f_1     ) ) { _nf1++;   }

     else if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::pi_plus) ) {
          _npi++;
          if( tauDaugVector[i]->getPDG() == PdtLund::pi_plus       ) _npip++;
          if( tauDaugVector[i]->getPDG() == PdtLund::pi_minus      ) _npim++;
          }

     else if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::K_plus) )  {
          _nK++;
          if( tauDaugVector[i]->getPDG() == PdtLund::K_plus        ) _nKp++;
          if( tauDaugVector[i]->getPDG() == PdtLund::K_minus       ) _nKm++;
          }

     else if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::K_S0    ) ) _nKS++;
     else if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::K_L0    ) ) _nKL++;

     else if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::K0      ) )     { _nK0++;     ptr = i; }
     else if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::K_star_plus ) ) { _nKstarp++; ptr = i; }
     else if( abs(tauDaugVector[i]->getPDG()) == abs(PdtLund::a_1_plus) )     { _na1++;     ptr = i; }

     else { std::cout << "===>RJS: CountParticle--unknown particle " << tauDaugVector[i]->getPDG() << std::endl; }
  }

  return ptr;
}
