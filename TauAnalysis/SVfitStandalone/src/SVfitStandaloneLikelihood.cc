#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneLikelihood.h"

#include "TauAnalysis/SVfitStandalone/interface/svFitStandaloneAuxFunctions.h"
#include "TauAnalysis/SVfitStandalone/interface/LikelihoodFunctions.h"

using namespace svFitStandalone;

/// indicate first iteration for integration or fit cycle for debugging
static bool FIRST = true;

SVfitStandaloneLikelihood::SVfitStandaloneLikelihood(const std::vector<MeasuredTauLepton>& measuredTauLeptons, const Vector& measuredMET, const TMatrixD& covMET, bool verbose) 
  : metPower_(1.0), 
    addLogM_(false), 
    powerLogM_(1.),
    addDelta_(true),
    addSinTheta_(false),
    addPhiPenalty_(true),
    verbose_(verbose), 
    idxObjFunctionCall_(0), 
    invCovMET_(2,2),
    errorCode_(0),
    requirePhysicalSolution_(false),
    marginalizeVisMass_(false),
    l1lutVisMass_(0),
    l2lutVisMass_(0),
    shiftVisMass_(false),
    l1lutVisMassRes_(0),
    l2lutVisMassRes_(0),
    shiftVisPt_(false),
    l1lutVisPtRes_(0),
    l2lutVisPtRes_(0)
{
  if ( verbose_ ) {
    std::cout << "<SVfitStandaloneLikelihood::SVfitStandaloneLikelihood>:" << std::endl;
  }
  measuredMET_ = measuredMET;
  // for integration mode the order of lepton or tau matters due to the choice of order in which 
  // way the integration boundaries are defined. In this case the lepton should always go before
  // the tau. For tau-tau or lep-lep the order is irrelevant.
  if ( measuredTauLeptons[0].type() == svFitStandalone::kTauToHadDecay ) {
    measuredTauLeptons_.push_back(measuredTauLeptons[1]);
    measuredTauLeptons_.push_back(measuredTauLeptons[0]);
  } else {
    measuredTauLeptons_= measuredTauLeptons;
  }
  if ( measuredTauLeptons_.size() != 2 ) {
    std::cout << " >> ERROR : the number of measured leptons must be 2 but is found to be: " << measuredTauLeptons_.size() << std::endl;
    errorCode_ |= LeptonNumber;
  }
  // determine transfer matrix for MET
  invCovMET_= covMET;
  covDet_ = invCovMET_.Determinant();
  if ( covDet_ != 0 ) { 
    invCovMET_.Invert(); 
  } else{
    std::cout << " >> ERROR: cannot invert MET covariance Matrix (det=0)." << std::endl;
    errorCode_ |= MatrixInversion;
  }
}

void 
SVfitStandaloneLikelihood::marginalizeVisMass(bool value, const TH1* l1lutVisMass, const TH1* l2lutVisMass)
{
  marginalizeVisMass_ = value;
  if ( marginalizeVisMass_ ) {
    l1lutVisMass_ = l1lutVisMass;
    l2lutVisMass_ = l2lutVisMass;
  }
}

void 
SVfitStandaloneLikelihood::shiftVisMass(bool value, const TH1* l1lutVisMassRes, const TH1* l2lutVisMassRes)
{
  shiftVisMass_ = value;
  if ( shiftVisMass_ ) {
    l1lutVisMassRes_ = l1lutVisMassRes;
    l2lutVisMassRes_ = l2lutVisMassRes;
  }
}

void 
SVfitStandaloneLikelihood::shiftVisPt(bool value, const TH1* l1lutVisPtRes, const TH1* l2lutVisPtRes)
{
  shiftVisPt_ = value;
  if ( shiftVisPt_ ) {
    l1lutVisPtRes_ = l1lutVisPtRes;
    l2lutVisPtRes_ = l2lutVisPtRes;
  }
}

const double*
SVfitStandaloneLikelihood::transform(double* xPrime, const double* x, bool fixToMtest, double mtest) const
{
  // if ( verbose_ ) {
  //   std::cout << "<SVfitStandaloneLikelihood:transform(double*, const double*)>:" << std::endl;
  // }
  LorentzVector fittedDiTauSystem;
  for ( size_t idx = 0; idx < measuredTauLeptons_.size(); ++idx ) {
    const MeasuredTauLepton& measuredTauLepton = measuredTauLeptons_[idx];

    // map to local variables to be more clear on the meaning of the individual parameters. The fit parameters are ayered 
    // for each tau decay
    double nunuMass, labframeXFrac, labframePhi;
    double visMass_unshifted = measuredTauLepton.mass();
    double visMass = visMass_unshifted; // visible momentum in lab-frame
    double labframeVisMom_unshifted = measuredTauLepton.momentum(); 
    double labframeVisMom = labframeVisMom_unshifted; // visible momentum in lab-frame
    double labframeVisEn  = measuredTauLepton.energy(); // visible energy in lab-frame    
    if ( measuredTauLepton.type() == kTauToElecDecay || measuredTauLepton.type() == kTauToMuDecay ) {
      labframeXFrac = x[idx*kMaxFitParams + kXFrac];
      nunuMass = x[idx*kMaxFitParams + kMNuNu];
      labframePhi = x[idx*kMaxFitParams + kPhi];
    } else {
      labframeXFrac = x[idx*kMaxFitParams + kXFrac];
      nunuMass = 0.;
      labframePhi = x[idx*kMaxFitParams + kPhi];
      if ( marginalizeVisMass_ || shiftVisMass_ ) {
	visMass = x[idx*kMaxFitParams + kVisMassShifted];
      } 
      if ( shiftVisPt_ ) {
	double shiftInv = 1. + x[idx*kMaxFitParams + kRecTauPtDivGenTauPt];
	double shift = ( shiftInv > 1.e-1 ) ?
	  (1./shiftInv) : 1.e+1;
	labframeVisMom *= shift;
	//visMass *= shift; // CV: take mass and momentum to be correlated
	//labframeVisEn = TMath::Sqrt(labframeVisMom*labframeVisMom + visMass*visMass);
	labframeVisEn *= shift;
      }
    }
    bool isValidSolution = true;
    // add protection against unphysical mass of visible tau decay products
    if ( visMass < electronMass || visMass > tauLeptonMass ) { 
      isValidSolution = false;
    }  
    // add protection against unphysical visible energy fractions
    if ( !(labframeXFrac >= 0. && labframeXFrac <= 1.) ) {
      isValidSolution = false;
    }
    // CV: do not spend time on unphysical solutions: returning 0 pointer will lead to 0 evaluation of prob
    if ( !isValidSolution ) {
      return 0;
    }
    double gjAngle_lab = gjAngleLabFrameFromX(labframeXFrac, visMass, nunuMass, labframeVisMom, labframeVisEn, tauLeptonMass, isValidSolution);
    double enTau_lab = labframeVisEn/labframeXFrac;
    if ( (enTau_lab*enTau_lab) < tauLeptonMass2 ) {
      enTau_lab = tauLeptonMass;
      isValidSolution = false;
    }  
    double pTau_lab = TMath::Sqrt(square(enTau_lab) - tauLeptonMass2);
    double gamma = enTau_lab/tauLeptonMass;
    double beta = TMath::Sqrt(1. - 1./(gamma*gamma));
    double pVis_parl_rf = -beta*gamma*labframeVisEn + gamma*TMath::Cos(gjAngle_lab)*labframeVisMom;
    double pVis_perp = labframeVisMom*TMath::Sin(gjAngle_lab);
    double gjAngle_rf = TMath::ATan2(pVis_perp, pVis_parl_rf);
    Vector p3Tau_unit = motherDirection(measuredTauLepton.direction(), gjAngle_lab, labframePhi);
    LorentzVector p4Tau_lab = motherP4(p3Tau_unit, pTau_lab, enTau_lab);
    fittedDiTauSystem += p4Tau_lab;
    // fill branch-wise nll parameters
    xPrime[ idx == 0 ? kNuNuMass1            : kNuNuMass2            ] = nunuMass;
    xPrime[ idx == 0 ? kVisMass1             : kVisMass2             ] = visMass;
    xPrime[ idx == 0 ? kDecayAngle1          : kDecayAngle2          ] = gjAngle_rf;
    xPrime[ idx == 0 ? kDeltaVisMass1        : kDeltaVisMass2        ] = ( marginalizeVisMass_ ) ? visMass : (visMass_unshifted - visMass);
    xPrime[ idx == 0 ? kRecTauPtDivGenTauPt1 : kRecTauPtDivGenTauPt2 ] = ( labframeVisMom > 0. ) ? (labframeVisMom_unshifted/labframeVisMom) : 1.e+3;
    xPrime[ idx == 0 ? kMaxNLLParams         : (kMaxNLLParams + 1)   ] = labframeXFrac;
    xPrime[ idx == 0 ? (kMaxNLLParams + 2)   : (kMaxNLLParams + 3)   ] = isValidSolution;
  }
 
  Vector fittedMET = fittedDiTauSystem.Vect() - (measuredTauLeptons_[0].p()+measuredTauLeptons_[1].p()); 
  // fill event-wise nll parameters
  xPrime[ kDMETx   ] = measuredMET_.x() - fittedMET.x(); 
  xPrime[ kDMETy   ] = measuredMET_.y() - fittedMET.y();
  if ( fixToMtest ) xPrime[ kMTauTau ] = mtest;       // CV: evaluate delta-function derrivate in case of VEGAS integration for nominal test mass,
  else xPrime[ kMTauTau ] = fittedDiTauSystem.mass(); //     not for fitted mass, to improve numerical stability of integration (this is what the SVfit plugin version does)

  if ( verbose_ && FIRST ) {
    std::cout << " >> input values for transformed variables: " << std::endl;
    std::cout << "    MET[x] = " <<  fittedMET.x() << " (fitted)  " << measuredMET_.x() << " (measured) " << std::endl; 
    std::cout << "    MET[y] = " <<  fittedMET.y() << " (fitted)  " << measuredMET_.y() << " (measured) " << std::endl; 
    std::cout << "    fittedDiTauSystem: [" 
	      << " px = " << fittedDiTauSystem.px() 
	      << " py = " << fittedDiTauSystem.py() 
	      << " pz = " << fittedDiTauSystem.pz() 
	      << " En = " << fittedDiTauSystem.energy() 
	      << " ]" << std::endl; 
    std::cout << " >> nll parameters after transformation: " << std::endl;
    std::cout << "    x[kNuNuMass1  ] = " << xPrime[kNuNuMass1  ] << std::endl;
    std::cout << "    x[kVisMass1   ] = " << xPrime[kVisMass1   ] << std::endl;
    std::cout << "    x[kDecayAngle1] = " << xPrime[kDecayAngle1] << std::endl;
    std::cout << "    x[kNuNuMass2  ] = " << xPrime[kNuNuMass2  ] << std::endl;
    std::cout << "    x[kVisMass2   ] = " << xPrime[kVisMass2   ] << std::endl;
    std::cout << "    x[kDecayAngle2] = " << xPrime[kDecayAngle2] << std::endl;
    std::cout << "    x[kDMETx      ] = " << xPrime[kDMETx      ] << std::endl;
    std::cout << "    x[kDMETy      ] = " << xPrime[kDMETy      ] << std::endl;
    std::cout << "    x[kMTauTau    ] = " << xPrime[kMTauTau    ] << std::endl;
  }
  return xPrime;
}

double
SVfitStandaloneLikelihood::prob(const double* x, bool fixToMtest, double mtest) const 
{
  // in case of initialization errors don't start to do anything
  if ( error() ) { 
    return 0.;
  }
  // if ( verbose_ ) {
  //   std::cout << "<SVfitStandaloneLikelihood:prob(const double*)>:" << std::endl;
  // }
  ++idxObjFunctionCall_;
  // if ( verbose_ && FIRST ) {
  //   std::cout << " >> ixdObjFunctionCall : " << idxObjFunctionCall_ << std::endl;  
  // }
  // prevent kPhi in the fit parameters (kFitParams) from trespassing the 
  // +/-pi boundaries
  double phiPenalty = 0.;
  if ( addPhiPenalty_ ) {
    for ( size_t idx = 0; idx < measuredTauLeptons_.size(); ++idx ) {
      if ( TMath::Abs(idx*kMaxFitParams + x[kPhi]) > TMath::Pi() ) {
	phiPenalty += (TMath::Abs(x[kPhi]) - TMath::Pi())*(TMath::Abs(x[kPhi]) - TMath::Pi());
      }
    }
  }
  // xPrime are the transformed variables from which to construct the nll
  // transform performs the transformation from the fit parameters x to the 
  // nll parameters xPrime. prob is the actual combined likelihood. The
  // phiPenalty prevents the fit to converge to unphysical values beyond
  // +/-pi 

//   double xPrime[kMaxNLLParams + 4];
  double xPrime[kMaxNLLParams + 4] = {};
  
//   double myxPrime = 0. ; 
//   double* xPrime = &myxPrime;


//   double myx[] = {0.10015, 0.461168, -2.59906, 0., 0., 0.252015, 0., -2.75104, 0., 0., 3.40126e-18, 54.6908, 344.251};
  
  const double* xPrime_ptr = transform(xPrime, x, fixToMtest, mtest);
//   const double* xPrime_ptr = transform(xPrime, myx, fixToMtest, mtest);

  // 3031041

//   if (x[12]>344 && x[12]<352){// && idxObjFunctionCall_ > 0 && idxObjFunctionCall_ < 5){
  if (x[12]>344 && x[12]<352 && idxObjFunctionCall_ > 3031040 && idxObjFunctionCall_ < 3031045){
    std::cout << "\n\n" << __LINE__ << "]\t" << __PRETTY_FUNCTION__  
              << "\nidxObjFunctionCall_ " << idxObjFunctionCall_
              << "\nkMaxNLLParams "       << kMaxNLLParams
              << "\nx[0] "                << x[0]
              << "\nx[1] "                << x[1]
              << "\nx[2] "                << x[2]
              << "\nx[3] "                << x[3]
              << "\nx[4] "                << x[4]
              << "\nx[5] "                << x[5]
              << "\nx[6] "                << x[6]
              << "\nx[7] "                << x[7]
              << "\nx[8] "                << x[8]
              << "\nx[9] "                << x[9]
              << "\nx[10] "               << x[10]
              << "\nx[11] "               << x[11]
              << "\nx[12] "               << x[12]
              << "\nxPrime[0] "           << xPrime[0]
              << "\nxPrime[1] "           << xPrime[1]
              << "\nxPrime[2] "           << xPrime[2]
              << "\nxPrime[3] "           << xPrime[3]
              << "\nxPrime[4] "           << xPrime[4]
              << "\nxPrime[5] "           << xPrime[5]
              << "\nxPrime[6] "           << xPrime[6]
              << "\nxPrime[7] "           << xPrime[7]
              << "\nxPrime[8] "           << xPrime[8]
              << "\nxPrime[9] "           << xPrime[9]
              << "\nxPrime[10] "          << xPrime[10]
              << "\nxPrime[11] "          << xPrime[11]
              << "\nxPrime[12] "          << xPrime[12]
              << "\nxPrime[13] "          << xPrime[13]
              << "\nxPrime[14] "          << xPrime[14]
              << "\nxPrime[15] "          << xPrime[15]
              << "\nxPrime[16] "          << xPrime[16]
//               << "\nxPrime_ptr[0] "       << xPrime_ptr[0]
//               << "\nxPrime_ptr[1] "       << xPrime_ptr[1]
//               << "\nxPrime_ptr[2] "       << xPrime_ptr[2]
//               << "\nxPrime_ptr[3] "       << xPrime_ptr[3]
//               << "\nxPrime_ptr[4] "       << xPrime_ptr[4]
//               << "\nxPrime_ptr[5] "       << xPrime_ptr[5]
//               << "\nxPrime_ptr[6] "       << xPrime_ptr[6]
//               << "\nxPrime_ptr[7] "       << xPrime_ptr[7]
//               << "\nxPrime_ptr[8] "       << xPrime_ptr[8]
//               << "\nxPrime_ptr[9] "       << xPrime_ptr[9]
//               << "\nxPrime_ptr[10] "      << xPrime_ptr[10]
//               << "\nxPrime_ptr[11] "      << xPrime_ptr[11]
//               << "\nxPrime_ptr[12] "      << xPrime_ptr[12]
//               << "\nxPrime_ptr[13] "      << xPrime_ptr[13]
//               << "\nxPrime_ptr[14] "      << xPrime_ptr[14]
//               << "\nxPrime_ptr[15] "      << xPrime_ptr[15]
//               << "\nxPrime_ptr[16] "      << xPrime_ptr[16]
              << "\nfixToMtest "          << fixToMtest
              << "\nmtest "               << mtest
              << std::endl;
  }

//   if (idxObjFunctionCall_== 3031041){  
//     xPrime[0]  = 0.461168;
//     xPrime[1]  = 0.10566;
//     xPrime[2]  = 2.48565;
//     xPrime[3]  = 0;
//     xPrime[4]  = 1;
//     xPrime[5]  = 0;
//     xPrime[6]  = 0.5585;
//     xPrime[7]  = 2.29173;
//     xPrime[8]  = 0;
//     xPrime[9]  = 1;
//     xPrime[10] = -3.38972;
//     xPrime[11] = 30.4944;
//     xPrime[12] = 344.251;
//     xPrime[13] = 0.10015;
//     xPrime[14] = 0.252015;
//     xPrime[15] = 1;
//     xPrime[16] = 1;
//     
//   }

  if ( xPrime_ptr ) {
    return prob(xPrime_ptr, phiPenalty);
  } else {
    return 0.;
  }
}

double 
SVfitStandaloneLikelihood::prob(const double* xPrime, double phiPenalty) const
{
  if ( verbose_ && FIRST ) {
    std::cout << "<SVfitStandaloneLikelihood:prob(const double*, double)>:" << std::endl;
  }
  if ( requirePhysicalSolution_ && (xPrime[ kMaxNLLParams + 2 ] < 0.5 || xPrime[ kMaxNLLParams + 3 ] < 0.5) ) return 0.;
  // start the combined likelihood construction from MET
  double prob = probMET(xPrime[kDMETx], xPrime[kDMETy], covDet_, invCovMET_, metPower_, (verbose_&& FIRST));
  if ( verbose_ && FIRST ) {
    std::cout << "probMET         = " << prob << std::endl;
  }
  // add likelihoods for the decay branches
  for ( size_t idx = 0; idx < measuredTauLeptons_.size(); ++idx ) {
    switch ( measuredTauLeptons_[idx].type() ) {
    case kTauToHadDecay :
      prob *= probTauToHadPhaseSpace(
                xPrime[idx == 0 ? kDecayAngle1 : kDecayAngle2], 
		xPrime[idx == 0 ? kNuNuMass1 : kNuNuMass2], 
		xPrime[idx == 0 ? kVisMass1 : kVisMass2], 
		xPrime[idx == 0 ? kMaxNLLParams : (kMaxNLLParams + 1)], 
		addSinTheta_, 
		(verbose_&& FIRST));
      assert(!(marginalizeVisMass_ && shiftVisMass_));
      if ( marginalizeVisMass_ ) {
	prob *= probVisMass(
                  xPrime[idx == 0 ? kDeltaVisMass1 : kDeltaVisMass2], 
		  idx == 0 ? l1lutVisMass_ : l2lutVisMass_,
		  (verbose_&& FIRST));
      }
      if ( shiftVisMass_ ) {
	prob *= probVisMassShift(
                  xPrime[idx == 0 ? kDeltaVisMass1 : kDeltaVisMass2], 
		  idx == 0 ? l1lutVisMassRes_ : l2lutVisMassRes_,
		  (verbose_&& FIRST));
      }
      if ( shiftVisPt_ ) {
	prob *= probVisPtShift(
		  xPrime[idx == 0 ? kRecTauPtDivGenTauPt1 : kRecTauPtDivGenTauPt2], 
		  idx == 0 ? l1lutVisPtRes_ : l2lutVisPtRes_, 
		  (verbose_&& FIRST));
      }
      if ( verbose_ && FIRST ) {
	std::cout << " *probTauToHad  = " << prob << std::endl;
      }
      break;
    case kTauToElecDecay :
    case kTauToMuDecay :
      prob *= probTauToLepMatrixElement(
		xPrime[idx == 0 ? kDecayAngle1 : kDecayAngle2], 
		xPrime[idx == 0 ? kNuNuMass1 : kNuNuMass2], 
		xPrime[idx == 0 ? kVisMass1 : kVisMass2], 
		xPrime[idx == 0 ? kMaxNLLParams : (kMaxNLLParams + 1)], 
		addSinTheta_, 
		(verbose_&& FIRST));
      if ( verbose_ && FIRST ) {
	std::cout << " *probTauToLep  = " << prob << std::endl;
      }
      break;
    default :
      break;
    }
  }
  // add additional logM term if configured such 
  if ( addLogM_ && powerLogM_ > 0. ) {
    if ( xPrime[kMTauTau] > 0. ) {
      prob *= TMath::Power(1.0/xPrime[kMTauTau], powerLogM_);
    }
    if ( verbose_ && FIRST ) {
      std::cout << " *1/mtautau     = " << prob << std::endl;
    }
  }
  if ( addDelta_ ) {
    prob *= (2.*xPrime[kMaxNLLParams + 1]/xPrime[kMTauTau]);
    if ( verbose_ && FIRST ) {
      std::cout << " *deltaDeriv.   = " << prob << std::endl;
    }
  }
  // add additional phiPenalty in case kPhi in the fit parameters 
  // (kFitParams) trespassed the physical boundaries from +/-pi 
  if ( phiPenalty > 0. ) {
    prob *= TMath::Exp(-phiPenalty);
    if ( verbose_ && FIRST ) {
      std::cout << "* phiPenalty   = " << prob << std::endl;
    }
  }
  // set FIRST to false after the first complete evaluation of the likelihood 
  FIRST = false;
  return prob;
}

void
SVfitStandaloneLikelihood::results(std::vector<LorentzVector>& fittedTauLeptons, const double* x) const
{
  if ( verbose_ ) {
    std::cout << "<SVfitStandaloneLikelihood:results(std::vector<LorentzVector>&, const double*)>:" << std::endl;
  }
  for ( size_t idx = 0; idx < measuredTauLeptons_.size(); ++idx ) {
    const MeasuredTauLepton& measuredTauLepton = measuredTauLeptons_[idx];

    // map to local variables to be more clear on the meaning of the individual parameters. The fit parameters are ayered 
    // for each tau decay
    double nunuMass                 = x[ idx*kMaxFitParams + kMNuNu ];       // nunu inv mass (can be const 0 for had tau decays) 
    double labframeXFrac            = x[ idx*kMaxFitParams + kXFrac ];       // visible energy fraction x in labframe
    double labframePhi              = x[ idx*kMaxFitParams + kPhi   ];       // phi in labframe 
    double visMass                  = measuredTauLepton.mass(); 
    double labframeVisMom_unshifted = measuredTauLepton.momentum(); 
    double labframeVisMom           = labframeVisMom_unshifted; // visible momentum in lab-frame
    double labframeVisEn            = measuredTauLepton.energy(); // visible energy in lab-frame    
    if ( measuredTauLepton.type() == kTauToHadDecay ) {
      if ( marginalizeVisMass_ || shiftVisMass_ ) {
	visMass = x[idx*kMaxFitParams + kVisMassShifted];
      }
      if ( shiftVisPt_ ) {
	double shiftInv = 1. + x[idx*kMaxFitParams + kRecTauPtDivGenTauPt];
        double shift = ( shiftInv > 1.e-1 ) ?
	  (1./shiftInv) : 1.e+1;
        labframeVisMom *= shift;
        //visMass *= shift; // CV: take mass and momentum to be correlated
        //labframeVisEn = TMath::Sqrt(labframeVisMom*labframeVisMom + visMass*visMass);
        labframeVisEn *= shift;
      }
    }
    if ( visMass < 5.1e-4 ) { 
      visMass = 5.1e-4; 
    } 
    bool isValidSolution = true;
    double gjAngle_lab = gjAngleLabFrameFromX(labframeXFrac, visMass, nunuMass, labframeVisMom, labframeVisEn, tauLeptonMass, isValidSolution);
    double enTau_lab = labframeVisEn/labframeXFrac;
    double pTau_lab = TMath::Sqrt(square(enTau_lab) - tauLeptonMass2);
    Vector p3Tau_unit = motherDirection(measuredTauLepton.direction(), gjAngle_lab, labframePhi);
    LorentzVector p4Tau_lab = motherP4(p3Tau_unit, pTau_lab, enTau_lab);
    // tau lepton four vector in labframe
    if ( idx < fittedTauLeptons.size() ) fittedTauLeptons[idx] = p4Tau_lab;
    else fittedTauLeptons.push_back(p4Tau_lab);
  }
}
