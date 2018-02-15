// ROOT stuff
#include <TRandom.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TText.h>
#include <TLine.h>
#include <TMath.h>
#include <TApplication.h>
#include <TTree.h>

// System stuff
#include <fstream>
#include <sys/time.h>
#include <sys/times.h>
#include <random>
#include <CLI/Timer.hpp>
#include <stdio.h>
#include <iostream>
#include <math.h>

// GooFit stuff
#include <goofit/Variable.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/PDFs/physics/DalitzPlotPdf.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/FitManager.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/BinnedDataSet.h>
#include <goofit/Application.h>
#include <goofit/PDFs/basic/SmoothHistogramPdf.h>

using namespace std;
using namespace GooFit;

const int PI = 3.14159265358979323846;

TCanvas* foo;
TCanvas* foodal;
timeval startTime, stopTime, totalTime;
clock_t startCPU, stopCPU;
tms startProc, stopProc;
UnbinnedDataSet* data = 0;
const unsigned int nbins = 1000;
TH2F* weightHistogram = 0;
TH2F* bkgHistogram    = 0;
TH2F* underlyingBins  = 0;

// How many events will be generated for Eff Bkg?
const int NevG = 1e7;

// PWA INPUT FILE NAME
const string  pwa_file = "files/PWA_COEFFS_50.txt";


// FIT OR JUST PLOT?
bool fit = true;

Observable m12("m12", 0.9, 2.0);
Observable m13("m13", 0.9, 2.0);
EventNumber eventNumber("eventNumber");
bool fitMasses = false;
Variable fixedPhiMass("phi_mass", 1.019461, 0.01, 0.7, 1.8);
Variable fixedPhiWidth("phi_width", 0.004266, 0.001, 1e-5, 1e-1);

const fptype _mDp = 1.86962;
const fptype KPlusMass = 0.493677;
double V = (m12.getUpperLimit() - m12.getLowerLimit())*(m13.getUpperLimit() - m13.getLowerLimit()); //Volume


const fptype D1Mass = KPlusMass;//
const fptype D2Mass = KPlusMass;
const fptype D3Mass = KPlusMass;
const fptype D1Mass2 = D1Mass*D1Mass;
const fptype D2Mass2 = D2Mass*D2Mass;
const fptype D3Mass2 = D3Mass*D3Mass;
const fptype MMass = _mDp;
const fptype MMass2 = MMass*MMass;


//const fptype MMass2inv = 1./MMass2;

// Constants used in more than one PDF component.
Variable  motherM("motherM", MMass);
Variable dau1M("dau1M", D1Mass);
Variable dau2M("dau2M", D2Mass);
Variable dau3M("dau3M", D3Mass);
Variable massSum("massSum", MMass2 + D1Mass2+D2Mass2+D3Mass2); // = 3.53481
Variable constantOne("constantOne", 1);
Variable constantZero("constantZero", 0);

std::vector<PdfBase*> comps;

// I don't like Globals! Henry
int verbosity = 3;

GooPdf* kzero_veto = 0;
char strbuffer[1000];
double mesonRad  = 1.5;
DalitzPlotPdf* signalDalitz;
bool doEffSwap = true;
bool saveEffPlot = true;
bool saveBkgPlot = true;

void makeToyDalitzData (GooPdf* overallSignal, const int iSeed = 0, string datadir = ".", const int nTotal = 1.e6 ) ;

DalitzPlotPdf* makeSignalPdf (GooPdf* eff = 0, bool fixAmps = false) ;

fptype cpuGetM23 (fptype massPZ, fptype massPM) {
	return (massSum.getValue() - massPZ - massPM);
}

bool cpuDalitz (fptype m_12, fptype m_13, fptype bigM = MMass, fptype dm1 = D1Mass, fptype dm2 = D2Mass, fptype dm3 = D3Mass) {
	if (m_12 < pow(dm1 + dm2, 2)) return false; // This m_12 cannot exist, it's less than the square of the (1,2) particle mass.
	if (m_12 > pow(bigM - dm3, 2)) return false;   // This doesn't work either, there's no room for an at-rest 3 daughter.

	// Calculate energies of 1 and 3 particles in m_12 rest frame.
	fptype e1star = 0.5 * (m_12 - dm2*dm2 + dm1*dm1) / sqrt(m_12);
	fptype e3star = 0.5 * (bigM*bigM - m_12 - dm3*dm3) / sqrt(m_12);

	// Bounds for m_13 at this value of m_12.
	fptype minimum = pow(e1star + e3star, 2) - pow(sqrt(e1star*e1star - dm1*dm1) + sqrt(e3star*e3star - dm3*dm3), 2);
	if (m_13 < minimum) return false;
	fptype maximum = pow(e1star + e3star, 2) - pow(sqrt(e1star*e1star - dm1*dm1) - sqrt(e3star*e3star - dm3*dm3), 2);
	if (m_13 > maximum) return false;

	return true;
}

void makeToyDalitzData (GooPdf* overallSignal, const int iSeed, string datadir, const int nTotal ) {
	std::vector<Observable> vars;
	vars.push_back(m12);
	vars.push_back(m13);
	vars.push_back(eventNumber);
	data = new UnbinnedDataSet(vars);
	UnbinnedDataSet currData(vars);
	std::vector<std::vector<double>> pdfValues;
	int ncount = 0;
	TRandom3 donram(iSeed);
	for (int i = 0; i < (m12.getNumBins()) ; ++i) {
		m12.setValue( m12.getLowerLimit() + (m12.getUpperLimit() - m12.getLowerLimit())*(i + 0.5) / m12.getNumBins() );
		for (int j = 0; j < m13.getNumBins(); ++j) {
			m13.setValue(m13.getLowerLimit() + (m13.getUpperLimit() - m13.getLowerLimit())*(j + 0.5) / m13.getNumBins());
			if (!cpuDalitz(m12.getValue(), m13.getValue(), MMass , D1Mass, D2Mass,D3Mass)) continue;
			eventNumber.setValue(ncount);
			ncount++;
			currData.addEvent();
		}
	}
	signalDalitz->setDataSize(currData.getNumEvents());
	overallSignal->setData(&currData);

	pdfValues = overallSignal->getCompProbsAtDataPoints();
	TH2F dalitzpp0_dat_hist("dalitzpp0_dat_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit(), m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
	dalitzpp0_dat_hist.GetXaxis()->SetTitle("m^{2}(K^{-}K^{+}) [GeV^{2}]");
	dalitzpp0_dat_hist.GetYaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV^{2}]");
	ncount = 0;
	ofstream writer;
	sprintf(strbuffer, "%s/dalitz_mytoyMC_%03d.txt", datadir.c_str(), iSeed);
	writer.open(strbuffer);
	vector<double> fIntegral;
	fIntegral.push_back(pdfValues[0][0]);
	Int_t ncells = pdfValues[0].size();
	for (unsigned int j = 1; j < ncells; ++j) {
		fIntegral.push_back(pdfValues[0][j]+fIntegral[j-1]);
	}
	for (unsigned int j = 0; j < ncells; ++j)  fIntegral[j] /= fIntegral[ncells-1];
	ncount = 0;
	int nEvents = donram.Poisson(nTotal);
	for (int iEvt = 0;iEvt<nEvents;iEvt++){
		double r = donram.Rndm();
		//Binary search for fIntegral[cell-1] < r < fIntegral[cell]
		int lo = 0, hi = ncells-1, mid = 0;
		while(lo <= hi){
			mid = lo + (hi-lo)/2;
			if( r<=fIntegral[mid]&&(mid==0||r>fIntegral[mid-1])) break;
			else if (r > fIntegral[mid] ) lo = mid+1;
			else hi = mid-1;
		}
		int j = mid;
		double currm12 = currData.getValue(m12, j);
		currm12 += (m12.getUpperLimit() - m12.getLowerLimit())*(donram.Rndm() - 0.5) / m12.getNumBins();
		double currm13 = currData.getValue(m13, j);
		currm13 += (m13.getUpperLimit() - m13.getLowerLimit())*(donram.Rndm() - 0.5) / m13.getNumBins();
		eventNumber.setValue(ncount++);
		dalitzpp0_dat_hist.Fill(currm12, currm13);
		data->addEvent();
		writer << ncount-1 << '\t'<<currm12 << '\t'<<currm13<<std::endl;
	}
	writer.close();
	std::cout<<"Entries generated: "<<data->getNumEvents()<<std::endl;
	foodal->cd();
	foodal->SetLogz(false);
	dalitzpp0_dat_hist.Rebin2D(10,10);
	dalitzpp0_dat_hist.Draw("colz");
	dalitzpp0_dat_hist.SetStats(0);
	//foodal->SaveAs("Dalitz_D2KKK_temp.root");
	foodal->SaveAs("D2KKK_Plots/Dalitz_D2KKK_temp.png");

}

void runToyGeneration(int numFile = 0){
	m12   = Observable("m12",   0.9, 2.0);
	m12.setNumBins(1500);

	m13   = Observable("m13",   0.9, 2.0);
	m13.setNumBins(1500);
	eventNumber = EventNumber("eventNumber", 0, INT_MAX);
	signalDalitz = makeSignalPdf(0,false);
	vector<PdfBase*> comps;
	comps.clear();
	comps.push_back(signalDalitz);

	std::cout << "Creating overall PDF\n";
	ProdPdf* overallSignal = new ProdPdf("overallSignal", comps);
	gettimeofday(&startTime, NULL);
	startCPU = times(&startProc);
	//  makeToyDalitzData (signalDalitz);
	makeToyDalitzData (overallSignal, numFile);
	stopCPU = times(&stopProc);
	gettimeofday(&stopTime, NULL);
}

void getToyData (std::string toyFileName) {
	 TH2F dalitzplot("dalitzplot", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit(), m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
  std::vector<Observable> vars;
  vars.push_back(m12);
  vars.push_back(m13);
  vars.push_back(eventNumber);
  data = new UnbinnedDataSet(vars);
//  const int MAXEVT = 1e4;

  const string suffix = ".root";
  if (toyFileName.rfind(suffix)+suffix.length() == toyFileName.length()){
      std::cout<<"Reading file "<<toyFileName<<std::endl;
      TFile*f = TFile::Open(toyFileName.c_str());
      TTree*t = (TTree*)f->Get("DecayTree");
      //TTree*t = (TTree*)f->Get("newTree");
      std::cout<<"Entries: "<<t->GetEntries()<<std::endl;
      assert(t);
      double m2_12, m2_13;
      //t->SetBranchAddress("s12_KK_DTF", &m2_12);
      //t->SetBranchAddress("s13_KK_DTF", &m2_13);
      t->SetBranchAddress("s12", &m2_12);
      t->SetBranchAddress("s13", &m2_13);
      for (int i=0;i<t->GetEntries()/*&&i<MAXEVT*/;i++){
      //for (int i=0;i<100000;i++){
          t->GetEntry(i);
          m12.setValue(m2_12);
          m13.setValue(m2_13);
          eventNumber.setValue(data->getNumEvents());
          data->addEvent();
          dalitzplot.Fill(m12.getValue(), m13.getValue());
      }
      f->Close();
  }
  else{
  std::ifstream reader;
  reader.open(toyFileName.c_str());
  std::string buffer;
  while (!reader.eof()) {
    reader >> buffer;
    if (buffer == "====") break;
    //std::cout << buffer;
  }

  double dummy = 0;
  while (!reader.eof()) {
    reader >> dummy;
    reader >> dummy;      // m23, m(pi+ pi-), called m12 in processToyRoot convention.
    reader >> m12; // Already swapped according to D* charge. m12 = m(pi+pi0)
    reader >> m13;

    // Errors on Dalitz variables
    reader >> dummy;
    reader >> dummy;
    reader >> dummy;

    reader >> dummy; // Decay time
    reader >> dummy; // sigma_t

    reader >> dummy; // Md0
    reader >> dummy; // deltaM
    reader >> dummy; // ProbSig
    reader >> dummy; // Dst charge
    reader >> dummy; // Run
    reader >> dummy; // Event
    reader >> dummy; // Signal and four bkg fractions.
    reader >> dummy;
    reader >> dummy;
    reader >> dummy;
    reader >> dummy;

    eventNumber.setValue(data->getNumEvents());
    data->addEvent();

    dalitzplot.Fill(m12.getValue(), m13.getValue());
  }}


  dalitzplot.SetStats(0);
  dalitzplot.Draw("colz");
  foodal->SaveAs("dalitzplot_D2KKK_gen.png"); }

/*GooPdf* makeKzeroVeto () {
	if (kzero_veto) return kzero_veto;


	Variable minimum("veto_min",0.475*0.475);
	Variable maximum("veto_max", 0.505*0.505);
	VetoInfo kVetoInfo(minimum,maximum,PAIR_23);

	vector<VetoInfo> vetos; vetos.push_back(kVetoInfo);
        kzero_veto = new DalitzVetoPdf("kzero_veto", m12, m13, motherM, dau1M, dau2M, dau3M, vetos);


	return kzero_veto;
}*/

void createWeightHistogram () {

  TFile*f = TFile::Open("Fit_Input/effspline300.root");
  weightHistogram = (TH2F*)f->Get("eff_spline");
  weightHistogram->SetStats(false);
}

void createBackgroundHistogram () {
  TFile*f = TFile::Open("Fit_Input/bkg_histo_300bins.root");
  bkgHistogram = (TH2F*)f->Get("bkgHist_acc");
  bkgHistogram->SetStats(false);
}

GooPdf* makeEfficiencyPdf () {
  vector<Observable> lvars;
  lvars.push_back(m12);
  lvars.push_back(m13);
  BinnedDataSet* binEffData = new BinnedDataSet(lvars);
  //createWeightHistogram();
  // Now testing your efficiency data by uniformly generating m12,m13 values
  TRandom3 donram(0);
  for (int i = 0; i < NevG; i++){
    do{
    m12.setValue(donram.Uniform(m12.getLowerLimit(), m12.getUpperLimit()));
    m13.setValue(donram.Uniform(m13.getLowerLimit(), m13.getUpperLimit()));
    }while(!cpuDalitz(m12.getValue(), m13.getValue(), MMass , D1Mass, D2Mass,D3Mass));
    //Weight will not be one if the physics boundary crosses the bin square.
    double weight = weightHistogram->GetBinContent(weightHistogram->FindBin(m12.getValue(), m13.getValue()));
    binEffData->addWeightedEvent(weight);
    //if (underlyingBins) underlyingBins->Fill(m12->value, m13->value, weight);
    // Imposing the requirement on efficiency symmetry along m12=m13 line
      if (doEffSwap){
      double swapmass = m12.getValue(); m12.setValue(m13.getValue()); m13.setValue(swapmass);
      weight = weightHistogram->GetBinContent(weightHistogram->FindBin(m12.getValue(), m13.getValue()));
      binEffData->addWeightedEvent(weight);
      //if (underlyingBins) underlyingBins->Fill(m12->value, m13->value, weight);
      //swapmass = m12->value; m12->value = m13->value; m13->value = swapmass;
      }
  }
  if (saveEffPlot) {
    foodal->cd();
    weightHistogram->Draw("colz");
    foodal->SaveAs("plots/efficiency_bins.png");
    foodal->SetLogz(true);
    foodal->SaveAs("plots/efficiency_bins_log.png");
    foo->cd();
  }
 // Smooth
  Variable effSmoothing("effSmoothing", 0);
  SmoothHistogramPdf* ret = new SmoothHistogramPdf("efficiency", binEffData, effSmoothing);
  return ret;
}

GooPdf* makeBackgroundPdf () {
  vector<Observable> lvars;
  lvars.push_back(m12);
  lvars.push_back(m13);
  BinnedDataSet* binBkgData = new BinnedDataSet(lvars);
  createBackgroundHistogram();
  // Now testing your efficiency data by uniformly generating m12,m13 values
  TRandom3 donram(0);
  for (int i = 0; i < NevG; i++){
    do{
    m12.setValue(donram.Uniform(m12.getLowerLimit(), m12.getUpperLimit()));
    m13.setValue(donram.Uniform(m13.getLowerLimit(), m13.getUpperLimit()));
    }while(!cpuDalitz(m12.getValue(), m13.getValue(), MMass , D1Mass, D2Mass,D3Mass));
    //Weight will not be one if the physics boundary crosses the bin square.
    double weight = bkgHistogram->GetBinContent(bkgHistogram->FindBin(m12.getValue(), m13.getValue()));
    binBkgData->addWeightedEvent(weight);
    // Imposing the requirement on efficiency symmetry along m12=m13 line
      if (doEffSwap){
      double swapmass = m12.getValue(); m12.setValue(m13.getValue()); m13.setValue(swapmass);
      weight = bkgHistogram->GetBinContent(bkgHistogram->FindBin(m12.getValue(), m13.getValue()));
      binBkgData->addWeightedEvent(weight);
      }
  }
  if (saveBkgPlot) {
    foodal->cd();
    bkgHistogram->Draw("colz");
    foodal->SetLogz(false);
    foodal->SaveAs("plots/background_bins.png");
    foodal->SetLogz(true);
    foodal->SaveAs("plots/background_bins_log.png");
    foo->cd();
  }
  Variable* effSmoothing = new Variable("effSmoothing",0);
  SmoothHistogramPdf* ret = new SmoothHistogramPdf("efficiency", binBkgData, *effSmoothing);
  return ret;
}

vector<fptype> HH_bin_limits;
vector<Variable> pwa_coefs_amp;
vector<Variable> pwa_coefs_phs;

ResonancePdf* loadPWAResonance(const string fname = pwa_file, bool fixAmp = false,unsigned int cyc=PAIR_12){
  std::ifstream reader;
  reader.open(fname.c_str());
  assert(reader.good());
  HH_bin_limits.clear();
  pwa_coefs_amp.clear();
  pwa_coefs_phs.clear();
  double e1,e2,e3,e4;
  double emag,ephs;
  int i = 0;
  while (reader >> e1 >> e2 >> e3 >> e4) {
      HH_bin_limits.push_back(e1*e1);

      emag = sqrt(e2*e2+e3*e3);
      //emag = e2;
      ephs = TMath::ATan2(e3,e2);
      //ephs = e3;
      sprintf(strbuffer, "pwa_coef_%d_mag", i);
      Variable va(strbuffer, emag, .000001, 0, 10000);//0.9*emag, 1.1*emag);
      sprintf(strbuffer, "pwa_coef_%d_phase", i);
      Variable vp(strbuffer, ephs, .000001, -360, 360);//0.9*ephs, 1.1*ephs);

      pwa_coefs_amp.push_back(va);
      pwa_coefs_phs.push_back(vp);
      i++;
      cout << "s12 = " << e1*e1 << ", mag = " << emag << ", phs = " << (180/PI)*ephs << endl;
  }
  //const fptype scale = 1;
  Variable swave_amp_real("swave_amp_real", 3.0,   0.001, 0, 0);
  Variable swave_amp_imag("swave_amp_imag", 0.0,   0.001, 0, 0);
  swave_amp_real.setFixed(true);
  swave_amp_imag.setFixed(true);

  if (fixAmp) { swave_amp_real.setValue(1.); swave_amp_imag.setValue(0.); swave_amp_real.setFixed(true); swave_amp_imag.setFixed(true); }
  cout<<"Numbers loaded: "<<HH_bin_limits.size()<<" / "<<i<<endl;

  ResonancePdf* swave = new Resonances::Spline("swave", swave_amp_real,swave_amp_imag, HH_bin_limits, pwa_coefs_amp, pwa_coefs_phs,cyc);
  return swave;
}


DalitzPlotPdf* makeSignalPdf (GooPdf* eff,bool fixAmps) {
	DecayInfo3 dtop0pp;
	dtop0pp.motherMass  = MMass;
	dtop0pp.daug1Mass  = D1Mass;
	dtop0pp.daug2Mass  = D2Mass;
	dtop0pp.daug3Mass  = D3Mass;
	dtop0pp.meson_radius  = 1.5;


  /*  // Make a random number generater heres

		random_device rd;

		mt19937 mt(rd());

		normal_distribution<double> rand_gen(0.0,0.1);


	auto rhop  = new Resonances::RBW("rhop",
			Variable("rhop_amp_real", 1),
			Variable("rhop_amp_imag", 0),
			fixedRhoMass,
			fixedRhoWidth,
			1,
			PAIR_12);






    auto var_func = [&](std::string name, double start, double err) -> Variable {
       return fixAmps ?
              Variable(name, start) :
              Variable(name, start + rand_gen(mt), err, 0, 0);
    };*/



  //phi
  Variable phi_amp_real("phi_amp_real", 1);
  Variable phi_amp_imag("phi_amp_imag", 0);
  fixedPhiMass.setFixed(true);
  fixedPhiWidth.setFixed(true);

ResonancePdf* phi  = new Resonances::RBW("phi",phi_amp_real,phi_amp_imag,fixedPhiMass,fixedPhiWidth,1,PAIR_12,true);
ResonancePdf* phi13  = new Resonances::RBW("phi13",phi_amp_real,phi_amp_imag,fixedPhiMass,fixedPhiWidth,1,PAIR_13,true);


  // f0(980)
  Variable f0_amp_real("f0_amp_real",    12.341*cos(-62.852*(PI/180)),   0.0001, -100, 100);
  Variable f0_amp_imag("f0_amp_imag",    12.341*sin(-62.852*(PI/180)),   0.0001, -100, 100);
  Variable f0Mass("f0Mass", 0.965);
  Variable f0g1("f0g1", 0.165);
  Variable rg1og2("rg1og2", 4.21);//,1.0,5.0);

  ResonancePdf* f0  = new Resonances::FLATTE("f0",f0_amp_real,f0_amp_imag,f0Mass,f0g1,rg1og2,PAIR_12, true); //Required to be symmetric
	ResonancePdf* f013  = new Resonances::FLATTE("f013",f0_amp_real,f0_amp_imag,f0Mass,f0g1,rg1og2,PAIR_13, true); //Required to be symmetric

  // f0(X)

  Variable f0X_amp_real("f0X_amp_real",  11.918*cos(20.248*(PI/180)),   0.0001, -100, 100);
  Variable f0X_amp_imag("f0X_amp_imag",  11.918*sin(20.248*(PI/180)),   0.0001, -100, 100);
  Variable f0XMass("f0XMass",    1.41478);//,   0.00001,    1.00, 3.00);
  Variable f0XWidth("f0XWidth",  0.309491);//,   0.00001, 0.00005, 3.00);

  ResonancePdf* f0X  = new Resonances::RBW("f0X",f0X_amp_real,f0X_amp_imag,f0XMass,f0XWidth,(unsigned int)0,PAIR_12, true); //Required to be symmetric
	ResonancePdf* f0X13  = new Resonances::RBW("f0X13",f0X_amp_real,f0X_amp_imag,f0XMass,f0XWidth,(unsigned int)0,PAIR_13, true); //Required to be symmetric

  // NR
  Variable nonr_amp_real("nonr_amp_real", 0.0,   0.001, -100, +100);
  Variable nonr_amp_imag("nonr_amp_imag", 0.0,   0.001, -100, +100);
  ResonancePdf* nonr  = new Resonances::NonRes("nonr",nonr_amp_real,nonr_amp_imag);

  //bool fixAmps = false;
  ResonancePdf* swave = loadPWAResonance(pwa_file, fixAmps,PAIR_12);
  ResonancePdf* swave13 = loadPWAResonance(pwa_file, fixAmps,PAIR_13);

  dtop0pp.resonances.push_back(phi);
	dtop0pp.resonances.push_back(phi13);

	dtop0pp.resonances.push_back(swave);
	dtop0pp.resonances.push_back(swave13);

	//dtop0pp.resonances.push_back(f0X);
	//dtop0pp.resonances.push_back(f0X13);

	//dtop0pp.resonances.push_back(f0);
	//dtop0pp.resonances.push_back(f013);

	dtop0pp.resonances.push_back(nonr);




  if (!eff) {
    // By default create a constant efficiency.
    vector<Variable> offsets;
    vector<Observable> observables;
    vector<Variable> coefficients;

    observables.push_back(m12);
    observables.push_back(m13);
    offsets.push_back(constantZero);
    offsets.push_back(constantZero);
    coefficients.push_back(constantOne);
    eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);
  }
  comps.clear();

  return new DalitzPlotPdf("signalPDF", m12, m13, eventNumber, dtop0pp, eff);
}

void DalitzNorm(GooPdf* overallSignal,int N){


		random_device rd;
		mt19937 mt(rd());
		uniform_real_distribution<double> xyvalues(0.9,2.0);

		std::vector<Observable> vars;
		vars.push_back(m12);
		vars.push_back(m13);
		vars.push_back(eventNumber);

		std::vector<fptype> rpdfValuesvec;

    UnbinnedDataSet data(vars);
		eventNumber = 0;


		for(int i=0; i<N ; i++){

			m12 = xyvalues(mt);
			m13 = xyvalues(mt);



				if(cpuDalitz(m12.getValue(), m13.getValue())==1 ){
        	data.addEvent();

					eventNumber.setValue(eventNumber.getValue()+1);

				}else{

					data.addEvent();
					m12 =0;
					m13 =0;
					eventNumber.setValue(eventNumber.getValue()+1);
				}

		}



		overallSignal->setData(&data);
		signalDalitz->setDataSize(data.getNumEvents());
		std::vector<std::vector<double>> pdfValues = overallSignal->getCompProbsAtDataPoints();


		double buffer = 0;

		cout<< "Sample Size = " << pdfValues[0].size() << endl;

		 for(int k=0; k < pdfValues[0].size();k++){

			 buffer += pdfValues[0][k];

		}

		double  mean = buffer/N;
		double diff = 0;

		for(int l=0;l<pdfValues[0].size();l++){

			diff += (pdfValues[0][l] - mean)*(pdfValues[0][l] - mean);

		}

		double variance = diff/(N-1);
		double sigma = sqrt(variance);
		double RMS = sigma*V/sqrt(N);

		double integral = V*mean;

		std::cout << "Integral: "<< integral << "\t Error: " << RMS << "\n"; ;

}


void runIntegration(int N = 10000){

	  //TApplication* rootapp = new TApplication("rootapp",&argc,argv);

	signalDalitz = makeSignalPdf(0,false);

	std::vector<PdfBase*> comps;
	comps.clear();
	comps.push_back(signalDalitz);

	ProdPdf* overallSignal = new ProdPdf("overallSignal", comps);



	DalitzNorm(overallSignal,N);

	std::cout  << '\n';


}



void drawFitPlotsWithPulls(TH1* hd, TH1* ht, string plotdir){
	const char* hname = hd->GetName();
	char obsname[10];
	for (int i=0;;i++) {
		if (hname[i]=='_') obsname[i] = '\0';
		else obsname[i] = hname[i];
		if (obsname[i] == '\0') break;
	}
	ht->Scale(hd->Integral()/ht->Integral());
	foo->cd();
	foo->Clear();
	ht->Draw("l");
	hd->Draw("epsame");
	sprintf(strbuffer, "%s/%s_fit.png", plotdir.c_str(), obsname);
	foo->SaveAs(strbuffer);
	sprintf(strbuffer, "%s/%s_fit.pdf", plotdir.c_str(), obsname);
	foo->SaveAs(strbuffer);
	/*    sprintf(strbuffer, "%s/%s_fit_log.pdf", plotdir.c_str(), obsname);
		  foo->SaveAs(strbuffer);*/
}

void makeToyDalitzPdfPlots (GooPdf* overallSignal, string plotdir = "plots") {
	TH1F m12_dat_hist("m12_dat_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
	m12_dat_hist.SetStats(false);
	m12_dat_hist.SetMarkerStyle(8);
	m12_dat_hist.SetMarkerSize(1);
	m12_dat_hist.GetXaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV]");
	sprintf(strbuffer, "Events / %.1f MeV", 1e3*m12_dat_hist.GetBinWidth(1));
	m12_dat_hist.GetYaxis()->SetTitle(strbuffer);
	TH1F m12_pdf_hist("m12_pdf_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
	m12_pdf_hist.SetStats(false);
	m12_pdf_hist.SetLineColor(kBlue);
	m12_pdf_hist.SetLineWidth(3);
	TH1F m13_dat_hist("m13_dat_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
	m13_dat_hist.SetStats(false);
	m13_dat_hist.SetMarkerStyle(8);
	m13_dat_hist.SetMarkerSize(1);
	m13_dat_hist.GetXaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV]");
	sprintf(strbuffer, "Events / %.1f MeV", 1e3*m13_dat_hist.GetBinWidth(1));
	m13_dat_hist.GetYaxis()->SetTitle(strbuffer);
	TH1F m13_pdf_hist("m13_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
	m13_pdf_hist.SetStats(false);
	m13_pdf_hist.SetLineColor(kBlue);
	m13_pdf_hist.SetLineWidth(3);
	TH1F m23_dat_hist("m23_dat_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
	m23_dat_hist.SetStats(false);
	m23_dat_hist.SetMarkerStyle(8);
	m23_dat_hist.SetMarkerSize(1.2);
	m23_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{-}) [GeV]");
	sprintf(strbuffer, "Events / %.1f MeV", 1e3*m13_dat_hist.GetBinWidth(1));
	m23_dat_hist.GetYaxis()->SetTitle(strbuffer);
	TH1F m23_pdf_hist("m23_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
	m23_pdf_hist.SetStats(false);
	m23_pdf_hist.SetLineColor(kBlue);
	m23_pdf_hist.SetLineWidth(3);
	double totalPdf = 0;
	double totalDat = 0;
	TH2F dalitzpp0_dat_hist("dalitzpp0_dat_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit(), m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
	dalitzpp0_dat_hist.SetStats(false);
	dalitzpp0_dat_hist.GetXaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV]");
	dalitzpp0_dat_hist.GetYaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV^{2}]");
	TH2F dalitzpp0_pdf_hist("dalitzpp0_pdf_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit(), m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
	/*  dalitzpp0_pdf_hist.GetXaxis()->SetTitle("m^{2}(K^{-} #pi^{0}) [GeV^{2}]");
		dalitzpp0_pdf_hist.GetYaxis()->SetTitle("m^{2}(K^{-} #pi^{+}) [GeV^{2}]");*/
	dalitzpp0_pdf_hist.GetXaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV^{2}]");
	dalitzpp0_pdf_hist.GetYaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV^{2}]");
	dalitzpp0_pdf_hist.SetStats(false);
	std::vector<Observable> vars;
	vars.push_back(m12);
	vars.push_back(m13);
	vars.push_back(eventNumber);
	UnbinnedDataSet currData(vars);
	int evtCounter = 0;

	for (int i = 0; i < m12.getNumBins(); ++i) {
		m12.setValue(m12.getLowerLimit() + (m12.getUpperLimit() - m12.getLowerLimit())*(i + 0.5) / m12.getNumBins());
		for (int j = 0; j < m13.getNumBins(); ++j) {
			m13.setValue(m13.getLowerLimit() + (m13.getUpperLimit() - m13.getLowerLimit())*(j + 0.5) / m13.getNumBins());
			if (!cpuDalitz(m12.getValue(), m13.getValue(), MMass , D1Mass, D2Mass,D3Mass)) continue;
			eventNumber.setValue(evtCounter);
			evtCounter++;
			currData.addEvent();
		}
	}
	overallSignal->setData(&currData);
	signalDalitz->setDataSize(currData.getNumEvents());
	std::vector<std::vector<double> > pdfValues = overallSignal->getCompProbsAtDataPoints();
	for (unsigned int j = 0; j < pdfValues[0].size(); ++j) {
		double currm12 = currData.getValue(m12, j);
		double currm13 = currData.getValue(m13, j);

		dalitzpp0_pdf_hist.Fill(currm12, currm13, pdfValues[0][j]);
		m12_pdf_hist.Fill(currm12, pdfValues[0][j]);
		m13_pdf_hist.Fill(currm13, pdfValues[0][j]);
		m23_pdf_hist.Fill(cpuGetM23(currm12, currm13), pdfValues[0][j]);
		totalPdf     += pdfValues[0][j];
	}
	foodal->cd();
	foodal->SetLogz(false);
	dalitzpp0_pdf_hist.Draw("colz");
    std::string command = "mkdir -p " + plotdir;
    if (system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making plot directory {} failed", plotdir);
	foodal->SaveAs((plotdir + "/dalitzpp0_pdf.png").c_str());
	/*  m12_pdf_hist.Draw("");
		foodal->SaveAs((plotdir + "/m12_pdf_hist.png").c_str());
		m13_pdf_hist.Draw("");
		foodal->SaveAs((plotdir + "/m13_pdf_hist.png").c_str());
		if (!data) return;*/
	for (unsigned int evt = 0; evt < data->getNumEvents(); ++evt) {
		double data_m12 = data->getValue(m12, evt);
		m12_dat_hist.Fill(data_m12);
		double data_m13 = data->getValue(m13, evt);
		m13_dat_hist.Fill(data_m13);
		dalitzpp0_dat_hist.Fill(data_m12, data_m13);
		m23_dat_hist.Fill(cpuGetM23(data_m12, data_m13));
		totalDat++;
	}
	dalitzpp0_dat_hist.Draw("colz");
	foodal->SaveAs((plotdir + "/dalitzpp0_dat.png").c_str());

	drawFitPlotsWithPulls(&m12_dat_hist, &m12_pdf_hist, plotdir);
	drawFitPlotsWithPulls(&m13_dat_hist, &m13_pdf_hist, plotdir);
	drawFitPlotsWithPulls(&m23_dat_hist, &m23_pdf_hist, plotdir);
}

void runToyFit (std::string toyFileName) {
	m12 = Observable("m12", 0.9, 2.0);
	m13 = Observable("m13", 0.9, 2.0);
	m12.setNumBins(nbins);
	m13.setNumBins(nbins);
	eventNumber = EventNumber("eventNumber", 0, INT_MAX);
	getToyData(toyFileName);


	signalDalitz = makeSignalPdf();
	comps.clear();
	comps.push_back(signalDalitz);
	ProdPdf* overallSignal = new ProdPdf("overallSignal", comps);
	overallSignal->setData(data);
	signalDalitz->setDataSize(data->getNumEvents());

	FitManager datapdf(overallSignal);

	for(int i=0;i<HH_bin_limits.size();i++){
      pwa_coefs_amp[i].setFixed(false);
      pwa_coefs_phs[i].setFixed(false);
      //pwa_coefs_amp[i]->error = pwa_coefs_phs[i]->error = 1.0;
  }

	gettimeofday(&startTime, NULL);
	startCPU = times(&startProc);
    datapdf.setVerbosity(verbosity);

		 // Maybe make optional? With a command line switch?
	datapdf.fit();
	stopCPU = times(&stopProc);
	gettimeofday(&stopTime, NULL);

	makeToyDalitzPdfPlots(overallSignal);
}

int main (int argc, char** argv) {

    GooFit::Application app{"D2K3_toy", argc, argv};
    app.add_option("-v,--verbose", verbosity, "Set the verbosity (to 0 for example", true);

    int fit_value;
    std::string name = "dalitz_mytoyMC_000.txt";

    auto fit = app.add_subcommand("fit");
    fit->add_option("-i,--int", fit_value, "A number to load");
    auto name_opt = fit->add_option("-n,--name,name", name, "The filename to load", true)
        ->excludes("--int");


		int N;
		auto run = app.add_subcommand("run");
		run->add_option("N",N, "")
		   ->required();


    int value;
    auto gen = app.add_subcommand("gen");
    gen->add_option("value", value, "The number to generate")
        ->required();

    app.require_subcommand(1);

    GOOFIT_PARSE(app);

    if(name_opt->count())
        name = fmt::format("dalitz_mytoyMC_{0:3}.txt", fit_value);

	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasColor(10);
	gStyle->SetFrameFillColor(10);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetPadColor(0);
	gStyle->SetTitleColor(1);
	gStyle->SetStatColor(0);
	gStyle->SetFillColor(0);
	gStyle->SetFuncWidth(1);
	gStyle->SetLineWidth(1);
	gStyle->SetLineColor(1);
	gStyle->SetPalette(kViridis, 0);
	gStyle->SetNumberContours(512);
	gStyle->SetOptStat("RMe");
	foo = new TCanvas();
	foodal = new TCanvas();
	foodal->Size(10, 10);


    if(*fit)
	    runToyFit(name);
    if(*gen)
	    runToyGeneration(value);
		if(*run) {
				CLI::AutoTimer timer("Integration");
				runIntegration(N);
		}

	// Print total minimization time
	double myCPU = stopCPU - startCPU;
	double totalCPU = myCPU;

	timersub(&stopTime, &startTime, &totalTime);
	std::cout << "Wallclock time  : " << totalTime.tv_sec + totalTime.tv_usec/1000000.0 << " seconds." << std::endl;
	std::cout << "CPU time: " << (myCPU / CLOCKS_PER_SEC) << std::endl;
	std::cout << "Total CPU time: " << (totalCPU / CLOCKS_PER_SEC) << std::endl;
	myCPU = stopProc.tms_utime - startProc.tms_utime;
	std::cout << "Processor time: " << (myCPU / CLOCKS_PER_SEC) << std::endl;

	return 0;
}
