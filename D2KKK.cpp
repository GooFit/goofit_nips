// ROOT stuff
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TROOT.h>
#include <TMinuit.h>


// System stuff
#include <CLI/Timer.hpp>
#include <fstream>



// GooFit stuff
#include <goofit/Application.h>
#include <goofit/BinnedDataSet.h>
#include <goofit/FitManager.h>
#include <goofit/fitting/FitManagerMinuit2.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/PDFs/basic/SmoothHistogramPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/PDFs/physics/DalitzPlotPdf.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>

#include <goofit/detail/Style.h>
#include <goofit/PDFs/physics/DalitzPlotter.h>

#include <Minuit2/MnScan.h>

using namespace std;
using namespace GooFit;
using namespace ROOT;

UnbinnedDataSet *data    = nullptr;
const unsigned int nbins = 1000;
TH2F *weightHistogram    = nullptr;
TH2F *bkgHistogram       = nullptr;
TH2F *underlyingBins     = nullptr;

// How many events will be generated for Eff Bkg?
const double NevG = 1e7;

// PWA INPUT FILE NAME
const string pwa_file = "files/PWA_COEFFS_50.txt";

// FIT OR JUST PLOT?
bool fit = true;

Observable m12("m12", 0.9, 2.0);
Observable m13("m13", 0.9, 2.0);
EventNumber eventNumber("eventNumber");
bool fitMasses = false;



const fptype _mDp      = 1.86962;
const fptype KPlusMass = 0.493677;
double V = (m12.getUpperLimit() - m12.getLowerLimit()) * (m13.getUpperLimit() - m13.getLowerLimit()); // Volume

const fptype D1Mass  = KPlusMass; //
const fptype D2Mass  = KPlusMass;
const fptype D3Mass  = KPlusMass;
const fptype D1Mass2 = D1Mass * D1Mass;
const fptype D2Mass2 = D2Mass * D2Mass;
const fptype D3Mass2 = D3Mass * D3Mass;
const fptype MMass   = _mDp;
const fptype MMass2  = MMass * MMass;



// Constants used in more than one PDF component.
Variable motherM("motherM", MMass);
Variable massSum("massSum", MMass2 + D1Mass2 + D2Mass2 + D3Mass2); // = 3.53481
Variable constantOne("constantOne", 1);
Variable constantZero("constantZero", 0);

std::vector<PdfBase *> comps;

// I don't like Globals! Henry
int verbosity = 3;

GooPdf *kzero_veto = nullptr;
double mesonRad = 1.5;
DalitzPlotPdf *signalDalitz;
bool doEffSwap   = true;
bool saveEffPlot = true;
bool saveBkgPlot = true;


// Declarations for (some of) the functions

DalitzPlotPdf *makeSignalPdf(GooPdf *eff = 0, bool fixAmps = false);
fptype cpuGetM23(fptype massPZ, fptype massPM) { return (massSum.getValue() - massPZ - massPM); }


void makeToyDalitzData(GooPdf *overallSignal, std::string name, size_t nTotal) {
    
    DalitzPlotter dp(overallSignal, signalDalitz);


    // Generate data
    data = new UnbinnedDataSet({m12, m13, eventNumber});
    
    { // Plotting block
        TCanvas foo;
        auto th1 = dp.make2D();
        th1->Rebin2D(5,5);
        th1->Draw("COLZ");
        foo.SaveAs("plots/plot1.png");
    }
        
    dp.fillDataSetMC(*data, nTotal);
    
    TH2F dalitzpp0_dat_hist("dalitzpp0_dat_hist",
                            "",
                            200,
                            m12.getLowerLimit(),
                            m12.getUpperLimit(),
                            200,
                            m13.getLowerLimit(),
                            m13.getUpperLimit());
    dalitzpp0_dat_hist.GetXaxis()->SetTitle("m^{2}(K^{-}K^{+}) [GeV^{2}]");
    dalitzpp0_dat_hist.GetYaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV^{2}]");
    
    // Make a writer
    {
        ofstream writer(name);
        
        // Fill histogram with generated data
        for(size_t i=0; i<data->getNumEvents(); i++) {
            data->loadEvent(i);
            dalitzpp0_dat_hist.Fill(m12, m13);
            
            writer << i << '\t' << m12.getValue() << '\t' << m13.getValue() << '\n';
        }
    }
    
   
    
    std::cout << "Entries generated: " << data->getNumEvents() << std::endl;
    TCanvas foo;
    foo.SetLogz(false);
    dalitzpp0_dat_hist.Draw("colz");
    dalitzpp0_dat_hist.SetStats(0);
    foo.SaveAs("plots/Dalitz_D2KKK_temp.png");
}


void runToyGeneration(std::string name, size_t events) {
    m12.setNumBins(1500);
    m13.setNumBins(1500);

    signalDalitz = makeSignalPdf(0, false);

    std::cout << "Creating overall PDF\n";
    ProdPdf *overallSignal = new ProdPdf("overallSignal", {signalDalitz});
    
    {
        CLI::AutoTimer timer{"makeToyDalitzData"};
        makeToyDalitzData(overallSignal, name, events);
    }

}

void getToyData(std::string toyFileName) {
    TH2F dalitzplot("dalitzplot",
                    "",
                    m12.getNumBins(),
                    m12.getLowerLimit(),
                    m12.getUpperLimit(),
                    m13.getNumBins(),
                    m13.getLowerLimit(),
                    m13.getUpperLimit());
    
    data = new UnbinnedDataSet({m12, m13, eventNumber});

    const string suffix = ".root";
    
    if(toyFileName.rfind(suffix) + suffix.length() == toyFileName.length()) {
        GOOFIT_INFO("Reading ROOT file: {}", toyFileName);
        
        TFile *f = TFile::Open(toyFileName.c_str());
        TTree *t = (TTree *)f->Get("DecayTree");

        std::cout << "Entries: " << t->GetEntries() << std::endl;
        assert(t);
        double m2_12, m2_13;

        t->SetBranchAddress("s12", &m2_12);
        t->SetBranchAddress("s13", &m2_13);
        for(int i = 0; i < t->GetEntries(); i++) {

            t->GetEntry(i);
            m12.setValue(m2_12);
            m13.setValue(m2_13);
            eventNumber.setValue(data->getNumEvents());
            data->addEvent();
            dalitzplot.Fill(m12.getValue(), m13.getValue());
        }
        f->Close();
    } else {
        GOOFIT_INFO("Reading 3 column TEXT file: {}", toyFileName);
        std::ifstream reader;
        reader.open(toyFileName.c_str());
        std::string buffer;

        size_t n = 0;

        while(!reader.eof()) {
            reader >> eventNumber;
            //eventNumber = n;
            reader >> m12;
            reader >> m13;
            data->addEvent();

            dalitzplot.Fill(m12.getValue(), m13.getValue());

            n++;
        }
		reader.close();
    }

    TCanvas foo;
    dalitzplot.SetStats(0);
    dalitzplot.Draw("colz");
    foo.SaveAs("plots/dalitzplot_D2KKK_gen.png");
}


void createWeightHistogram() {
    TFile *f        = TFile::Open("Fit_Input/effspline300.root");
    weightHistogram = (TH2F *)f->Get("eff_spline");
    weightHistogram->SetStats(false);
}

void createBackgroundHistogram() {
    TFile *f     = TFile::Open("Fit_Input/bkg_histo_300bins.root");
    bkgHistogram = (TH2F *)f->Get("bkgHist_acc");
    bkgHistogram->SetStats(false);
}

GooPdf *makeEfficiencyPdf() {
    vector<Observable> lvars;
    lvars.push_back(m12);
    lvars.push_back(m13);
    BinnedDataSet *binEffData = new BinnedDataSet(lvars);
    // createWeightHistogram();
    // Now testing your efficiency data by uniformly generating m12,m13 values
    TRandom3 donram(0);
    for(int i = 0; i < NevG; i++) {
        do {
            m12.setValue(donram.Uniform(m12.getLowerLimit(), m12.getUpperLimit()));
            m13.setValue(donram.Uniform(m13.getLowerLimit(), m13.getUpperLimit()));
        } while(!inDalitz(m12.getValue(), m13.getValue(), MMass, D1Mass, D2Mass, D3Mass));
        // Weight will not be one if the physics boundary crosses the bin square.
        double weight = weightHistogram->GetBinContent(weightHistogram->FindBin(m12.getValue(), m13.getValue()));
        binEffData->addWeightedEvent(weight);
        // if (underlyingBins) underlyingBins->Fill(m12->value, m13->value, weight);
        // Imposing the requirement on efficiency symmetry along m12=m13 line
        if(doEffSwap) {
            double swapmass = m12.getValue();
            m12.setValue(m13.getValue());
            m13.setValue(swapmass);
            weight = weightHistogram->GetBinContent(weightHistogram->FindBin(m12.getValue(), m13.getValue()));
            binEffData->addWeightedEvent(weight);
            // if (underlyingBins) underlyingBins->Fill(m12->value, m13->value, weight);
            // swapmass = m12->value; m12->value = m13->value; m13->value = swapmass;
        }
    }
    if(saveEffPlot) {
        TCanvas foo;
        foo.cd();
        weightHistogram->Draw("colz");
        foo.SaveAs("plots/efficiency_bins.png");
        foo.SetLogz(true);
        foo.SaveAs("plots/efficiency_bins_log.png");
    }
    // Smooth
    Variable effSmoothing("effSmoothing", 0);
    SmoothHistogramPdf *ret = new SmoothHistogramPdf("efficiency", binEffData, effSmoothing);
    return ret;
}

GooPdf *makeBackgroundPdf() {
    BinnedDataSet *binBkgData = new BinnedDataSet({m12, m13});
    createBackgroundHistogram();
    // Now testing your efficiency data by uniformly generating m12,m13 values
    TRandom3 donram(0);
    for(int i = 0; i < NevG; i++) {
        do {
            m12.setValue(donram.Uniform(m12.getLowerLimit(), m12.getUpperLimit()));
            m13.setValue(donram.Uniform(m13.getLowerLimit(), m13.getUpperLimit()));
        } while(!inDalitz(m12.getValue(), m13.getValue(), MMass, D1Mass, D2Mass, D3Mass));
        // Weight will not be one if the physics boundary crosses the bin square.
        double weight = bkgHistogram->GetBinContent(bkgHistogram->FindBin(m12.getValue(), m13.getValue()));
        binBkgData->addWeightedEvent(weight);
        // Imposing the requirement on efficiency symmetry along m12=m13 line
        if(doEffSwap) {
            double swapmass = m12.getValue();
            m12.setValue(m13.getValue());
            m13.setValue(swapmass);
            weight = bkgHistogram->GetBinContent(bkgHistogram->FindBin(m12.getValue(), m13.getValue()));
            binBkgData->addWeightedEvent(weight);
        }
    }
    if(saveBkgPlot) {
        TCanvas foo;
        bkgHistogram->Draw("colz");
        foo.SetLogz(false);
        foo.SaveAs("plots/background_bins.png");
        foo.SetLogz(true);
        foo.SaveAs("plots/background_bins_log.png");
    }
    Variable *effSmoothing  = new Variable("effSmoothing", 0);
    SmoothHistogramPdf *ret = new SmoothHistogramPdf("efficiency", binBkgData, *effSmoothing);
    return ret;
}

vector<fptype> HH_bin_limits;
vector<Variable> pwa_coefs_amp;
vector<Variable> pwa_coefs_phs;


ResonancePdf *loadPWAResonance(const string fname = pwa_file, bool fixAmp = false) {
    std::ifstream reader;
	GOOFIT_INFO("LOADING FILE {}",fname);
    reader.open(fname.c_str());
    assert(reader.good());
    HH_bin_limits.clear();
    pwa_coefs_amp.clear();
    pwa_coefs_phs.clear();
    double e1, e2, e3, e4;
    double emag, ephs;
    int i = 0;
    while(reader >> e1 >> e2 >> e3 >> e4) {
        HH_bin_limits.push_back(e1 * e1);

        emag = sqrt(e2 * e2 + e3 * e3);
        // emag = e2;
        ephs = TMath::ATan2(e3, e2);
        // ephs = e3;

        Variable va(fmt::format("pwa_coef_{}_mag", i), emag, .000001, 0, 10000); // 0.9*emag, 1.1*emag);
        Variable vp(fmt::format("pwa_coef_{}_phase", i), ephs, .000001, -360, 360); // 0.9*ephs, 1.1*ephs);

        pwa_coefs_amp.push_back(va);
        pwa_coefs_phs.push_back(vp);
        i++;
        cout << "s12 = " << e1 * e1 << ", mag = " << emag << ", phs = " << (180 / M_PI) * ephs << endl;
    }
    // const fptype scale = 1;
    Variable swave_amp_real("swave_amp_real", 3.0, 0.001, 0, 0);
    Variable swave_amp_imag("swave_amp_imag", 0.0, 0.001, 0, 0);
    swave_amp_real.setFixed(true);
    swave_amp_imag.setFixed(true);

    if(fixAmp) {
        swave_amp_real.setValue(1.);
        swave_amp_imag.setValue(0.);
        swave_amp_real.setFixed(true);
        swave_amp_imag.setFixed(true);
    }
    cout << "Numbers loaded: " << HH_bin_limits.size() << " / " << i << endl;

    ResonancePdf *swave_12 = new Resonances::Spline(
        "swave_12", swave_amp_real, swave_amp_imag, HH_bin_limits, pwa_coefs_amp, pwa_coefs_phs, PAIR_12, true);

    return swave_12;
}

DalitzPlotPdf *makeSignalPdf(GooPdf *eff, bool fixAmps) {
    DecayInfo3 dtop0pp;
    dtop0pp.motherMass   = MMass;
    dtop0pp.daug1Mass    = D1Mass;
    dtop0pp.daug2Mass    = D2Mass;
    dtop0pp.daug3Mass    = D3Mass;
    dtop0pp.meson_radius = 1.5;

    // phi

	Variable fixedPhiMass("phi_mass", 1.019461, 0.01, 0.7, 1.8);
	Variable fixedPhiWidth("phi_width", 0.004266, 0.001, 1e-5, 1e-1);
    Variable phi_amp_real("phi_amp_real", 1);
    Variable phi_amp_imag("phi_amp_imag", 0);
    fixedPhiMass.setFixed(true);
    fixedPhiWidth.setFixed(true);

    ResonancePdf *phi
        = new Resonances::RBW("phi", phi_amp_real, phi_amp_imag, fixedPhiMass, fixedPhiWidth, 1, PAIR_12, true);

    // f0

	Variable f0_amp_real("f0_amp_real", 12.341 * cos(-62.852 * (M_PI / 180)), 0.0001, -100, 100);
	Variable f0_amp_imag("f0_amp_imag", 12.341 * sin(-62.852 * (M_PI / 180)), 0.0001, -100, 100);
    Variable f0Mass("f0Mass", 0.965);
    Variable f0g1("f0g1", 0.165);
    Variable rg1og2("rg1og2", 4.21*0.165);

    ResonancePdf *f0
            = new Resonances::FLATTE("f0", f0_amp_real, f0_amp_imag, f0Mass, f0g1, rg1og2, PAIR_12, true); // Required to be symmetric

    // f0(X)

    Variable f0X_amp_real("f0X_amp_real", 11.918 * cos(20.248 * (M_PI / 180)), 0.0001, -100, 100);
    Variable f0X_amp_imag("f0X_amp_imag", 11.918 * sin(20.248 * (M_PI / 180)), 0.0001, -100, 100);
    Variable f0XMass("f0XMass", 1.41478);
    Variable f0XWidth("f0XWidth", 0.309491);
    ;
    ResonancePdf *f0X = new Resonances::RBW("f0X",
                                            f0X_amp_real,
                                            f0X_amp_imag,
                                            f0XMass,
                                            f0XWidth,
                                            (unsigned int)0,
                                            PAIR_12,
                                            true); // Required to be symmetric

    // NR
    Variable nonr_amp_real("nonr_amp_real", 1.0, 0.001, -100, +100);
    Variable nonr_amp_imag("nonr_amp_imag", 0.0, 0.001, -100, +100);
    ResonancePdf *nonr = new Resonances::NonRes("nonr", nonr_amp_real, nonr_amp_imag);

    ResonancePdf *swave_12 = loadPWAResonance(pwa_file, fixAmps);

    dtop0pp.resonances.push_back(phi);
    //dtop0pp.resonances.push_back(swave_12);
    dtop0pp.resonances.push_back(f0X);
    dtop0pp.resonances.push_back(f0);
	dtop0pp.resonances.push_back(nonr);

    if(!eff) {
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

std::tuple<double, double> DalitzNorm(GooPdf *overallSignal, int N) {
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> xyvalues(0.9, 2.0);

    std::vector<Observable> vars;
    vars.push_back(m12);
    vars.push_back(m13);
    vars.push_back(eventNumber);

    std::vector<fptype> rpdfValuesvec;

    UnbinnedDataSet data(vars);
    eventNumber = 0;

    for(int i = 0; i < N; i++) {
        m12 = xyvalues(mt);
        m13 = xyvalues(mt);

        if(inDalitz(m12.getValue(), m13.getValue(), MMass, D1Mass, D2Mass, D3Mass) == 1) {
            eventNumber.setValue(eventNumber.getValue() + 1);
            data.addEvent();
        }
    }

    //overallSignal->normalize();
    overallSignal->setData(&data);
    signalDalitz->setDataSize(data.getNumEvents());
    std::vector<std::vector<double>> pdfValues = overallSignal->getCompProbsAtDataPoints();

    double buffer = 0;

    ofstream wr("res.txt");
    // cout<< "Sample Size = " << pdfValues[0].size() << endl;
    for(int k = 0; k < pdfValues[0].size(); k++) {
        buffer += pdfValues[0][k];

        wr <<  pdfValues[0][k] << '\t' << pdfValues[1][k] << '\n';
    }

    wr.close();

    double mean = buffer / N;
    double diff = 0;

    for(int l = 0; l < pdfValues[0].size(); l++) {
        diff += (pdfValues[0][l] - mean) * (pdfValues[0][l] - mean);
    }

    double variance_f = diff / (N - 1);
    double variance   = V * V * variance_f / N;
    double sigma      = sqrt(variance);
    double integral   = V * mean;

    return std::make_tuple(integral, sigma);
}



void runIntegration(int N = 10000,int Nint=100) {
    signalDalitz = makeSignalPdf(0, false);
    std::vector<PdfBase *> comps;
    comps.clear();
    comps.push_back(signalDalitz);

    ProdPdf *overallSignal = new ProdPdf("overallSignal", comps);

    double arr[Nint];
    double arr_error[Nint];
    std::fill(arr, arr + Nint, 0);
    std::fill(arr_error, arr_error + Nint, 0);

    TH1D integral_hist("integral", "integral", 30, arr[0] * (0.95), arr[Nint-1] * (1.05));

    for(int i = 0; i < Nint; i++) {
        auto integral2 = DalitzNorm(overallSignal, N);
        arr[i]         = std::get<0>(integral2);
        arr_error[i]   = std::get<1>(integral2);
    }

    std::sort(arr, arr + Nint);

    for(int l = 0; l < Nint; l++) {
        integral_hist.Fill(arr[l]);
    }

    integral_hist.GetXaxis()->SetTitle("Integral");
    integral_hist.GetYaxis()->SetTitle("Frequency");

    TCanvas integral_Canvas("integral", "integral", 800, 800);
    integral_hist.Draw("E");
    integral_Canvas.SaveAs("plots/D2KKK_Plots_Integral.png");

    std::cout << '\n';
    std::cout << "<E>_{N_Integrations=100}: " << integral_hist.GetMean() << '\t'
              << "stdError: " << integral_hist.GetMeanError() << "\t\t"
              << "stdDev: " << integral_hist.GetStdDev() << "\n\n";

    double integral, sigma2;
    std::tie(integral, sigma2) = DalitzNorm(overallSignal, N);

    std::cout << "<E>_{N_Integrations=1}: " << integral << '\t' << "<delta_E>: " << sigma2 << "\n\n";
    std::cout << "|stdDev - <delta_E>|= " << abs(integral_hist.GetStdDev() - sigma2) << "\n\n";

	//delete integral_Canvas;
	//delete integral_hist;
}

void drawFitPlotsWithPulls(TH1 *hd, TH1 *ht, string plotdir) {
    const char *hname = hd->GetName();
    char obsname[10];
    for(int i = 0;; i++) {
        if(hname[i] == '_')
            obsname[i] = '\0';
        else
            obsname[i] = hname[i];
        if(obsname[i] == '\0')
            break;
    }
    ht->Scale(hd->Integral() / ht->Integral()*5);
	ht->SetLineColor(kRed);
    ht->SetMarkerStyle(0);

	hd->SetMarkerColor(kBlack);
	hd->Rebin(5);


    TCanvas foo;
    
	hd->Draw("E");
    ht->Draw("HIST C same");
    
    
    foo.SaveAs(TString::Format("plots/%s_fit.png",obsname));
    //foo.SaveAs(TString::Format("plots/%s_fit.pdf",obsname));

}

void makeToyDalitzPdfPlots(GooPdf *overallSignal, string plotdir = "plots") {
    TH1F m12_dat_hist("m12_dat_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    m12_dat_hist.GetXaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV]");
    m12_dat_hist.GetYaxis()->SetTitle(TString::Format("Events / %.1f MeV", 1e3 * m12_dat_hist.GetBinWidth(1)));
    
    TH1F m12_pdf_hist("m12_pdf_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    
    TH1F m13_dat_hist("m13_dat_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m13_dat_hist.GetXaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV]");
    m13_dat_hist.GetYaxis()->SetTitle(TString::Format("Events / %.1f MeV", 1e3 * m13_dat_hist.GetBinWidth(1)));
    
    TH1F m13_pdf_hist("m13_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    
    TH1F m23_dat_hist("m23_dat_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m23_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{-}) [GeV]");
    m23_dat_hist.GetYaxis()->SetTitle(TString::Format("Events / %.1f MeV", 1e3 * m13_dat_hist.GetBinWidth(1)));
    
    TH1F m23_pdf_hist("m23_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    
    double totalPdf = 0;
    double totalDat = 0;
    TH2F dalitzpp0_dat_hist("dalitzpp0_dat_hist",
                            "",
                            m12.getNumBins(),
                            m12.getLowerLimit(),
                            m12.getUpperLimit(),
                            m13.getNumBins(),
                            m13.getLowerLimit(),
                            m13.getUpperLimit());
    dalitzpp0_dat_hist.SetStats(false);
    dalitzpp0_dat_hist.GetXaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV]");
    dalitzpp0_dat_hist.GetYaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV^{2}]");
    TH2F dalitzpp0_pdf_hist("dalitzpp0_pdf_hist",
                            "",
                            m12.getNumBins(),
                            m12.getLowerLimit(),
                            m12.getUpperLimit(),
                            m13.getNumBins(),
                            m13.getLowerLimit(),
                            m13.getUpperLimit());

    dalitzpp0_pdf_hist.GetXaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV^{2}]");
    dalitzpp0_pdf_hist.GetYaxis()->SetTitle("m^{2}(K^{-} K^{+}) [GeV^{2}]");
    dalitzpp0_pdf_hist.SetStats(false);
    std::vector<Observable> vars;
    vars.push_back(m12);
    vars.push_back(m13);
    vars.push_back(eventNumber);
    UnbinnedDataSet currData(vars);
    int evtCounter = 0;

    for(int i = 0; i < m12.getNumBins(); ++i) {
        m12.setValue(m12.getLowerLimit() + (m12.getUpperLimit() - m12.getLowerLimit()) * (i + 0.5) / m12.getNumBins());
        for(int j = 0; j < m13.getNumBins(); ++j) {
            m13.setValue(m13.getLowerLimit()
                         + (m13.getUpperLimit() - m13.getLowerLimit()) * (j + 0.5) / m13.getNumBins());
            if(!inDalitz(m12.getValue(), m13.getValue(), MMass, D1Mass, D2Mass, D3Mass))
                continue;
            eventNumber.setValue(evtCounter);
            evtCounter++;
            currData.addEvent();
        }
    }
    overallSignal->setData(&currData);
    signalDalitz->setDataSize(currData.getNumEvents());
    std::vector<std::vector<double>> pdfValues = overallSignal->getCompProbsAtDataPoints();
    for(unsigned int j = 0; j < pdfValues[0].size(); ++j) {
        double currm12 = currData.getValue(m12, j);
        double currm13 = currData.getValue(m13, j);

        dalitzpp0_pdf_hist.Fill(currm12, currm13, pdfValues[0][j]);
        m12_pdf_hist.Fill(currm12, pdfValues[0][j]);
        m13_pdf_hist.Fill(currm13, pdfValues[0][j]);
        m23_pdf_hist.Fill(cpuGetM23(currm12, currm13), pdfValues[0][j]);
        totalPdf += pdfValues[0][j];
    }
    
    TCanvas foo;
    foo.SetLogz(false);
    dalitzpp0_pdf_hist.Draw("colz");

    foo.SaveAs("plots/dalitzpp0_pdf.png");

    for(unsigned int evt = 0; evt < data->getNumEvents(); ++evt) {
        double data_m12 = data->getValue(m12, evt);
        m12_dat_hist.Fill(data_m12);
        double data_m13 = data->getValue(m13, evt);
        m13_dat_hist.Fill(data_m13);
        dalitzpp0_dat_hist.Fill(data_m12, data_m13);
        m23_dat_hist.Fill(cpuGetM23(data_m12, data_m13));
        totalDat++;
    }
    dalitzpp0_dat_hist.Draw("colz");
    foo.SaveAs("plots/dalitzpp0_dat.png");

    drawFitPlotsWithPulls(&m12_dat_hist, &m12_pdf_hist, plotdir);
    drawFitPlotsWithPulls(&m13_dat_hist, &m13_pdf_hist, plotdir);
    drawFitPlotsWithPulls(&m23_dat_hist, &m23_pdf_hist, plotdir);
}

void fit_fractions(DalitzPlotPdf* signal, std::vector<fpcomplex> coefs){

    size_t n_res = signal->getDecayInfo().resonances.size();

    size_t nEntries = signal->getCachedWave(0).size();

    fptype cache_results[n_res][n_res];

    fptype buffer_all = 0;
    fptype buffer =0 ;
    double norm = 0;

    for(size_t i = 0; i < n_res ; i++){

        for(size_t j = 0; j < n_res ; j++){

            norm =0;
            buffer =0;


            auto cached_i = signal->getCachedWave(i);
            auto cached_j = signal->getCachedWave(j);
            fpcomplex coef_i = coefs[i];
            fpcomplex coef_j = coefs[j];

            for(size_t k = 0 ; k < nEntries ; k++){


                fpcomplex cached_i_val = cached_i[k];
                fpcomplex cached_j_val = cached_j[k];



                if(i!=j) {
                    norm = (coef_i * conj<fptype>(coef_j) * cached_i_val * conj<fptype>(cached_j_val)).real();
                }
                else {
                    fpcomplex multiply = coef_i * conj<fptype>(coef_j) * cached_i_val * conj<fptype>(cached_j_val);
                    norm = thrust::norm(multiply);
                }

                buffer += norm;


            }

            buffer_all += buffer;

            cache_results[i][j] = (buffer/nEntries);

        }



    }

    double den = buffer_all/nEntries;


    double sum = 0;

    std::cout << "nEntries= " << nEntries << '\n';
    for(size_t i = 0; i < n_res ; i++){

        for(size_t j = 0; j< n_res ; j++){
            std::cout << "FF[" << i << "," << j <<"]= " << cache_results[i][j]/den << std::endl;

        }

        sum+=cache_results[i][i]/den;
    }

    std::cout << "Sum[i,i]= " << sum << std::endl;
}

void runToyFit(std::string toyFileName) {
    m12.setNumBins(nbins);
    m13.setNumBins(nbins);
    getToyData(toyFileName);

    GOOFIT_INFO("Number of events in dataset: {}", data->getNumEvents());


    signalDalitz = makeSignalPdf();
    comps.clear();
    comps.push_back(signalDalitz);
    ProdPdf *overallSignal = new ProdPdf("overallSignal", comps);
    overallSignal->setData(data);
    signalDalitz->setDataSize(data->getNumEvents());

    FitManagerMinuit2 fitter(overallSignal);

    fitter.setVerbosity(verbosity);


    for(int i = 0; i < HH_bin_limits.size(); i++) {
        pwa_coefs_amp[i].setFixed(false);
        pwa_coefs_phs[i].setFixed(false);
        // pwa_coefs_amp[i]->error = pwa_coefs_phs[i]->error = 1.0;
    }

    // Maybe make optional? With a command line switch?
    auto func_min = fitter.fit();
    makeToyDalitzPdfPlots(overallSignal);

    auto param = fitter.getParams()->Parameters();

    fpcomplex phi_coef(param[0].Value()*cos(param[1].Value()), param[0].Value()*sin(param[1].Value()));
    fpcomplex f0X_coef(param[2].Value()*cos(param[3].Value()), param[2].Value()*sin(param[3].Value()));
    fpcomplex f0_coef(param[4].Value()*cos(param[5].Value()), param[4].Value()*sin(param[5].Value()));
    fpcomplex NR_coef(param[6].Value()*cos(param[7].Value()), param[6].Value()*sin(param[7].Value()));
    std::vector<fpcomplex> coefs = {phi_coef,f0X_coef,f0_coef,NR_coef};

    fit_fractions(signalDalitz,coefs);




}



int main(int argc, char **argv) {
    int fit_value = 0;
    std::string name;
    
    GooFit::Application app{"D2K3_toy", argc, argv};
    app.add_option("-v,--verbose", verbosity, "Set the verbosity (to 0 for example", true);
    app.add_option("-i,--int", fit_value, "A number to load", true);

    app.add_option("-n,--name", name, "The filename to load")->excludes("--int");
    
    

    auto fit = app.add_subcommand("fit");
    
    int N,e;
    auto run = app.add_subcommand("run");
    run->add_option("-e", e, "sample size")->required();
	run->add_option("-N", N, "Number of integrations")->required();

    size_t events = 1000000;

    auto gen = app.add_subcommand("gen");
    gen->add_option("-e,--events", events, "The number of events to generate", true);

    /// Must get 1 or more subcommands
    app.require_subcommand();

    GOOFIT_PARSE(app);

    if(name.empty())
        name = fmt::format("dalitz_mytoyMC_{0:03}.txt", fit_value);
    
    /// Make the plot directory if it does not exist
    std::string command = "mkdir -p plots";
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making `plots` directory failed");

    GooFit::setROOTStyle();

    if(*gen)
        runToyGeneration(name, events);
    if(*fit)
        runToyFit(name);
    if(*run) {
        CLI::AutoTimer timer("Integration");
        runIntegration(e,N);
        std::cout << "\n\n";
    }

 
    return 0;
}
