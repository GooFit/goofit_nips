// ROOT stuff
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TText.h>
#include <TTree.h>

// System stuff
#include <CLI/Timer.hpp>
#include <fstream>
#include <iostream>
#include <random>
#include <stdio.h>
#include <sys/time.h>
#include <sys/times.h>

// GooFit stuff
#include <goofit/Application.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/PDFs/physics/DalitzPlotPdf.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>

using namespace std;
using namespace GooFit;

TCanvas *foo;
TCanvas *foodal;
timeval startTime, stopTime, totalTime;
clock_t startCPU, stopCPU;
tms startProc, stopProc;
UnbinnedDataSet *data = 0;

Observable m12("m12", 0.0, 3.0);
Observable m13("m13", 0.0, 3.0);

EventNumber eventNumber("eventNumber");
bool fitMasses = false;
Variable fixedRhoMass("rho_mass", 0.7758, 0.01, 0.7, 0.8);
Variable fixedRhoWidth("rho_width", 0.1503, 0.01, 0.1, 0.2);

const fptype _mD0       = 1.86484;
const fptype piPlusMass = 0.13957018;
const fptype piZeroMass = 0.1349766;
const fptype D1Mass     = piZeroMass;
const fptype D2Mass     = piPlusMass;
const fptype D3Mass     = piPlusMass;
const fptype D1Mass2    = D1Mass * D1Mass;
const fptype D2Mass2    = D2Mass * D2Mass;
const fptype D3Mass2    = D3Mass * D3Mass;
const fptype MMass      = _mD0;
const fptype MMass2     = MMass * MMass;
const fptype MMass2inv  = 1. / MMass2;

// Constants used in more than one PDF component.
Variable motherM("motherM", MMass);
Variable dau1M("dau1M", D1Mass);
Variable dau2M("dau2M", D2Mass);
Variable dau3M("dau3M", D3Mass);
Variable massSum("massSum", MMass2 + D1Mass2 + D2Mass2 + D3Mass2); // = 3.53481
Variable constantOne("constantOne", 1);
Variable constantZero("constantZero", 0);

std::vector<PdfBase *> comps;

// I don't like Globals! Henry
int verbosity = 3;

GooPdf *kzero_veto = 0;
char strbuffer[1000];
double mesonRad = 1.5;
DalitzPlotPdf *signalDalitz;

void makeToyDalitzData(GooPdf *overallSignal, const int iSeed = 0, string datadir = ".", const int nTotal = 1e5);

DalitzPlotPdf *makeSignalPdf(GooPdf *eff = 0, bool fixAmps = false);

void makeToyDalitzData(GooPdf *overallSignal, const int iSeed, string datadir, const int nTotal) {
    std::vector<Observable> vars;
    vars.push_back(m12);
    vars.push_back(m13);
    vars.push_back(eventNumber);
    data = new UnbinnedDataSet(vars);
    UnbinnedDataSet currData(vars);
    std::vector<std::vector<double>> pdfValues;
    int ncount = 0;
    TRandom3 donram(iSeed);
    for(int i = 0; i < (m12.getNumBins()); ++i) {
        m12.setValue(m12.getLowerLimit() + (m12.getUpperLimit() - m12.getLowerLimit()) * (i + 0.5) / m12.getNumBins());
        for(int j = 0; j < m13.getNumBins(); ++j) {
            m13.setValue(m13.getLowerLimit()
                         + (m13.getUpperLimit() - m13.getLowerLimit()) * (j + 0.5) / m13.getNumBins());
            if(!inDalitz(m12.getValue(), m13.getValue(), _mD0, piZeroMass, piPlusMass, piPlusMass))
                continue;
            eventNumber.setValue(ncount);
            ncount++;
            currData.addEvent();
        }
    }
    signalDalitz->setDataSize(currData.getNumEvents());
    overallSignal->setData(&currData);

    pdfValues = overallSignal->getCompProbsAtDataPoints();
    TH2F dalitzpp0_dat_hist("dalitzpp0_dat_hist",
                            "",
                            m12.getNumBins(),
                            m12.getLowerLimit(),
                            m12.getUpperLimit(),
                            m13.getNumBins(),
                            m13.getLowerLimit(),
                            m13.getUpperLimit());
    dalitzpp0_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+}#pi^{0}) [GeV^{2}]");
    dalitzpp0_dat_hist.GetYaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV^{2}]");
    ncount = 0;
    ofstream writer;
    sprintf(strbuffer, "%s/dalitz_mytoyMC_%03d.txt", datadir.c_str(), iSeed);
    writer.open(strbuffer);
    vector<double> fIntegral;
    fIntegral.push_back(pdfValues[0][0]);
    Int_t ncells = pdfValues[0].size();
    for(unsigned int j = 1; j < ncells; ++j) {
        fIntegral.push_back(pdfValues[0][j] + fIntegral[j - 1]);
    }
    for(unsigned int j = 0; j < ncells; ++j)
        fIntegral[j] /= fIntegral[ncells - 1];
    ncount      = 0;
    int nEvents = donram.Poisson(nTotal);
    for(int iEvt = 0; iEvt < nEvents; iEvt++) {
        double r = donram.Rndm();
        // Binary search for fIntegral[cell-1] < r < fIntegral[cell]
        int lo = 0, hi = ncells - 1, mid = 0;
        while(lo <= hi) {
            mid = lo + (hi - lo) / 2;
            if(r <= fIntegral[mid] && (mid == 0 || r > fIntegral[mid - 1]))
                break;
            else if(r > fIntegral[mid])
                lo = mid + 1;
            else
                hi = mid - 1;
        }
        int j          = mid;
        double currm12 = currData.getValue(m12, j);
        currm12 += (m12.getUpperLimit() - m12.getLowerLimit()) * (donram.Rndm() - 0.5) / m12.getNumBins();
        double currm13 = currData.getValue(m13, j);
        currm13 += (m13.getUpperLimit() - m13.getLowerLimit()) * (donram.Rndm() - 0.5) / m13.getNumBins();
        eventNumber.setValue(ncount++);
        dalitzpp0_dat_hist.Fill(currm12, currm13);
        data->addEvent();
        writer << ncount - 1 << '\t' << currm12 << '\t' << currm13 << std::endl;
    }
    writer.close();
    std::cout << "Entries generated: " << data->getNumEvents() << std::endl;
    foodal->cd();
    foodal->SetLogz(false);
    dalitzpp0_dat_hist.Draw("colz");
    foodal->SaveAs("dalitzpp0_dat_temp.png");
}

void runToyGeneration(int numFile = 0) {
    m12 = Observable("m12", 0.0, 3.0);
    m12.setNumBins(1500);
    //  m12   = Variable("m12",   0.4, 3.0);
    m13 = Observable("m13", 0.0, 3.0);
    m13.setNumBins(1500);
    eventNumber  = EventNumber("eventNumber", 0, INT_MAX);
    signalDalitz = makeSignalPdf(0, true);
    vector<PdfBase *> comps;
    comps.clear();
    comps.push_back(signalDalitz);
    //  comps.push_back(sig0_jsugg);
    std::cout << "Creating overall PDF\n";
    ProdPdf *overallSignal = new ProdPdf("overallSignal", comps);
    gettimeofday(&startTime, NULL);
    startCPU = times(&startProc);
    //  makeToyDalitzData (signalDalitz);
    makeToyDalitzData(overallSignal, numFile);
    stopCPU = times(&stopProc);
    gettimeofday(&stopTime, NULL);
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
    std::vector<Observable> vars;
    vars.push_back(m12);
    vars.push_back(m13);
    vars.push_back(eventNumber);
    data = new UnbinnedDataSet(vars);

    std::ifstream reader;
    reader.open(toyFileName.c_str());
    assert(reader.good());
    std::string buffer;
    double dummy = 0;
    reader >> buffer;
    bool oldtype = false;
    if(buffer == "====")
        oldtype = true;
    else
        reader.seekg(0); // Set to the beginning
    while(!reader.eof()) {
        reader >> dummy;
        if(oldtype)
            reader >> dummy;
        reader >> m12;
        reader >> m13;
        if(oldtype) {
            for(int i = 0; i < 16; i++)
                reader >> dummy;
        }

        eventNumber.setValue(data->getNumEvents());
        data->addEvent();

        dalitzplot.Fill(m12.getValue(), m13.getValue());
    }
    reader.close();

    dalitzplot.SetStats(false);
    dalitzplot.Draw("colz");
    foodal->SaveAs("dalitzplot.png");
    foodal->SaveAs("dalitzplot.root");
}

GooPdf *makeKzeroVeto() {
    if(kzero_veto)
        return kzero_veto;

    Variable minimum("veto_min", 0.475 * 0.475);
    Variable maximum("veto_max", 0.505 * 0.505);
    VetoInfo kVetoInfo(minimum, maximum, PAIR_23);

    /*VetoInfo* kVetoInfo = new VetoInfo();
    kVetoInfo->cyclic_index = PAIR_23;
    kVetoInfo->minimum = Variable("veto_min", 0.475*0.475);
    kVetoInfo->maximum = Variable("veto_max", 0.505*0.505);*/

    vector<VetoInfo> vetos;
    vetos.push_back(kVetoInfo);
    kzero_veto = new DalitzVetoPdf("kzero_veto", m12, m13, motherM, dau1M, dau2M, dau3M, vetos);

    return kzero_veto;
}

DalitzPlotPdf *makeSignalPdf(GooPdf *eff, bool fixAmps) {
    DecayInfo3 dtop0pp;
    dtop0pp.motherMass   = MMass;
    dtop0pp.daug1Mass    = D1Mass;
    dtop0pp.daug2Mass    = D2Mass;
    dtop0pp.daug3Mass    = D3Mass;
    dtop0pp.meson_radius = 1.5;

    // Make a random number generater heres

    random_device rd;

    mt19937 mt(rd());

    normal_distribution<double> rand_gen(0.0, 0.1);

    auto rhop = new Resonances::RBW(
        "rhop", Variable("rhop_amp_real", 1), Variable("rhop_amp_imag", 0), fixedRhoMass, fixedRhoWidth, 1, PAIR_12);

    auto var_func = [&](std::string name, double start, double err) -> Variable {
        return fixAmps ? Variable(name, start) : Variable(name, start + rand_gen(mt), err, 0, 0);
    };

    ResonancePdf *rhom = new Resonances::RBW("rhom",
                                             var_func("rhom_amp_real", 0.714, 0.001),
                                             var_func("rhom_amp_imag", -0.025, 0.1),
                                             fixedRhoMass,
                                             fixedRhoWidth,
                                             1,
                                             PAIR_13);

    ResonancePdf *rho0 = new Resonances::RBW("rho0",
                                             var_func("rho0_amp_real", 0.565, 0.001),
                                             var_func("rho0_amp_imag", 0.164, 0.1),
                                             fixedRhoMass,
                                             fixedRhoWidth,
                                             1,
                                             PAIR_23);

    Variable sharedMass("rhop_1450_mass", 1.465, 0.01, 1.0, 2.0);
    Variable shareWidth("rhop_1450_width", 0.400, 0.01, 0.01, 5.0);

    ResonancePdf *rhop_1450 = new Resonances::RBW("rhop_1450",
                                                  var_func("rhop_1450_amp_real", -0.174, 0.001),
                                                  var_func("rhop_1450_amp_imag", -0.117, 0.1),
                                                  sharedMass,
                                                  shareWidth,
                                                  1,
                                                  PAIR_12);

    ResonancePdf *rho0_1450 = new Resonances::RBW("rho0_1450",
                                                  var_func("rho0_1450_amp_real", 0.325, 0.001),
                                                  var_func("rho0_1450_amp_imag", 0.057, 0.1),
                                                  sharedMass,
                                                  shareWidth,
                                                  1,
                                                  PAIR_23);

    ResonancePdf *rhom_1450 = new Resonances::RBW("rhom_1450",
                                                  var_func("rhom_1450_amp_real", 0.788, 0.001),
                                                  var_func("rhom_1450_amp_imag", 0.226, 0.1),
                                                  sharedMass,
                                                  shareWidth,
                                                  1,
                                                  PAIR_13);

    sharedMass = Variable("rhop_1700_mass", 1.720, 0.01, 1.6, 1.9);
    shareWidth = Variable("rhop_1700_width", 0.250, 0.01, 0.1, 1.0);

    ResonancePdf *rhop_1700 = new Resonances::RBW("rhop_1700",
                                                  var_func("rhop_1700_amp_real", 2.151, 0.001),
                                                  var_func("rhop_1700_amp_imag", -0.658, 0.1),
                                                  sharedMass,
                                                  shareWidth,
                                                  1,
                                                  PAIR_12);

    ResonancePdf *rho0_1700 = new Resonances::RBW("rho0_1700",
                                                  var_func("rho0_1700_amp_real", 2.400, 0.001),
                                                  var_func("rho0_1700_amp_imag", -0.734, 0.1),
                                                  sharedMass,
                                                  shareWidth,
                                                  1,
                                                  PAIR_23);

    ResonancePdf *rhom_1700 = new Resonances::RBW("rhom_1700",
                                                  var_func("rhom_1700_amp_real", 1.286, 0.001),
                                                  var_func("rhom_1700_amp_imag", -1.532, 0.1),
                                                  sharedMass,
                                                  shareWidth,
                                                  1,
                                                  PAIR_13);

    auto f0_980 = new Resonances::FLATTE("f0_980",
                                         var_func("f0_980_amp_real", 0.008 * (-MMass2), 0.001),
                                         var_func("f0_980_amp_imag", -0.013 * (-MMass2), 0.1),
                                         Variable("f0_980_mass", 0.9399 /*0.980*/, 0.01, 0.8, 1.2),
                                         Variable("f0_980_width", 0.199 /*0.044*/, 0.001, 0.001, 0.08),
                                         Variable("f0_980_rg2og1", 3.0, 0.1, 1e-3, 10),
                                         PAIR_23,
                                         false);

    ResonancePdf *f0_1370 = new Resonances::RBW("f0_1370",
                                                var_func("f0_1370_amp_real", -0.058 * (-MMass2), 0.001),
                                                var_func("f0_1370_amp_imag", 0.026 * (-MMass2), 0.1),
                                                Variable("f0_1370_mass", 1.434, 0.01, 1.2, 1.6),
                                                Variable("f0_1370_width", 0.173, 0.01, 0.01, 0.4),
                                                (unsigned int)0,
                                                PAIR_23);

    ResonancePdf *f0_1500 = new Resonances::RBW("f0_1500",
                                                var_func("f0_1500_amp_real", 0.057 * (-MMass2), 0.001),
                                                var_func("f0_1500_amp_imag", 0.012 * (-MMass2), 0.1),
                                                Variable("f0_1500_mass", 1.507, 0.01, 1.3, 1.7),
                                                Variable("f0_1500_width", 0.109, 0.01, 0.01, 0.3),
                                                (unsigned int)0,
                                                PAIR_23);

    ResonancePdf *f0_1710 = new Resonances::RBW("f0_1710",
                                                var_func("f0_1710_amp_real", 0.070 * (-MMass2), 0.001),
                                                var_func("f0_1710_amp_imag", 0.087 * (-MMass2), 0.1),
                                                Variable("f0_1710_mass", 1.714, 0.01, 1.5, 2.9),
                                                Variable("f0_1710_width", 0.140, 0.01, 0.01, 0.5),
                                                (unsigned int)0,
                                                PAIR_23);

    ResonancePdf *f2_1270 = new Resonances::RBW("f2_1270",
                                                var_func("f2_1270_amp_real", -1.027 * (-MMass2inv), 0.001),
                                                var_func("f2_1270_amp_imag", -0.162 * (-MMass2inv), 0.1),
                                                Variable("f2_1270_mass", 1.2754, 0.01, 1.0, 1.5),
                                                Variable("f2_1270_width", 0.1851, 0.01, 0.01, 0.4),
                                                2,
                                                PAIR_23);

    ResonancePdf *f0_600 = new Resonances::RBW("f0_600",
                                               var_func("f0_600_amp_real", 0.068 * (-MMass2), 0.001),
                                               var_func("f0_600_amp_imag", 0.010 * (-MMass2), 0.1),
                                               Variable("f0_600_mass", 0.500, 0.01, 0.3, 0.7),
                                               Variable("f0_600_width", 0.400, 0.01, 0.2, 0.6),
                                               (unsigned int)0,
                                               PAIR_23);

    ResonancePdf *nonr = new Resonances::NonRes(
        "nonr", var_func("nonr_amp_real", 0.5595 * (-1), 0.001), var_func("nonr_amp_imag", -0.108761 * (-1), 0.1));

    dtop0pp.resonances.push_back(nonr);
    dtop0pp.resonances.push_back(rho0);
    dtop0pp.resonances.push_back(rhom);
    dtop0pp.resonances.push_back(rhop_1450);
    dtop0pp.resonances.push_back(rho0_1450);
    dtop0pp.resonances.push_back(rhom_1450);
    dtop0pp.resonances.push_back(rhop_1700);
    dtop0pp.resonances.push_back(rho0_1700);
    dtop0pp.resonances.push_back(rhom_1700);
    dtop0pp.resonances.push_back(f0_980);
    dtop0pp.resonances.push_back(f0_1370);
    dtop0pp.resonances.push_back(f0_1500);
    dtop0pp.resonances.push_back(f0_1710);
    dtop0pp.resonances.push_back(f2_1270);
    dtop0pp.resonances.push_back(f0_600);

    if(!fitMasses) {
        for(vector<ResonancePdf *>::iterator res = dtop0pp.resonances.begin(); res != dtop0pp.resonances.end(); ++res) {
            (*res)->setParameterConstantness(true);
        }
    }

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
    comps.push_back(eff);
    if(!kzero_veto)
        makeKzeroVeto();
    comps.push_back(kzero_veto);
    ProdPdf *effWithVeto = new ProdPdf("effWithVeto", comps);

    return new DalitzPlotPdf("signalPDF", m12, m13, eventNumber, dtop0pp, effWithVeto);
}

double DalitzNorm(GooPdf *overallSignal, int N, double phi) {
    double max_pdf_value = phi * 1.1;

    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> xyvalues(0.0, 3.0);
    uniform_real_distribution<double> rpdfValues(0.0, max_pdf_value);

    std::vector<Observable> vars;
    vars.push_back(m12);
    vars.push_back(m13);
    vars.push_back(eventNumber);

    std::vector<fptype> rpdfValuesvec;

    UnbinnedDataSet data(vars);
    eventNumber = 0;

    std::ofstream wt2("toyData.txt");

    for(int i = 0; i < N; i++) {
        m12 = xyvalues(mt);
        m13 = xyvalues(mt);

        if(inDalitz(m12.getValue(), m13.getValue(), _mD0, piZeroMass, piPlusMass, piPlusMass) == 1) {
            data.addEvent();
            eventNumber.setValue(eventNumber.getValue() + 1);
        }
    }

    overallSignal->setData(&data);
    signalDalitz->setDataSize(data.getNumEvents());
    std::vector<std::vector<double>> pdfValues = overallSignal->getCompProbsAtDataPoints();

    int hit  = 0;
    int miss = 0;

    double buffer = 0;

    for(int k = 0; k < pdfValues[0].size(); k++) {
        buffer += pdfValues[0][k];

        if(rpdfValues(mt) < pdfValues[0][k]) {
            hit++;
            data.loadEvent(k);
            wt2 << m12.getValue() << " " << m13.getValue() << std::endl;

        } else {
            miss++;
        }
    }

    wt2.close();

    double H = (m12.getUpperLimit() - m12.getLowerLimit()) * (m13.getUpperLimit() - m13.getLowerLimit()); // Area

    //		hit = data.getNumEvents();

    double integral = H * hit / N;

    return integral;
}

void runIntegration(int n = 100) {
    // TApplication* rootapp = new TApplication("rootapp",&argc,argv);

    signalDalitz = makeSignalPdf(0, true);

    std::vector<PdfBase *> comps;
    comps.clear();
    comps.push_back(signalDalitz);

    ProdPdf *overallSignal = new ProdPdf("overallSignal", comps);

    int N = 100000;

    ofstream wt3("integral.txt");

    double integral = 0;

    std::cout << "start ! " << std::endl;

    for(int i = 0; i < n; i++) {
        integral = DalitzNorm(overallSignal, N, 1.8);

        wt3 << integral << std::endl;

        std::cout << "integral " << i << " = " << integral << std::endl;
    }

    wt3.close();

    std::cout << "end !\n\n" << std::endl;

    TCanvas *c = new TCanvas("c", "", 800, 800);
    TH1F *hist = new TH1F();
    TTree tree("integral.txt", "x");
    tree.ReadFile("integral.txt", "Integral");
    hist->SetTitle("Integral");
    hist->GetXaxis()->SetTitle("Values");
    hist->GetYaxis()->SetTitle("Frequency");
    tree.Draw("Integral>>hist");
    c->SaveAs("histogram.png");

    TCanvas *d  = new TCanvas("d", "", 800, 800);
    TH1F *hist2 = new TH1F();
    TTree tree2("toyData.txt", "m12:m13");
    tree2.ReadFile("toyData.txt", "m12:m13");
    tree2.Draw("m12:m13>>hist2", "", "colz");
    d->SaveAs("toyDalitz.png");

    // std::cout << "mean: " << hist->GetMean() << "\n";
    // std::cout << "rms: " << hist->GetRMS() << "\n";

    // rootapp->Run();
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
    ht->Scale(hd->Integral() / ht->Integral());
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

void makeToyDalitzPdfPlots(GooPdf *overallSignal, string plotdir = "plots") {
    TH1F m12_dat_hist("m12_dat_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    m12_dat_hist.SetStats(false);
    m12_dat_hist.SetMarkerStyle(8);
    m12_dat_hist.SetMarkerSize(1.2);
    m12_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
    sprintf(strbuffer, "Events / %.1f MeV", 1e3 * m12_dat_hist.GetBinWidth(1));
    m12_dat_hist.GetYaxis()->SetTitle(strbuffer);
    TH1F m12_pdf_hist("m12_pdf_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
    m12_pdf_hist.SetStats(false);
    m12_pdf_hist.SetLineColor(kBlue);
    m12_pdf_hist.SetLineWidth(3);
    TH1F m13_dat_hist("m13_dat_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m13_dat_hist.SetStats(false);
    m13_dat_hist.SetMarkerStyle(8);
    m13_dat_hist.SetMarkerSize(1.2);
    m13_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV]");
    sprintf(strbuffer, "Events / %.1f MeV", 1e3 * m13_dat_hist.GetBinWidth(1));
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
    sprintf(strbuffer, "Events / %.1f MeV", 1e3 * m13_dat_hist.GetBinWidth(1));
    m23_dat_hist.GetYaxis()->SetTitle(strbuffer);
    TH1F m23_pdf_hist("m23_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
    m23_pdf_hist.SetStats(false);
    m23_pdf_hist.SetLineColor(kBlue);
    m23_pdf_hist.SetLineWidth(3);
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
    dalitzpp0_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
    dalitzpp0_dat_hist.GetYaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV]");
    TH2F dalitzpp0_pdf_hist("dalitzpp0_pdf_hist",
                            "",
                            m12.getNumBins(),
                            m12.getLowerLimit(),
                            m12.getUpperLimit(),
                            m13.getNumBins(),
                            m13.getLowerLimit(),
                            m13.getUpperLimit());
    /*  dalitzpp0_pdf_hist.GetXaxis()->SetTitle("m^{2}(K^{-} #pi^{0}) [GeV^{2}]");
        dalitzpp0_pdf_hist.GetYaxis()->SetTitle("m^{2}(K^{-} #pi^{+}) [GeV^{2}]");*/
    dalitzpp0_pdf_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV^{2}]");
    dalitzpp0_pdf_hist.GetYaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV^{2}]");
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
            if(!inDalitz(m12.getValue(), m13.getValue(), _mD0, piZeroMass, piPlusMass, piPlusMass))
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
    foodal->cd();
    foodal->SetLogz(false);
    dalitzpp0_pdf_hist.Draw("colz");
    std::string command = "mkdir -p " + plotdir;
    if(system(command.c_str()) != 0)
        throw GooFit::GeneralError("Making plot directory {} failed", plotdir);
    foodal->SaveAs((plotdir + "/dalitzpp0_pdf.png").c_str());
    /*  m12_pdf_hist.Draw("");
        foodal->SaveAs((plotdir + "/m12_pdf_hist.png").c_str());
        m13_pdf_hist.Draw("");
        foodal->SaveAs((plotdir + "/m13_pdf_hist.png").c_str());
        if (!data) return;*/
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
    foodal->SaveAs((plotdir + "/dalitzpp0_dat.png").c_str());

    drawFitPlotsWithPulls(&m12_dat_hist, &m12_pdf_hist, plotdir);
    drawFitPlotsWithPulls(&m13_dat_hist, &m13_pdf_hist, plotdir);
    drawFitPlotsWithPulls(&m23_dat_hist, &m23_pdf_hist, plotdir);
}

void runToyFit(std::string toyFileName) {
    m12 = Observable("m12", 0, 3);
    m13 = Observable("m13", 0, 3);
    m12.setNumBins(300);
    m13.setNumBins(300);
    eventNumber = EventNumber("eventNumber", 0, INT_MAX);
    getToyData(toyFileName);

    // EXERCISE 1 (real part): Create a PolynomialPdf which models
    // the efficiency you imposed in the preliminary, and use it in constructing
    // the signal PDF.

    // EXERCISE 2: Create a K0 veto function and use it as the efficiency.

    // EXERCISE 3: Make the efficiency a product of the two functions
    // from the previous exercises.

    signalDalitz = makeSignalPdf();
    comps.clear();
    comps.push_back(signalDalitz);
    ProdPdf *overallSignal = new ProdPdf("overallSignal", comps);
    overallSignal->setData(data);
    signalDalitz->setDataSize(data->getNumEvents());

    FitManager datapdf(overallSignal);

    gettimeofday(&startTime, NULL);
    startCPU = times(&startProc);
    datapdf.setVerbosity(verbosity);

    // Maybe make optional? With a command line switch?
    datapdf.fit();
    stopCPU = times(&stopProc);
    gettimeofday(&stopTime, NULL);

    // Get the fractions w/ uncertainties
    // vector<double> fracList;
    // signalDalitz->getFractions(fracList);
    /*  const int nRes = fracList.size();
        vector <float> fractions[nRes];
        float mean[nRes];
        float rms[nRes];
        for (int ii=0;ii<nRes;ii++) mean[ii] = rms[ii] = 0;
        for (int ii=0;ii<nSamples;ii++){
        datapdf.loadSample(ii);
        signalDalitz->getFractions(fracList);
        for (int jj=0;jj<nRes; jj++) {
        fractions[jj].push_back(fracList[jj]);
        mean[jj] += fracList[jj];
        rms[jj] += fracList[jj]*fracList[jj];
        }
        }
        TH1F* hFracs[nRes];
        TFile * froot = new TFile("fractionHists.root", "recreate");
        for (int ii=0;ii<nRes;ii++) {
        mean[ii] /= nSamples;
        rms[ii] = sqrt(rms[ii]/nSamples-mean[ii]*mean[ii]);
        sprintf(strbuffer, "hfrac_res%d", ii);
        hFracs[ii] = new TH1F(strbuffer, "", 100, mean[ii]-4*rms[ii], mean[ii]+4*rms[ii]);
        for (int jj=0;jj<nSamples;jj++)
        hFracs[ii]->Fill(fractions[ii][jj]);
        hFracs[ii]->Write();
        }
        froot->Close();*/

    makeToyDalitzPdfPlots(overallSignal);
}

int main(int argc, char **argv) {
    GooFit::Application app{"D2K3_toy", argc, argv};
    app.add_option("-v,--verbose", verbosity, "Set the verbosity (to 0 for example", true);

    int fit_value;
    std::string name = "dalitz_mytoyMC_000.txt";

    auto fit = app.add_subcommand("fit");
    fit->add_option("-i,--int", fit_value, "A number to load");
    auto name_opt = fit->add_option("-n,--name,name", name, "The filename to load", true)->excludes("--int");

    int number;
    auto run = app.add_subcommand("run");
    run->add_option("N", number, "How many integrations do you want? N=100")->required();

    int value;
    auto gen = app.add_subcommand("gen");
    gen->add_option("value", value, "The number to generate")->required();

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
    gStyle->SetPalette(1, 0);
    gStyle->SetOptStat("RM");
    foo    = new TCanvas();
    foodal = new TCanvas();
    foodal->Size(10, 10);

    if(*fit)
        runToyFit(name);
    if(*gen)
        runToyGeneration(value);
    if(*run) {
        CLI::AutoTimer timer("Integration");
        runIntegration(number);
    }

    // Print total minimization time
    double myCPU    = stopCPU - startCPU;
    double totalCPU = myCPU;

    timersub(&stopTime, &startTime, &totalTime);
    std::cout << "Wallclock time  : " << totalTime.tv_sec + totalTime.tv_usec / 1000000.0 << " seconds." << std::endl;
    std::cout << "CPU time: " << (myCPU / CLOCKS_PER_SEC) << std::endl;
    std::cout << "Total CPU time: " << (totalCPU / CLOCKS_PER_SEC) << std::endl;
    myCPU = stopProc.tms_utime - startProc.tms_utime;
    std::cout << "Processor time: " << (myCPU / CLOCKS_PER_SEC) << std::endl;

    return 0;
}
