{
#include <iostream>
#include <fstream>

std::ifstream reader("fitResults.txt");

double amp_real_f0,amp_real_f0_e,amp_img_f0,amp_img_f0_e,amp_real_f0X,amp_real_f0X_e,amp_img_f0X,amp_img_f0X_e,fcn,norm,bins;

size_t nEvents = 100;

std::vector<double> vec_amp_real_f0;
std::vector<double> vec_amp_real_f0_e;
std::vector<double> vec_amp_img_f0;
std::vector<double> vec_amp_img_f0_e;
std::vector<double> vec_amp_real_f0X;
std::vector<double> vec_amp_real_f0X_e;
std::vector<double> vec_amp_img_f0X;
std::vector<double> vec_amp_img_f0X_e;
std::vector<double> vec_fcn;
std::vector<double> vec_norm;
std::vector<double> vec_bins;



while(!reader.eof()){

reader >> amp_real_f0 >> amp_real_f0_e >> amp_img_f0 >> amp_img_f0_e >> amp_real_f0X >> amp_real_f0X_e >> amp_img_f0X >> amp_img_f0X_e >> fcn >> norm >> bins; 

vec_amp_real_f0.push_back(amp_real_f0);
vec_amp_real_f0_e.push_back(amp_real_f0_e);
vec_amp_img_f0.push_back(amp_img_f0);
vec_amp_img_f0_e.push_back(amp_img_f0_e);

vec_amp_real_f0X.push_back(amp_real_f0X);
vec_amp_real_f0X_e.push_back(amp_real_f0X_e);
vec_amp_img_f0X.push_back(amp_img_f0X);
vec_amp_img_f0X_e.push_back(amp_img_f0X_e);

vec_fcn.push_back(fcn);
vec_norm.push_back(norm);
vec_bins.push_back(bins);
}

auto amp_real_f0_min = std::min_element(vec_amp_real_f0.begin(),vec_amp_real_f0.end());
auto amp_real_f0_max = std::max_element(vec_amp_real_f0.begin(),vec_amp_real_f0.end());
auto amp_img_f0_min = std::min_element(vec_amp_img_f0.begin(),vec_amp_img_f0.end());
auto amp_img_f0_max = std::max_element(vec_amp_img_f0.begin(),vec_amp_img_f0.end());

auto amp_real_f0X_min = std::min_element(vec_amp_real_f0X.begin(),vec_amp_real_f0X.end());
auto amp_real_f0X_max = std::max_element(vec_amp_real_f0X.begin(),vec_amp_real_f0X.end());
auto amp_img_f0X_min = std::min_element(vec_amp_img_f0X.begin(),vec_amp_img_f0X.end());
auto amp_img_f0X_max = std::max_element(vec_amp_img_f0X.begin(),vec_amp_img_f0X.end());

auto fcn_min = std::min_element(vec_fcn.begin(),vec_fcn.end());
auto fcn_max = std::max_element(vec_fcn.begin(),vec_fcn.end());

auto norm_min = std::min_element(vec_norm.begin(),vec_norm.end());
auto norm_max = std::max_element(vec_norm.begin(),vec_norm.end());

auto bins_min = std::min_element(vec_bins.begin(),vec_bins.end());
auto bins_max = std::max_element(vec_bins.begin(),vec_bins.end());

TFile f("normStudies.root","RECREATE");
TTree s("Tree","Tree");
s.ReadFile("fitResults.txt","amp_real_f0:amp_real_f0_e:amp_img_f0:amp_img_f0_e:amp_real_f0X:amp_real_f0X_e:amp_img_f0X:amp_img_f0X_e:fcn:norm:bins");
s.Write();

//some histos
TH2F amp_r_vs_bins_f0("amp_r_vs_bins_f0","amp_real_vs_bins_f0",nEvents,*bins_min,*bins_max,nEvents,*amp_real_f0_min,*amp_real_f0_max);
TH2F amp_i_vs_bins_f0("amp_i_vs_bins_f0","amp_img_vs_bins_f0",nEvents,*bins_min,*bins_max,nEvents,*amp_img_f0_min,*amp_img_f0_max);
TH2F amp_r_vs_bins_f0X("amp_r_vs_bins_f0X","amp_real_vs_bins_f0X",nEvents,*bins_min,*bins_max,nEvents,*amp_real_f0X_min,*amp_real_f0X_max);
TH2F amp_i_vs_bins_f0X("amp_i_vs_bins_f0X","amp_img_vs_bins_f0X",nEvents,*bins_min,*bins_max,nEvents,*amp_img_f0X_min,*amp_img_f0X_max);
TH2F fcn_vs_bins("fcn_vs_bins","fcn_vs_bins",nEvents,*bins_min,*bins_max,nEvents,*fcn_min,*fcn_max);
TH2F norm_vs_bins("norm_vs_bins","norm_vs_bins",nEvents,*bins_min,*bins_max,nEvents,*norm_min,*norm_max);


amp_r_vs_bins_f0.SetMarkerStyle(3);
s.Draw("amp_real_f0:bins>>amp_r_vs_bins_f0");
amp_r_vs_bins_f0.Write();



amp_i_vs_bins_f0.SetMarkerStyle(3);
s.Draw("amp_img_f0:bins>>amp_i_vs_bins_f0");
amp_i_vs_bins_f0.Write();


amp_r_vs_bins_f0X.SetMarkerStyle(3);
s.Draw("amp_real_f0X:bins>>amp_r_vs_bins_f0X");
amp_r_vs_bins_f0X.Write();


amp_i_vs_bins_f0X.SetMarkerStyle(3);
s.Draw("amp_img_f0X:bins>>amp_i_vs_bins_f0X");
amp_i_vs_bins_f0X.Write();


fcn_vs_bins.SetMarkerStyle(3);
s.Draw("FCN:bins>>fcn_vs_bins");
fcn_vs_bins.Write();


norm_vs_bins.SetMarkerStyle(3);
s.Draw("norm:bins>>norm_vs_bins");
norm_vs_bins.Write();

f.Close();

}
