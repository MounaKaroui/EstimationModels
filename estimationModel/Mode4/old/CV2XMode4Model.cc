#include "CV2XMode4Model.h"
#include<vector>
#include <algorithm>    // std::find_if

namespace CV2XMode4Model {

//void CV2XMode4_common(double lambda, double Pt, std::vector<double> distance, double Psen, double step_dB, double noise, double coding, double& deltaHD, double& deltaSEN, double& deltaPRO)
//{
//  double D, PDF_SNR, PL_E_R, SNR, std_dev_E_R ;
//  D = distance.size() ;
//
//  deltaHD= lambda/1000.0 ;
//
//  get_PL_SH(distance, PL_E_R, std_dev_E_R) ;
//  deltaSEN = 0.5*(1-erf((Pt-PL_E_R-Psen)/(std_dev_E_R*sqrt(2)))) ;
//  //get_SINRdistribution(Pt-PL_E_R, -180, std_dev_E_R, 3, noise, Psen, step_dB, SNR, PDF_SNR) ;
//
//  //deltaPRO = get_BLER(SNR, PDF_SNR, coding, step_dB) ;
//}


//void get_PL_SH(std::vector<double> distance, std::vector<double>& PL, std::vector<double>& std_dev)
//{
//  double  PLfree,dBP, i ;
//  for(int j=0; j<distance.size();j++)
//      distance.at(j) = abs(distance.at(j)) ;
//
//  dBP = 4*(hBS-environmentHeight)*(hMS-environmentHeight)*fc/c ;
//  auto it = std::find_if(std::begin(distance), std::end(distance), [](int k){return k <3 ;});
////  while (it != std::end(distance)) {
////     results.emplace_back(std::distance(std::begin(distance), it));
////     it = std::find_if(std::next(it), std::end(v), [](int i){return i > 5;});
////  }
//  i = *it + 1 ;
//  distance.at(i) = 3 ;
//  auto it = std::find_if(std::begin(distance), std::end(distance), [](int k){return k <dBP ;});
//  i = *it + 1 ;
//
//  PL.at(i) = 22.7*arma::log10(distance.at(i))+27+20*arma::log10(fc/1e9) ;
//  std_dev.at(i) = 3 ;
//  auto it = std::find_if(std::begin(distance), std::end(distance), [](int k){return k >=dBP ;});
//  i = *it + 1 ;
//  PL.at(i) = 40*arma::log10(distance.at(i))+7.56-17.3*arma::log10(hBS-environmentHeight)-17.3*arma::log10(hMS-environmentHeight)+2.7*arma::log10(fc/1e9) ;
//  std_dev.at(i) = 3 ;
//
//  PLfree = 20*arma::log10(distance.at(i))+46.4+20*arma::log10(fc*1e-9/5.0) ;
//  auto it = std::find_if(std::begin(PLfree), std::end(PLfree), [](double k){return k>PL ;});
//  i = *it + 1 ;
//  PL.at(i) = PLfree.at(i) ;
//}



}
