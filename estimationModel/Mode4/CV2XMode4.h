#ifndef CV2XMODE4_H_
#define CV2XMODE4_H_
#include <stdio.h>
#include <cstdlib>
#include<string>
#include <vector>
#define ARMA_64BIT_WORD
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#include <iostream>
using namespace std;
using namespace arma;


namespace CV2XMode4{};
namespace CV2XMode4{


    // Math tools xcorr
    void vec_push(arma::vec & v, double value) ;
    vec xCorr(vec x,vec y);
    vec xcorr(vec x, vec y,bool autoflag);
    int nextPower_2 ( unsigned  int  x);


    inline int pow2i(int x) { return ((x < 0) ? 0 : (1 << x)); }


    void get_SINRdistribution(vec Pr_dBm_avg, double Pi_dBm_avg, vec std_dev_Pr, double std_dev_Pi,
     double noise_dBm, double sensingThreshold, double step_dB, mat& SINR, mat& PDF_SINR);

    void get_PL_SH(vec d, vec & PL, vec &  std_dev);
    vec get_BLER(mat SINR, mat PDF, double coding, double step_dB);

    void CV2XMode4_common(double lambda, double Pt, vec distance, double Psen, double step_dB, double noise,
    double coding, vec & deltaHD, vec & deltaSEN, vec & deltaPRO);

    void CV2XMode4_Step2(double beta, double lambda, double Pt, double S, vec distance, double Psen, double step_dB,
    double noise, double coding, vec deltaPRO, vec& deltaCOL, double& CBR);

    vec CV2XMode4_Step3(double beta, double lambda, double Pt, double S, vec distance, double Psen,
    double step_dB, double noise, double coding, vec deltaPRO);

    void CV2XMode4Principal(double beta, double lamda, double Pt, double S, double B,
 vec distance, vec& PDR, vec & deltaHD, vec & deltaSEN, vec & deltaPRO,vec& deltaCOL, double& CBR);

    void recordOutput(vec var);
    double estimateDelay(double PDR, double lamda, double Ttr);
    vec calculateEWA(vec dist,double beta, double newSample);
   // double const Psen = -90.5;             //  % Sensing threshold (dBm)
   // double const step_dB = 0.1; //% Discrete steps to compute the PDF of the SNR and SINR (dB)



}
#endif
