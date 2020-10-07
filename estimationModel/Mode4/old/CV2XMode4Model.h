#ifndef CV2XMODE4_H_
#define CV2XMODE4_H_
#include <armadillo>
#include <stdio.h>
#include <cstdlib>
#include<string>
#include <iostream>
#include <vector>
#include <algorithm>    // std::find_if

using namespace arma ;
namespace CV2XMode4Model{};
namespace CV2XMode4Model{
//
//    void CV2XMode4_common(double lambda, double Pt, std::vector<double> distance, double Psen,
//            double step_dB, double noise, double coding, double& deltaHD, double& deltaSEN, double& deltaPRO);

    // void get_PL_SH(std::vector<double> distance, std::vector<double>& PL, std::vector<double>& std_dev);

//    double get_BLER(double SINR, double PDF, double coding, double step_dB);
//
//
//    TYPE CV2XMode4_Step3(TYPE beta, TYPE lambda, TYPE Pt, TYPE S, TYPE distance,
//            TYPE Psen, TYPE step_dB, TYPE noise, TYPE coding, TYPE deltaPRO);
//
//    void CV2XMode4_Step2(TYPE beta, TYPE lambda, TYPE Pt, TYPE S, TYPE distance, TYPE Psen,
//            TYPE step_dB, TYPE noise, TYPE coding, TYPE deltaPRO, TYPE& deltaCOL, TYPE& CBR);
//
//    void get_SINRdistribution(TYPE Pr_dBm_avg, TYPE Pi_dBm_avg, TYPE std_dev_Pr, TYPE std_dev_Pi,
//            TYPE noise_dBm, TYPE sensingThreshold, TYPE step_dB, TYPE& SINR, TYPE& PDF_SINR);

       double const c = 3e8;
       double const fc = 5.91e9 ; // % Carrier frequency (Hz)
       double const  hBS = 1.5 ; //   % Transmitter antenna height (m)
       double const  hMS = 1.5 ;  //   % Receiver antenna height (m)
       double const  environmentHeight = 0 ;//   % Average environmental height (m)

       double const Psen = -90.5;             //  % Sensing threshold (dBm)
       double const step_dB = 0.1; //% Discrete steps to compute the PDF of the SNR and SINR (dB)


}
#endif
