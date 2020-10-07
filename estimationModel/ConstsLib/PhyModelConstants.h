//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see http://www.gnu.org/licenses/.
//

#ifndef PhyModelConstants_H_
#define PhyModelConstants_H_

#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include<string>
#include <iostream>
#include <cmath>
#include "inet/physicallayer/ieee80211/packetlevel/errormodel/Ieee80211NistErrorModel.h"

using namespace inet::physicallayer;
namespace PhyModelConstants{};
namespace PhyModelConstants {

    //###METHODS
    double computeTwoRayInterference(double distance , double z);
    double getBER(double datarate,double snr);
    double estimateSNR(double datarate, double pt, double distance, std::string propModel,double z);
    double computeReceivePower(double pt,double d ,std::string propModel, double z);
    //###CONST
    // light of speed
    double const lightSpeed=2.99792*pow(10,8);
    // carrier frequency
    double const f=5.9*pow(10,9);
    //length wave
    double const lambda=(lightSpeed/f);
    // Antenna receive Gain
    double const GR=1;
    // Anntena transmit Gain
    double const GT=1;
    // const of boltzman
    double const Kb=1.38064852*pow(10,-23);
    // Kelvin
    double const T=335;
    // spectral density
    double const  N0=Kb*T;
    double const Mbps=pow(10,6);



}
#endif
