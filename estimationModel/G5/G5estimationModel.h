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

#ifndef G5ESTIMATIONMODEL_H_
#define G5ESTIMATIONMODEL_H_

#include <omnetpp.h>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include<string>
#include <iostream>
#include <cmath>
#include "artery/estimationModel/ConstsLib/PhyModelConstants.h"

using namespace omnetpp;
using namespace PhyModelConstants;


namespace G5estimationModel{};
namespace G5estimationModel
{

    double estimateTransmissionTime(double payload, double datarate, std::string type);
    double estimateTau();
    /// Busy probability
    double estimateProb(double n, std::string type);
    /// Losses
    double estimateHiddenTerminals(double nh, double lamda, double payload, double dataRate);
    /// for error model computation
    double estimatePe(double dataRate,double payload,double distance,double pt, std::string propModel,double z);

    void estimateAccessDelay(double n, double lamda,double payload, double dataRate,
    		double distance,double pt, std::string propModel,double z, double& meanA, double& varA,double& totDelay);

    double estimateDelayWithoutChannelError(double n, double lamda,double payload, double dataRate, double pb);


    double estimatePacketDelivery(double n, double dataRate, double payload,
    		double distance, double pt, double propModel, double z, double lamda, double nh);


    void estimateThroughput(double Lp, double Lh, double Rh, double R,double distance,
    		double pt, std::string propModel, double z, double n, double nh,double lamda,double& th);

    double const CW=1023;
    double const AIFS=58*pow(10,-6);
    double const Th=40*pow(10,-6);
    double const EIFS=188*pow(10,-6);
    double const DIFS=64*pow(10,-6); /// From paper: Performance Analysis of the IEEE 802.11p Multichannel MAC Protocol in vehicular Ad Hoc network
    double const sigma=13*pow(10,-6);


}; //namespace

#endif
