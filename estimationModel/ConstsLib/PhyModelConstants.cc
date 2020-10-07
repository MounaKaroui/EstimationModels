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

#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include<string>
#include <iostream>
#include <cmath>
#include "PhyModelConstants.h"

namespace PhyModelConstants {
constexpr double squared(double x) { return x * x; }


double computeTwoRayInterference(double distance , double z)
{
    const double h_sum = 2*z;
    const double h_dif =0;
    const double epsilon_r=1.02;
    const double d_los = sqrt(distance * distance + h_dif * h_dif);
    const double  d_ref = sqrt(distance * distance + h_sum * h_sum);
    const double phi = 2 * M_PI * (d_los - d_ref) /lambda;

    const double sin_theta = h_sum / d_ref;
    const double cos_theta = distance / d_ref;
    const double gamma = sqrt(epsilon_r - cos_theta * cos_theta);
    const double Gamma = (sin_theta - gamma) / (sin_theta + gamma);

    const double space = squared(4 * M_PI * distance /  lambda);
    const double interference = squared(1 + Gamma * cos(phi)) + Gamma * Gamma * squared(sin(phi));
    return interference / space;
}

double getBER(double datarate, double snr)
{
    Ieee80211NistErrorModel* errorModel=new Ieee80211NistErrorModel();
    double BER=0;

    if(datarate==3*Mbps)
    {
        // BPSK
        BER=errorModel->getBpskBer(snr);
    }

    if((datarate==6*Mbps)||(datarate==9*Mbps))
    {
        // QPSK
        BER=errorModel->getQpskBer(snr);
    }

    if((datarate==12*Mbps)||(datarate==18*Mbps))
    {
        // 16QAM
        BER=errorModel->get16QamBer(snr);
    }

    if(datarate==24*Mbps)
    {
        // 64 QAM
        BER=errorModel->get64QamBer(snr) ;
    }

    return BER;
}

double computeReceivePower(double pt,double d ,std::string propModel, double z)
{
    double ht=z;
    double hr=z;
    double pr=0;
    if(propModel=="freespace")
    {
        // free space
        pr=pt*GT*GR*squared(lambda/(4*M_PI*d));
    }
    if(propModel=="tworays")
    {
        double pathloss=computeTwoRayInterference(d,z);
        double l=1;
        pr=std::min(l,pt*GT*GR*pathloss);
    }
    return pr;
}

double estimateSNR(double datarate, double pt, double distance, std::string propModel,double z)
{

        // received power
        double pr=computeReceivePower(pt, distance, propModel,z);
        double m=0;
        double snr=0;
        if(datarate==3*Mbps)
        {
            // modulation BPSK
            m=2;
            snr=(2*log(m))*pr/(datarate*N0);
        }

        if((datarate==6*Mbps)||(datarate==9*Mbps))
        {
            // modulation QPSK
            m=4;
            snr=(2*log(m))*pr/(datarate*N0);
        }

        if((datarate==12*Mbps)||(datarate==18*Mbps))
        {
            // modulation 16-QAM
            m=16;
            snr=(2*log(m))*pr/(datarate*N0);
        }

        if(datarate==24*Mbps)
        {
            // modulation 64-QAM
            m=64;
            snr=(2*log(m))*pr/(datarate*N0);
        }

        return snr;
}

};
