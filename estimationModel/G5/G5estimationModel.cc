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

#include "G5estimationModel.h"



namespace G5estimationModel
{


double estimateProb(double n, std::string  type)
{

    double tau=estimateTau();

    if(type=="ps")
    {
        return n*tau*pow(1-tau,n-1);
    }

    if(type=="pi")
    {
        return pow(1-tau,n);
    }

    if(type=="pc")
    {
        return 1-tau*pow(1-tau,n-1);
    }

}

double estimateTau()
{
	return 1/(CW+1);
}

double estimateTransmissionTime(double payload, double dataRate, std::string type)
{
    if(type=="T")
    {
	return (payload/dataRate + AIFS);
    }

    // Ts
    if(type=="T_s")
    {
    return (payload/dataRate + DIFS);
    }
    // Tc
    if(type=="T_c")
    {
        return (payload/dataRate + EIFS);
    }
}


double estimatePe(double dataRate,double payload,double distance,double pt, std::string propModel,double z)
{
    /// BER header

    // datarate fix header=3 Mbps
    double header=0 ; // TODO search for the header value

    double snr_h=PhyModelConstants::estimateSNR(3*Mbps, pt, distance, propModel, z);
    double BER_h=PhyModelConstants::getBER(3*Mbps, snr_h); // BER header

    /// BER payload
    double snr_p=PhyModelConstants::estimateSNR(dataRate, pt, distance, propModel, z);
    double BER_p=PhyModelConstants::getBER(dataRate, snr_p); // BER payload

    /// channel error probability pe
    double pe=pow(1-BER_h,header)*pow(1-BER_p,payload);

    return pe;
}


double estimateHiddenTerminals(double nh, double lamda, double payload, double dataRate)
{
    double pi=estimateProb(nh, "pi");
    double Tdata=payload/dataRate;
    return (pi*exp(-1*lamda*(Tdata-DIFS)));
}





double estimateDelayWithoutChannelError(double n, double lamda,double payload, double dataRate, double pb)
{
    double tau=estimateTau();
    double T=estimateTransmissionTime(payload,dataRate,"T");
    double B=(CW-1)/2*(sigma+(1-pow(1-tau,n-1))*T); // mean Backoff durations
    double T_res=T/(1-exp(-1*lamda*T))-1/lamda; // residual life time of ongoing transmission
    double A=(B+T_res)*pb + (1-pb)*DIFS; // mean access delay
    return A;
}




void estimateAccessDelay(double n, double lamda,double payload, double dataRate,
		double distance,double pt, std::string propModel,double z, double& meanA, double& varA,double& totDelay)
{

	double ps=estimateProb(n,"ps");
	double pi=estimateProb(n,"pi");
	double Ts=estimateTransmissionTime(payload, dataRate,"T_s");
	double Tc=estimateTransmissionTime(payload, dataRate,"T_c");
	double pc=estimateProb(n,"pc");
	double ph=estimateHiddenTerminals(n-1,lamda,payload,dataRate); // nh=n-1
	double pe=estimatePe(dataRate, payload, distance, pt, propModel, z);
	double pl=1-(1-pc)*(1-pe)*ph; // system losses
	double pb=(Ts*ps+Tc*pl)/(sigma*pi+Ts*ps+Tc*pl); // busy probability

	double tau=estimateTau(); // tau value
	double T=estimateTransmissionTime(payload,dataRate,"T"); // Transmission time
	double B=(CW-1)/2*(sigma+(1-pow(1-tau,n-1))*T); // mean Backoff durations
	double T_res=T/(1-exp(-1*lamda*T))-1/lamda; // residual life time of ongoing transmission
	meanA=(B+T_res)*pb + (1-pb)*DIFS; // mean access delay

	double Y=(1-pow((1-tau),n-1))*T;
	double varY=(1-pow((1-tau),n-1))*(pow((1-tau),n-1))*T*T;

	double varB=varY*(CW-1)/2+(sigma+Y)*(sigma+Y)*(CW*CW-1)/12;
	double varTr=1/(lamda*lamda) -T*T*exp(-1*lamda*T)/((1-exp(-1*lamda*T))*(1-exp(-1*lamda*T)));

	varA= varB + (meanA-B)*(meanA-B);
	double s=meanA+T;
	double meanQ=lamda*(varA+s*s)/2*(1-lamda*s);
	totDelay=meanA+meanQ;
}





void estimateThroughput(double Lp, double Lh, double Rh, double R,
		double distance, double pt, std::string propModel, double z, double n, double nh,double lamda, double& th)
{


	double pc=estimateProb(n,"pc");
	double pe=estimatePe(R, Lp, distance, pt,propModel, z);
	double ph=estimateHiddenTerminals(nh, lamda, Lp, R);

	double pd=(1-pc)*(1-pe)*ph;
	double meanA, varA, totDelay;
	estimateAccessDelay(n, lamda, Lp, R, distance, pt, propModel, z,meanA,varA,totDelay);
	th=(pd*Lp)/(Lp/R + Lh/Rh + meanA);



}

};






