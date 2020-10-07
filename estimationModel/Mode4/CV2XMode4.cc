#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include<string>
#include <iostream>
#include <cmath>
#include <vector>
#include "CV2XMode4.h"
#include<iomanip>
#include<limits>
#include <fstream>
#include<fenv.h>

namespace CV2XMode4
{

void vec_push(vec & v, double value) {
    vec av(1);
    av.at(0) = value;
    v.insert_rows(v.n_rows, av.row(0));
}

 int nextPower_2 ( unsigned  int  x)
{

double nextnum = ceil(log2(x)) ;

double result = pow  (2.0  ,  nextnum) ;

return  int ( result ) ;
}

vec xCorr(vec x,vec y)
{

    return arma::conv(x,arma::reverse(y));

}


vec xcorr(vec x, vec y,bool autoflag)
{
  int maxlag=0;
  int N = std::max(x.size(), y.size());
  //Compute the FFT size as the "next power of 2" of the input vector's length (max)
  int fftsize=0;

  fftsize=nextPower_2(2*N-1);
  int e = fftsize -1;
  cx_vec temp2;

  if (autoflag == true) {
    //Take FFT of input vector
    cx_vec X = cx_vec(x,zeros(x.size()));
    X= fft(X,fftsize);
    //cout<<"fft"<< X(span(0,3))<<endl;
    //Compute the abs(X).^2 and take the inverse FFT.
    temp2=ifft(X%conj(X));
    //cout<<"ifft"<< temp2(span(0,3))<<endl;
  }
  else
  {
   //Take FFT of input vectors
    cx_vec X=cx_vec(x,zeros(x.size()));
    cx_vec Y=cx_vec(y,zeros(y.size()));
    X = fft(X,fftsize);
    Y = fft(Y,fftsize);
    temp2 =ifft(X%conj(Y));
   }
   maxlag=N-1;
   return real(join_cols(temp2(span(e - maxlag + 1, e)),temp2(span(0,maxlag))));
}



void get_PL_SH(vec d, vec & PL, vec &  std_dev)
{
      vec PLfree;
      uvec ind;
      //PLfree.resize(d.size());
      PL=zeros(d.size());
      double i,p;
      double const fc = 5.91e9 ; // % Carrier frequency (Hz)
      double const  hBS = 1.5 ; //   % Transmitter antenna height (m)
      double const  hMS = 1.5 ;  //   % Receiver antenna height (m)
      double const  environmentHeight = 0 ;//   % Average environmental height (m)
      double const c = 3e8;
      vec dist=abs(d);

      //verifOutput(dist);
      double dBP = 4*(hBS-environmentHeight)*(hMS-environmentHeight)*fc/c ;
      // avoid errors for very small distances
      ind = arma::find(dist<3);
      for(int j=0; j<ind.size();j++)
      {
        //cout<<"index value= "<<ind(j);
        i=ind(j);
        dist(i)= 3 ;
        //cout<<"\n d= "<<d(i)<<" \n";
      }
      //dist.for_each( [](double& val) { if(val<3) return 3; } );  // NOTE: the '&' is crucial!
      //calculate pathloss for distances lower than the breakpoint distance

      ind= arma::find(dist<dBP);
      //cout<<"\n ind dBP = "<<ind <<" \n";
      for(int j=0; j<ind.size();j++)
      {
          i=ind(j);
          //cout<<"\n d tes = "<<d(i) <<" \n";
          //cout<<" i = "<<i <<" \n";
          p= 22.7*log10(dist(i))+27+20*log10(fc/1e9) ;
          PL(i)=p;
          //vec_push(PL,p);
          vec_push(std_dev,3);
       }


      // calculate pathloss for distances higher than the breakpoint distance
      ind=find(dist>=dBP);
       for(int j=0; j<ind.size();j++)
      {
          i=ind(j);
          p = 40*log10(dist(i))+7.56-17.3*log10(hBS-environmentHeight)-17.3*log10(hMS-environmentHeight)+2.7*log10(fc/1e9) ;

          PL(i)=p;
          //vec_push(PL,p);
          vec_push(std_dev,3) ;
       }


       //Compare obtained pathloss with free pathloss
       PLfree= 20*log10(dist)+46.4+20*log10(fc*1e-9/5) ;
       ind = find(PLfree>PL);


       for(int j=0; j<ind.size();j++)
       {
          i=ind(j);
          //cout<<" i = "<<i <<" \n";
          PL(i) = PLfree(i) ;
       }

}






void get_SINRdistribution(vec Pr_dBm_avg, double Pi_dBm_avg, vec std_dev_Pr, double std_dev_Pi, double noise_dBm, double sensingThreshold,
 double step_dB, mat& SINR, mat& PDF_SINR)
{
  vec  distrib_Pr,distrib_Pi_noise;
  vec x;
  uvec ind, aux;
  double noise;
  x = regspace(-200,step_dB,200) ;
  vec t=regspace(min(x)*2, step_dB, max(x)*2);

  for (int i=0; i<Pr_dBm_avg.size(); i++)
  {

    distrib_Pr = normpdf(x, Pr_dBm_avg(i), std_dev_Pr(i)) ;
    //cout << "\n distrib_Pr size " << distrib_Pr.size() << endl;
    ind=arma::find(x < sensingThreshold);
    //cout << "\n ind" << ind << endl;
    for(int i=0;i<ind.size();i++)
    {
     double j=ind(i);
     //cout << "\n j" << j<< endl;
     distrib_Pr(j) = 0 ;
     //cout << "\n distrib" << distrib_Pr(j)<< endl;
    }

    distrib_Pr = distrib_Pr/arma::sum(distrib_Pr)/step_dB ;
    //cout << "distrib tes " << distrib_Pr<< endl;
    distrib_Pi_noise.resize(x.size());
    if (Pi_dBm_avg==-datum::inf)
    {
      distrib_Pi_noise = arma::zeros<vec>(x.size()) ;
      double k=round((noise_dBm-x(1))/step_dB)+1;
      distrib_Pi_noise.at(k) = 1/step_dB ;
    }
    else
    {
      aux = find(x<=noise_dBm);
      double s=aux.size();
      noise = pow(10, (noise_dBm/10.0)) ;
      vec x1=x(span(s, x.size()-1));
      x1/=10;
      /// TODO: optimize this part
      vec t1=exp10(x1);
      vec t_aux=t1;
      t_aux-=noise;
      vec t_aux1=t_aux;
      t_aux1*=std_dev_Pi*sqrt(2*datum::pi);
      t_aux=log10(t_aux);
      t_aux*=10;
      t_aux-=Pi_dBm_avg;
      t_aux=arma::square(t_aux);
      t_aux*=-1;
      t_aux/=(2*pow(std_dev_Pi, 2));
      t_aux=arma::exp(t_aux);
      vec t2=t1/t_aux1;
      ///
      vec distrib_Pi_noise_aux =t2%t_aux;
      distrib_Pi_noise=join_cols(zeros<vec>(s),distrib_Pi_noise_aux);
      distrib_Pi_noise = distrib_Pi_noise/arma::sum(distrib_Pi_noise)/step_dB;
     }
}
     rowvec temp, temp1;
     temp.resize(t.size());
     temp1.resize(t.size());
     vec vecInter=xcorr(distrib_Pr, distrib_Pi_noise,false);
     for(int j=0; j<t.size();j++)
     {
        temp(j)=t(j);
        temp1(j)=vecInter(j);
     }

     for (int i=0; i<Pr_dBm_avg.size(); i++)
     {
     SINR.insert_rows(i,temp);
     PDF_SINR.insert_rows(i,temp1);
     PDF_SINR.row(i)=PDF_SINR.row(i)/sum(PDF_SINR.row(i))/step_dB;
     }
}


void recordOutput(vec var)
{
  ofstream myfile;
  myfile.open("output.txt");
  myfile << var <<endl;
  myfile.close();

}

vec get_BLER(mat SINR, mat PDF, double coding, double step_dB)
{
  vec BLER_interp, vector_BLER_paper, vector_SNR_paper ;
  vec avg_BLER;
  double a,b;
  a = min(-3.0, SINR(1)) ;
  b = max(21.0, SINR(SINR.size()-1)) ;

  if (coding==1)
  {
    // % 190 Bytes, QPSK r=0.7, Vr = 280 km/h - From R1-160284, DMRS enhancement of V2V in 3GPP
    vector_SNR_paper = {a, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 20+1e-6, b} ; ;
    vector_BLER_paper = {1, 1, 0.9, 0.7, 0.4, 0.13, 0.045, 0.017, 0.007, 1e-3, 1e-3, 1e-3, 1e-4, 1e-4} ;
  }
  else if (coding==2)
  {
    // % 190 Bytes, QPSK r=0.5, Vr = 280 km/h - From R1-160284, DMRS enhancement of V2V in 3GPP
    vector_SNR_paper = {a, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 20+1e-6, b} ;
    vector_BLER_paper = {1, 1, 0.9, 0.7, 0.3, 0.09, 0.02, 0.002, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-4, 1e-4} ;
  }else
  {
    cout<<" Coding is not valid " <<endl;
  }
    avg_BLER.resize(SINR.n_rows);
    for (int i=0; i<SINR.n_rows; i++)
    {
     rowvec rs=SINR.row(i);
     vec sv=conv_to<vec>::from(rs);
     arma::interp1(vector_SNR_paper, vector_BLER_paper, sv, BLER_interp, "linear") ;
     rowvec pdfrs=PDF.row(i);
     vec pdfv=conv_to<vec>::from(pdfrs);
     vec tes= pdfrs*BLER_interp*step_dB;
     avg_BLER(i)=sum(tes);
    }
     return avg_BLER ;
}

// Verify with matlab example
void CV2XMode4_common(double lambda, double Pt, vec distance, double Psen, double step_dB,
 double noise, double coding, vec & deltaHD, vec & deltaSEN, vec & deltaPRO)
{
  vec PL_E_R, std_dev_E_R ;
  mat PDF_SNR, SNR;
  int D=distance.size();
  deltaHD.resize(D);
  for(int i=0; i<D; i++)
  {
  deltaHD(i)=lambda/1000.0;
  }
  get_PL_SH(distance, PL_E_R, std_dev_E_R) ;
  deltaSEN = 0.5*(1-erf((Pt-PL_E_R-Psen)/(std_dev_E_R*sqrt(2)))) ;
  get_SINRdistribution(Pt-PL_E_R, -180, std_dev_E_R, 3, noise, Psen, step_dB, SNR, PDF_SNR) ;
  deltaPRO=get_BLER(SNR, PDF_SNR, coding, step_dB);
}



// verify with  matlab example
void CV2XMode4_Step2(double beta, double lambda, double Pt, double S, vec distance, double Psen, double step_dB,
    double noise, double coding, vec deltaPRO, vec& deltaCOL, double& CBR)
{
    vec PL_E_R,std_dev_E_R,PL_I_R, std_dev_I_R,PL_I_E, std_dev_I_E,p_int,
    pSsensing_I_E,p_s, p_sim_2,p_col,
    p_SINR,distance_int_to_rx,PSR,d_aux,PL,std_dev,R_PSR,Ce,Ca,Cc,distance_int_to_tx;
    mat SINR, PDF_SINR;

    double Pi_dB,Spsr,Ne,Na,Lint_max,R0,N,Nc;
    N = S*1000/lambda;
    Nc = 0.2*N ;
    d_aux = regspace(-1500, 1, 1500) ;
    get_PL_SH(d_aux, PL, std_dev) ;
    PSR = 0.5*(1+erf((Pt-PL-Psen)/(std_dev*sqrt(2)))) ; // OK
    //cout<<"PSR size " << PSR.size() <<endl; // OK
    //cout<<"PSR " << PSR(span(0,4)) <<endl; // OK

    Spsr = beta*sum(PSR) ;
    Ne = Spsr/2.0+arma::sum(1-(regspace(1, 1, (Spsr/2.0)))/(N-Spsr/2.0)) ;
    Na = N-Ne ;
    CBR = Ne/N ; // OK
    Lint_max = round(1000*beta)/beta ;
    distance_int_to_rx =regspace(-1*Lint_max,1/beta,Lint_max);
    distance_int_to_rx((distance_int_to_rx.size()+1)/2.0)={} ;
    R_PSR =  xcorr(PSR,PSR,true);
    R_PSR = R_PSR(span(2*max(d_aux)+1, R_PSR.size()-1)) ;
    R0 = R_PSR(1) ;
    Ce = (R_PSR/R0)*(beta*Ne*R0/Spsr-Ne*Ne/N)+Ne*Ne/N ;
    Ca = N-2*Ne+Ce ;
    Cc = Ca*pow((Nc/Na), 2) ;

    double tau = lambda ;
    get_PL_SH(distance, PL_E_R, std_dev_E_R) ;

    //vec dist1=distance;
    //dist1.resize(distance_int_to_rx.size());
    //distance_int_to_tx=distance_int_to_rx+dist1;

    get_PL_SH(distance_int_to_rx, PL_I_R, std_dev_I_R) ;

    distance_int_to_tx=regspace(-1000,10,1000); // replaced to resolve bug//
    get_PL_SH(distance_int_to_tx, PL_I_E, std_dev_I_E) ; //

    pSsensing_I_E = 0.5*(1+erf((Pt-PL_I_E-Psen)/(pow(2, 0.5)*std_dev_I_E))) ;
    p_s = 1-(1-(1/tau))*pSsensing_I_E ; // Bug rounding mode to check later


    p_sim_2=zeros(p_s.size());

    //cout<< "loop1"<<endl;
    for (int i=0; i<distance_int_to_rx.size(); i++)
    {
    //cout<<"loop 2"<<endl;
      Pi_dB = Pt-PL_I_R(i) ;
      //cout<<"Pi_dB " << Pi_dB <<endl; // OK
      get_SINRdistribution(Pt-PL_E_R, Pi_dB, std_dev_E_R, std_dev_I_R(0), noise, Psen, step_dB, SINR, PDF_SINR) ;
      double temp=p_s(i)*Cc(round(abs(distance_int_to_tx(i)))+1)/(Nc*Nc) ;
      p_sim_2(i)=temp;
    }
      p_SINR= get_BLER(SINR, PDF_SINR, coding, step_dB) ; // bug PDF,
     // cout<< "p_sinr " <<p_SINR<<endl;

     for (int d=0; d<distance.size(); d++)
     {
      if (deltaPRO(d)==1)
      {
       vec_push(p_int,0);
      }
      else
      {
         double tes= (p_SINR(d)-deltaPRO(d))/(1-deltaPRO(d));
         vec_push(p_int,tes);
      }

       double tem1=p_sim_2(d)*p_int(d);
       vec_push(deltaCOL,1-prod(1-tem1));
      }

}






vec CV2XMode4_Step3(double beta, double lambda, double Pt, double S, vec distance, double Psen, double step_dB,
 double noise, double coding, vec deltaPRO)
{
    vec PL,std_dev,d_aux,PSRn,distance_int_to_rx, PSR,R_PSR,
    Ce,Cc,Ca,  PL_E_R, std_dev_E_R, distance_int_to_tx,
    PL_I_R, std_dev_I_R,PL_I_E, std_dev_I_E,p_int,p_SINR,
    pSensar_I_E, p_s,p_sim_3,p_col,deltaCOL;

    mat SINR, PDF_SINR;
    double Spsr_n, Ne, Na,N,Nc,Psen_Step3,Lint_max,Spsr,R0,Pi_dB ;

    N = S*1000*1.0/lambda ;
    Nc = 0.2*N ;
    Psen_Step3 = Psen-step_dB ;
    Na = 0 ;
    d_aux = regspace(-1500, 1, 1500) ;
    get_PL_SH(d_aux, PL, std_dev) ;
    double s=0;
    while (Na<Nc)
    {
    //cout<<"while "<< endl;
    Psen_Step3 += step_dB ;
    PSRn = 0.5*(1+erf((Pt-PL-Psen_Step3)/(std_dev*sqrt(2)))) ;
    Spsr_n = 2*beta*sum(PSRn) ;
    Ne = Spsr_n/2.0+arma::sum(1-(regspace(1, 1, (Spsr_n/2.0)))/(N-Spsr_n/2.0)) ;
    Na = N-Ne ;
    }
    //cout<<"after while" <<endl;
    Lint_max = round(1000*beta)/beta ;
    distance_int_to_rx =regspace(-Lint_max,1/beta,Lint_max);
    uvec x=find(abs(distance_int_to_rx)<1e-6);

    for(int i=0;i<x.size();i++)
    {
    distance_int_to_rx(x(i))={};
    }

    PSR = 0.5*(1+erf((Pt-PL-Psen)/(std_dev*sqrt(2)))) ;
    R_PSR = xcorr(PSR,PSR,true) ;
    R_PSR = R_PSR(span(2*max(d_aux), R_PSR.size()-1)) ;
    Spsr = beta*sum(PSR) ;
    R0 = R_PSR(1) ;
    Ce = (R_PSR/R0)*(beta*Ne*R0/Spsr-Ne*Ne/N)+Ne*Ne/N ;
    Ca = N-2*Ne+Ce ;
    Cc = Ca*pow((Nc/Na), 2) ;

    get_PL_SH(distance, PL_E_R, std_dev_E_R) ;
    //vec dist1=distance;
    //dist1.resize(distance_int_to_rx.size());
    //distance_int_to_tx=distance_int_to_rx+dist1;

    get_PL_SH(distance_int_to_rx, PL_I_R, std_dev_I_R) ;
    distance_int_to_tx=regspace(-1000,10,1000); // replaced to resolve bug//
    get_PL_SH(distance_int_to_tx, PL_I_E, std_dev_I_E) ;


    double tau = lambda ;
    pSensar_I_E = 0.5*(1+erf((Pt-PL_I_E-Psen_Step3)/(pow(2, 0.5)*std_dev_I_E))) ;
    p_s = 1-(1-1*1.0/tau)*pSensar_I_E ;

    p_sim_3=zeros(p_s.size());


    for (int i=0; i<distance_int_to_rx.size(); i++)
    {
       Pi_dB = Pt-PL_I_R(i) ;
       get_SINRdistribution(Pt-PL_E_R, Pi_dB, std_dev_E_R, std_dev_I_R(0), noise, Psen, step_dB, SINR, PDF_SINR) ;
       double temp=p_s(i)*Cc(round(abs(distance_int_to_tx(i)))+1)/(Nc*Nc) ;
       vec_push(p_sim_3,temp);
    }
     p_SINR= get_BLER(SINR, PDF_SINR, coding, step_dB) ;
     //p_int.resize(p_SINR.size());
    for (int d=0; d<distance.size(); d++)
    {
      if (deltaPRO(d)==1)
      {
        vec_push(p_int,0);
      }
      else
      {
        // TO fix the problem is the expression below:
         double tes= (p_SINR(d)-deltaPRO(d))/(1-deltaPRO(d));
         vec_push(p_int,tes);
      }

      double tem1=p_sim_3(d)*p_int(d);
      vec_push(deltaCOL,1-prod(1-tem1));

     }

  return deltaCOL ;
}






void CV2XMode4Principal(double beta, double lamda, double Pt, double S, double B,
 vec distance, vec& PDR, vec & deltaHD, vec & deltaSEN, vec & deltaPRO,vec& deltaCOL, double& CBR)
{

    cout<<"principal " <<endl;
    double Psen=-90.5;
    double step_dB=0.1;
    double coding;
    double RBs;

        //if(B==190){
         //   {
    /// this is true for B=190 B
            if(S==4)
            {
                coding=1;
                RBs=10;

            }
            else if(S==2)
            {
                coding=2;
                RBs=12;
            }

           // }
           // }




    double noise=-95-10*log10(50/RBs);

    vec deltaHD_pre, deltaSEN_pre, deltaCOL2_pre, deltaCOL3_pre,deltaPRO_pre;
    double alpha;

    CV2XMode4_common(lamda,Pt, distance, Psen, step_dB, noise, coding,deltaHD_pre,deltaSEN_pre,deltaPRO_pre);


    CV2XMode4_Step2(beta, lamda,  Pt,  S, distance,  Psen, step_dB,
        noise,  coding, deltaPRO_pre, deltaCOL2_pre, CBR);


    if(CBR<0.2)
    {
        alpha=0;

    }
    else if(CBR<=0.7)
    {
        alpha=2*CBR-0.4;
    }
    else
    {
        alpha=1;
    }

    if(alpha<1)
    {
        deltaCOL3_pre=CV2XMode4_Step3(beta,lamda, Pt, S, distance, Psen, step_dB,
      noise, coding, deltaPRO_pre);
    }else
    {
     deltaCOL3_pre=zeros(distance.size());
    }

        //% Calculate final probabilities for each type of error:
        deltaHD   = deltaHD_pre;                              // % Equation (6.1)
        deltaSEN  = deltaSEN_pre % (1 - deltaHD);             //% Equation (6.2)

        deltaPRO  = deltaPRO_pre % (1 - deltaHD_pre) % (1 - deltaSEN_pre); //% Equation (6.3)

        vec deltaCOL2 = deltaCOL2_pre % (1 - deltaHD_pre) % (1 - deltaSEN_pre) % (1 - deltaPRO_pre); //% Equation (6.4)
        vec deltaCOL3 = deltaCOL3_pre % (1 - deltaHD_pre) % (1 - deltaSEN_pre) % (1 - deltaPRO_pre); //% Equation (6.5)
        deltaCOL = alpha*deltaCOL2 + (1-alpha)*deltaCOL3; //% Equation (21)
        //% Calculate PDR:
        PDR = 1 - deltaHD - deltaSEN - deltaPRO - deltaCOL; //% Equation (6)
        cout<<" PDR= " <<  PDR <<endl;
        recordOutput(PDR);
    }

double estimateDelay(double PDR, double lamda, double Ttr)
{

	double PL=1-PDR;
	return Ttr*(0.5+(lamda*PL/(1-PL)));


}

vec calculateEWA(vec dist, double beta, double newSample)
{

	for (int i=1; i<dist.size(); i++)
	{
		dist(i)=dist(i-1)*beta +(1-beta)*newSample;
	}

	return dist;
}



};
