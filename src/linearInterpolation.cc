#include <iostream>
#include "../include/linearInterpolation.h"
#include "TMath.h"
double Lerp(double a, double b, double t)
{
    return a + t * (b - a);
}

int searchInterval(double energy, const double energylist[], int len){
    if(energy < energylist[0]) return -1;
    for(int i = 0; i < (len-1); i++){
      if(energy < energylist[i+1]) return i;
    }
    return len - 1;
}

double getLerp(double energy, const double energylist[], const double values[], int len){
    int interval = searchInterval(energy, energylist, len);
    double begin;
    double end;
    double value_begin;
    double value_end;
    if(interval == -1) {
        begin = energylist[0];
        end = energylist[1];
        value_begin = values[0];
        value_end = values[1];
    }
    else if (interval == (len-1)){
        begin = energylist[len-2];
        end = energylist[len-1];
        value_begin = values[len-2];
        value_end = values[len-1];
    }
    else{
        begin = energylist[interval];
        end = energylist[interval+1];
        value_begin = values[interval];
        value_end = values[interval+1];
    }
    return Lerp(value_begin, value_end, (energy-begin)/(end-begin));
}

void getIntegral(double Integral[], const double energylist[], const double values[], int len){
     Integral[0] = 0;
     for(int i = 0; i < len-1; i++) Integral[i+1] = Integral[i] + (values[i] + values[i+1]) * (energylist[i+1] - energylist[i]) / 2;
}

void getCDF(double CDF[], const double energylist[], const double values[], int len){
     getIntegral(CDF, energylist, values, len);
     for(int i = 0; i < len-1; i++) CDF[i+1] = CDF[i+1] / CDF[len-1]; 
}

void getPDF(double PDF[], const double energylist[], const double values[], int len){
    double sum = 0;
    for(int i = 0; i < len-1; i++) sum += (values[i] + values[i+1]) * (energylist[i+1] - energylist[i]) / 2;
    for(int i = 0; i < len; i++) PDF[i] = values[i]/sum;
}

double getInterpoCDF(double energy,const double energylist[], const double CDF[], const double PDF[], int len){
    int interval = searchInterval(energy, energylist, len);
    if(interval == -1) return 0;
    else if(interval == (len-1)) return 1;
    else{
        double slope = (PDF[interval+1] - PDF[interval]) / (energylist[interval+1] - energylist[interval]);
        return (2 * PDF[interval] + slope * (energy - energylist[interval])) * (energy - energylist[interval]) / 2 + CDF[interval];
    }

}

double invertCDF(double valCDF, const double energylist[], const double CDF[], const double PDF[], int len){
    int interval = searchInterval(valCDF, CDF, len);
    if(interval == -1) return energylist[0];
    else if(interval == (len-1)) return energylist[len-1];
    else{
        //std::cout<<"valCDF: "<<valCDF<<std::endl;
        //std::cout<<"interval: "<<interval<<std::endl;
        //std::cout<<"CDF[interval]: "<<CDF[interval]<<std::endl;
        double slope = (PDF[interval+1] - PDF[interval]) / (energylist[interval+1] - energylist[interval]);
        //std::cout<<"slope: "<<slope<<std::endl;
        //std::cout<<"PDF[interval] * PDF[interval] + 8 * slope * (valCDF - CDF[interval]) : "<<PDF[interval] * PDF[interval] + 8 * slope * (valCDF - CDF[interval])<<std::endl;
        //std::cout<<"sqrt: "<<TMath::Sqrt(PDF[interval] * PDF[interval] + 8 * slope * (valCDF - CDF[interval]))<<std::endl;
        return energylist[interval] + (TMath::Sqrt(4 * PDF[interval] * PDF[interval] + 8 * slope * (valCDF - CDF[interval])) - 2 * PDF[interval]) / (2 * slope);
    }

}


