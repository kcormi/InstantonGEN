#include <iostream>
#include "../include/linearInterpolation.h"

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


