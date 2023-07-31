#include "../include/configBuilder.h"
#include <iostream>

using namespace std;

float expectNg(float Q){
  return 6.02;
}

int sampleNg(float Q2,TRandom3 rand){
  float Q = sqrt(Q2);
  return rand.Poisson(6.02);
}

configBuilder::configBuilder()
{
    //Seed random number generator
    rand.SetSeed(0);

    //Build Fermion Doublets, the order of them is SUPER important
    particle partBuf;
    partBuf.mass = DQ_MASS;
    partBuf.pid = DQ_PID;
    partBuf.color[0] = 1;
    partBuf.color[1] = 0;
    partBuf.q3 = DQ_Q3;
    parts.push_back(partBuf);  // 2

    partBuf.mass = UQ_MASS;
    partBuf.pid = UQ_PID;
    partBuf.color[0] = 1;
    partBuf.color[1] = 0;
    partBuf.q3 = UQ_Q3;
    parts.push_back(partBuf);  // 3

    partBuf.mass = SQ_MASS;
    partBuf.pid = SQ_PID;
    partBuf.color[0] = 1;
    partBuf.color[1] = 0;
    partBuf.q3 = SQ_Q3;
    parts.push_back(partBuf);  // 6

    partBuf.mass = CQ_MASS;
    partBuf.pid = CQ_PID;
    partBuf.color[0] = 1;
    partBuf.color[1] = 0;
    partBuf.q3 = CQ_Q3;
    parts.push_back(partBuf);  // 7

    partBuf.mass = BQ_MASS;
    partBuf.pid = BQ_PID;
    partBuf.color[0] = 1;
    partBuf.color[1] = 0;
    partBuf.q3 = BQ_Q3;
    parts.push_back(partBuf);  // 10

    partBuf.mass = TQ_MASS;
    partBuf.pid = TQ_PID;
    partBuf.color[0] = 1;
    partBuf.color[1] = 0;
    partBuf.q3 = TQ_Q3;
    parts.push_back(partBuf);  // 11

    partBuf.mass = GLU_MASS;
    partBuf.pid = GLU_PID;
    partBuf.color[0] = 1;
    partBuf.color[1] = 1;
    partBuf.q3 = GLU_Q3;
    parts.push_back(partBuf);  // 11

}

configBuilder::~configBuilder()
{
}

vector<particle> configBuilder::build(int iQ1, int iQ2, float Q2, int seed, int Nf)
{
    vector<particle> conf;
    if( Nf > 6 || Nf < 1){
        cout<<"Nf = "<< Nf<<", not accepted. Nf choice: 1, 2, 3, 4, 5, 6"<<endl;
        return conf;
    }

    if(iQ1 != 21 || iQ2 != 21){
        cout<<"Inital state pid "<<iQ1<<", "<<iQ2<<" not accepted. Currently I can only generate instantons from gluon inital states."<<endl;
        return conf;
    }


    for(int i = 0; i < Nf; i++){
        particle part = parts[i];
        conf.push_back(part);
        part.pid = -part.pid;
        part.color[0] = 0;
        part.color[1] = 1;
        conf.push_back(part);
    }
    rand.SetSeed(seed);
    for(int j = 0; j < sampleNg(Q2,rand); j++){
        particle part = parts[6];
        conf.push_back(part);
    }

    return conf;
}






