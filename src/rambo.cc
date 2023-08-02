#include "../include/rambo.h"
#include "TMath.h"
#include <iostream>
#include <string>
#include <stdexcept>

using namespace std;



Rambo::Rambo(int n, TLorentzVector parent, std::vector<double> masses, int seed){
    num_daughter = n;
    m_parentMomenta = parent; 
    m_daughterMassContainer.clear();
    if(num_daughter != masses.size()) throw std::invalid_argument("number of daughters and length of the mass vector do not match");
    //else m_daughterMassContainer = masses;
    for(int i = 0; i < num_daughter; i++ ) m_daughterMassContainer.push_back(masses[i]);
    rand.SetSeed(seed);
}

Rambo::~Rambo(){};

void Rambo::Generate(){
    bool accept = false;
    do{
        this->m_daughterMomenta_random.resize(0);
        this->m_daughterMomenta_massless.clear();
        this->m_daughterMomenta_massive.clear();
        this->GenerateQ();
        this->GenerateP();
        this->GenerateMass();
        //this->CalculateWeightsMassless();
        this->CalculateWeightsMassive();
        this->CalculateAcceptance();
        //cout<<"Acceptance: "<<(this->acceptance)<<endl;
        if( rand.Uniform(0,1) < (this->acceptance)){ 
            accept = true;
            //cout<<"Accepted"<<endl;
        }
        //else cout<<"Rejected"<<endl;
    }
    while(!accept);
    this->BoostDaughter();
}

TLorentzVector * Rambo::GetDecay(int i){
    if(i < this->num_daughter) return &(this->m_daughterMomenta_massive[i]);
    else{
        
        throw std::invalid_argument("Error: the requested index "+ std::to_string(i)+" goes out of the scope 0-"+std::to_string(this->num_daughter));
    }
}

void Rambo::GenerateQ(){
    m_daughterMomenta_random.resize(num_daughter);
    for(int i = 0; i < num_daughter; i++){
        double rho1 = rand.Uniform(0,1);
        double rho2 = rand.Uniform(0,1);
        double rho3 = rand.Uniform(0,1);
        double rho4 = rand.Uniform(0,1);
        double c_i = 2*rho1 - 1;
        double phi_i = 2*TMath::Pi()*rho2;
        double q0 = -TMath::Log(rho3*rho4);
        double qx = q0 * sqrt(1 - c_i * c_i) * TMath::Cos(phi_i);
        double qy = q0 * sqrt(1 - c_i * c_i) * TMath::Sin(phi_i);
        double qz = q0 * c_i;

        TLorentzVector Qi;
        Qi.SetPxPyPzE(qx, qy, qz, q0);
        m_daughterMomenta_random[i] = Qi;
    }
}

void Rambo::GenerateP(){
    TLorentzVector Q_mu = m_daughterMomenta_random.sum();
    double M = Q_mu.M();
    TVector3 b = -Q_mu.Vect()*(1/M);
    double gamma = Q_mu.E()/M;
    double a = 1/( 1 + gamma );
    double x = m_parentMomenta.M()/M;
    m_daughterMomenta_massless.clear();

    for(int i = 0; i < num_daughter; i++){
        double qiE = m_daughterMomenta_random[i].E();
        TVector3 qi_3Vec = m_daughterMomenta_random[i].Vect();
        double b_dot_q = b.Dot(qi_3Vec);
        double p0 = x * (gamma * qiE + b_dot_q);
        TVector3 pi_3Vec = x * ( qi_3Vec + b * qiE + a * ( b.Dot(qi_3Vec) ) * b );
        TLorentzVector aMomenta(pi_3Vec, p0);
        m_daughterMomenta_massless.push_back(aMomenta);
    }
}

void Rambo::GenerateMass(){
    //if(m_daughterMomenta_massless.size() == m_daughterMassContainer.size()){
    //    const int ContSize = m_daughterMassContainer.size();
    //    this->m_daughterMomenta_massive.resize(ContSize);
    //}

    this->m_daughterMomenta_massive.resize(num_daughter);
    const int iterMax = 1000;
    const double Accuracy = 1e-7;

    std::valarray<double> E(num_daughter);
    std::vector<double> XM2(num_daughter);
    std::valarray<double> P2(num_daughter);

    for(int i = 0; i < num_daughter; i++){
        XM2[i] = m_daughterMassContainer[i] * m_daughterMassContainer[i];
        P2[i] = m_daughterMomenta_massless[i].E() * m_daughterMomenta_massless[i].E();
    }

    double XMT = 0;
    double ET = m_parentMomenta.M();

    for(int i = 0; i < num_daughter; i++){
        XMT = XMT + m_daughterMassContainer[i];
    }

    double XMAX = std::sqrt(1 - (XMT/ET) * (XMT/ET) );
    double X = XMAX;
    int n_iter = 0;

    while(n_iter < iterMax){
        double F0 = -ET;
        double G0 = 0;
        double X2 = X*X;

        for(int i = 0; i < num_daughter; i++){
            E[i] = std::sqrt(XM2[i] + X2 * P2[i]);
            G0 = G0 + P2[i] / E[i];
        }

        F0 = F0 + E.sum();

        if(fabs(F0) <= Accuracy){
           for(int i = 0; i < num_daughter; i++){
               TLorentzVector massive(X*m_daughterMomenta_massless[i].Vect(),E[i]);
               this->m_daughterMomenta_massive[i] = massive;
           }
           break;
        }
        else{
            n_iter++;
            X = X - F0/(X*G0);
        }
        if (n_iter == iterMax){
            std::cout<<"iterations: "<<n_iter<<endl;
            std::cout<<"RAMBO DID not CONVERGE"<<endl;
        }
    }
}

int Rambo::fak_n(int num){
    if(num > 1) return fak_n(num-1)*num;
    return 1;
}

void Rambo::CalculateWeightsMassless(){
    double w = m_parentMomenta.M();
    double prod1 = pow(TMath::Pi() / 2, num_daughter - 1);
    double prod2 = pow(w, 2* num_daughter - 4);
    double W_0 = (prod1 * prod2) / ((this->fak_n(num_daughter - 1)) * (this->fak_n(num_daughter - 2)));
    m_weightMassless = W_0;
}

void Rambo::CalculateWeightsMassive(){
    double w = m_parentMomenta.M();
    double sum1 = 0;
    double sum1a = 0;
    double sum1b = 0;
    double sum2 = 0;
    double sum2a = 0;
    double prod1 = 1;
    for( int i = 0; i < num_daughter; i++ ){
        double k_mag = m_daughterMomenta_massive[i].Vect().Mag();
        sum1a = sum1a + k_mag;

        double k_mag2 = m_daughterMomenta_massive[i].Vect().Mag2();
        double k_i0 = m_daughterMomenta_massive[i].E();
        double k_ratio2 = k_mag2/k_i0;
        sum2a = sum2a + k_ratio2;

        double k_ratio = k_mag/k_i0;
        prod1 = prod1 * k_ratio;
    }

    sum1b = (1/w) * sum1a;
    sum1 = pow(sum1b, 2 * num_daughter - 3);
    sum2 = pow(sum2a, -1);

    double W_m = sum1 * sum2 * prod1;

    this->m_weightMassive = W_m;
}

void Rambo::CalculateAcceptance(){
  this->acceptance = (this->m_weightMassive) * m_parentMomenta.M();
}

void Rambo::BoostDaughter(){
    TVector3 v_boost = m_parentMomenta.BoostVector();
    for(int i = 0; i < num_daughter; i++){
        m_daughterMomenta_massive[i].Boost((v_boost));
    }
}
