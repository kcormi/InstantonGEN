#include <iostream>
#include <stdlib.h> 
#include <cmath>
#include <vector>
#include <string>

#include "LHAPDF/LHAPDF.h"

#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TString.h"

#include "../include/configBuilder.h"
#include "../include/InstantonGEN.h"
#include "../include/LHEWriter.h"
#include "../include/rambo.h"
#define ARGS 6
#define MPW 5.0e-2

using namespace LHAPDF;
using namespace std;

int main(int argc, char* argv[])
{

    //read arguements and do conversions and checks on them
    if(argc != ARGS+1)
    {
        cout << "InstantonGEN sqrtS threshold maxweight Nevents Nf Filename" << endl;
        return -99;
    }

    float SQRTS = atof(argv[1]);
    float thr = atof(argv[2]);
    float MCW = atof(argv[3]);
    if(MCW < 0.0)
    {
        cout << "maxweight must be greater than zero" << endl;
        return -99;
    }
    int Nevt = atoi(argv[4]);
    int Nf = atoi(argv[5]);
    if ( Nf > 6 || Nf < 1){
        cout << "Nf must be between 1 and 6" << endl;
        return -99;
    }

    string ofName = string(argv[6]);

    //Init structures needed during generatation
    double maxwt = 0;
    double maxMCtot = 0.0;

    int NF = 0;
    TRandom3 rand;
    rand.SetSeed(0);
    particleBase *partBase = new particleBase();
    configBuilder confBuild;
    vector<double> daughtersPt;
    vector<double> daughtersEta;
    vector<double> daughtersPhi;
    vector<double> daughtersE;
    vector<double> daughtersID;
    TLorentzVector mom(0.0,0.0,0.0,thr);
    TLorentzVector u1(0.0,0.0,0.0,0.0);
    TLorentzVector u2(0.0,0.0,0.0,0.0);
    vector<double> masses;

    double weight = 0.0;

    //TGenPhaseSpace gen;
    double x1 = 0.0;
    double x2 = 0.0;
    int iq1 = 0;
    int iq2 = 0;
    double Q2 = 0.0;
    double momM = 0.0;
    double pz = 0.0;

    //Init output files
    TFile *myF = new TFile((ofName+".root").c_str(),"RECREATE","Holds daughters from sphaleron decay");
    LHEWriter lheF(ofName, SQRTS);

    double minx = thr*thr/(SQRTS*SQRTS);

    //Init parton distribution functions
    const PDF* LHApdf = mkPDF("CT10",0);
    //cout << LHApdf->xfxQ2(2, 0.5, SQRTS*SQRTS) << endl;

    //Initialize histograms for debugging
    TH1D *x1_h = new TH1D("x1_h","x1 inclusive",1000,0.0,1.0);
    TH1D *mcTot_h = new TH1D("mcTot_h","Monte Carlo Probabilities",100,0.0,MCW);
    //TH1D *sumInterQ3_h = new TH1D("sumInterQ3_h","Intermediate particle charges",21,-10.5,10.5);
    TH1D *sphM_h = new TH1D("sphM_h","Sphaleron Transition Energy;Invariant Mass [TeV];Events / 100 GeV",int(SQRTS/100.0),0.0,SQRTS/1000.0);
    TH1D *st_h = new TH1D("st_h","S_{T};S_{T} [TeV];Events / 100 GeV", int(SQRTS/100.0), 0.0, SQRTS/1000.0);
    TH1D *sphPz_h = new TH1D("sphPz_h","Sphaleron p_{z};p_{z} [GeV];Events / 100 GeV",80,-4000.0,4000.0);
    TH1D *outID_h = new TH1D("outID_h","Outgoing PDG IDs;PDG ID;Entries",33,-16.5,16.5);
    TH1D *p1x_h = new TH1D("p1x_h","Parton 1 Momentum Fraction;x_{1};Events / 0.01",60,0.4,1.0);
    TH1D *p1id_h = new TH1D("p1id_h","Parton 1 PDG ID;PDG ID;Events",13,-6.5,6.5);

    TH2D *inQid_h = new TH2D("inQid_h","Colliding Parton Species;Parton 2 PDG ID;Parton 1 PDG ID",11,-5.5,5.5,11,-5.5,5.5);
    TH2D *frac2D_h = new TH2D("frac2D_h","Parton Momentum Fractions;x_{2};x_{1}",60,0.4,1.0,60,0.4,1.0);

    vector<TH1D> pt_hv;
    vector<TH1D> eta_hv;
    vector<TH1D> phi_hv;
    for(int iq = -6; iq < 7; iq++)
    {
        string nameBuf;
        string titBuf;
        titBuf = Form("Quark %i p_{T};p_{T} [GeV];Entries / 20 GeV",iq);
        if(iq < 0) nameBuf = Form("pt_iqm%i_h",int(fabs(iq)));
        else nameBuf = Form("pt_iqp%i_h",int(fabs(iq)));
        TH1D hBufpt(nameBuf.c_str(),titBuf.c_str(),250,0.0,5000.0);
        pt_hv.push_back(hBufpt);

        titBuf = Form("Quark %i #eta;#eta;Entries / 0.1",iq);
        if(iq < 0) nameBuf = Form("eta_iqm%i_h",int(fabs(iq)));
        else nameBuf = Form("eta_iqp%i_h",int(fabs(iq)));
        TH1D hBufeta(nameBuf.c_str(),titBuf.c_str(),100,-5,5);
        eta_hv.push_back(hBufeta);

        titBuf = Form("Quark %i #phi;#phi;Entries / 0.1 #pi",iq);
        if(iq < 0) nameBuf = Form("phi_iqm%i_h",int(fabs(iq)));
        else nameBuf = Form("phi_iqp%i_h",int(fabs(iq)));
        TH1D hBufphi(nameBuf.c_str(),titBuf.c_str(),20,-1.0*TMath::Pi(),TMath::Pi());
        phi_hv.push_back(hBufphi);
    }

    //Initialize TTree to save event info
    TTree *myT = new TTree("mcTree","mcTree");
    myT->Branch("daughtersPt",&daughtersPt); 
    myT->Branch("daughtersEta",&daughtersEta); 
    myT->Branch("daughtersPhi",&daughtersPhi); 
    myT->Branch("daughtersE",&daughtersE); 
    myT->Branch("daughtersID",&daughtersID); 
    myT->Branch("weight",&weight); 
    myT->Branch("pz",&pz); 
    myT->Branch("x1",&x1); 
    myT->Branch("x2",&x2); 
    myT->Branch("iq1",&iq1); 
    myT->Branch("iq2",&iq2); 
    myT->Branch("Q2",&Q2); 
    myT->Branch("momM",&momM); 

    int pdfN = 0;

    while(NF < Nevt)//Main event generation loop
    {
        //Make sure everything is reset from the last event
        daughtersPt.clear();
        daughtersEta.clear();
        daughtersPhi.clear();
        daughtersE.clear();
        daughtersID.clear();

        Q2 = 0.0;
        double mcP = 1.1;
        bool mcPass = false;
        //Collide protons
        while(!mcPass)
        {
            //Choose x1 and x2 in proper range
            Q2 = 0.0;
            while(Q2 < thr*thr)
            {
                x1 = (1.0-minx)*rand.Uniform()+minx;
                x2 = (1.0-minx)*rand.Uniform()+minx;
                Q2 = x1*x2*SQRTS*SQRTS;
            }
            mcP = MCW*rand.Uniform();
            double mcTot = 0.0;
            for(int i1 = 21; i1 < 22; i1++)
            {
                if(i1 == 0) continue;
                iq1 = i1;
                double x1p = LHApdf->xfxQ2(iq1, x1, Q2)/x1;
                for(int i2 = 21; i2 < 22; i2++)
                {
                    if(i2 == 0) continue;
                    iq2 = i2;
                    //double x2p = pdf->parton(iq2,x2,SQRTS);
                    double x2p = LHApdf->xfxQ2(iq2, x2, Q2)/x2;
                    mcTot += x1p*x2p;
                    if(mcTot > mcP) {mcPass = true; break;}
                }
                if(mcPass) break;
            }
            pdfN++;
            x1_h->Fill(x1);
            mcTot_h->Fill(mcTot);
            if(maxMCtot < mcTot) {maxMCtot = mcTot; cout << "Max MC Total: " << maxMCtot << endl;}
        }

        //Build incoming particles, sphaleron, and prepare to decay
        particle partBuf1 = partBase->getParticle(iq1);
        if(iq1 != partBuf1.pid) cout << "iq1 = " << iq1 << " != partBuf.pid = " << partBuf1.pid << endl;
        partBuf1.p4v.SetXYZM(0.0,0.0,x1*SQRTS/2.0,partBuf1.mass);
        particle partBuf2 = partBase->getParticle(iq2);
        partBuf2.p4v.SetXYZM(0.0,0.0,x2*SQRTS/-2.0,partBuf2.mass);

        //Generate vector of outgoing fermionic configuration
        vector<particle> confBuf = confBuild.build(iq1,iq2,Q2,rand.GetSeed(),Nf);
        int Nline = confBuf.size()-Nf+2;

        bool color_assigned = false;
        do{
            vector<int> ipcolor;
            vector<int> ipanticolor;
            ipcolor.push_back(-1);
            ipcolor.push_back(-2);
            ipanticolor.push_back(-1);
            ipanticolor.push_back(-2);
            for(int i = 0; i < confBuf.size(); i++){
                if((confBuf[i].pid > 0)&&(confBuf[i].pid < 7)) ipcolor.push_back(i);
                else if ((confBuf[i].pid < 0)&&(confBuf[i].pid > -7)) ipanticolor.push_back(i);
                else if (confBuf[i].pid == 21){
                    ipcolor.push_back(i);
                ipanticolor.push_back(i);
                }
            }
            int iline = 501;
            for(int i = 0; i < Nline; i++){
                int ipc = rand.Integer(ipcolor.size());
                int ipac;
                do{
                  ipac = rand.Integer(ipanticolor.size());
                  if(ipanticolor.size() == 1){
                      if(ipanticolor[ipac]==ipcolor[ipc])  break;
                      else color_assigned = true;
                  }
                }
                while(ipanticolor[ipac]==ipcolor[ipc]);
                if(ipcolor[ipc] < 0){
                    if( ipcolor[ipc] == -1 ) partBuf1.color[1] = iline;
                    else partBuf2.color[1] = iline;
                }
                else{
                    confBuf[ipcolor[ipc]].color[0] = iline;
                }
                if(ipanticolor[ipac] < 0){
                    if( ipanticolor[ipac] == -1 ) partBuf1.color[0] = iline;
                    else partBuf2.color[0] = iline;
                }
                else{
                    confBuf[ipanticolor[ipac]].color[1] = iline;
                }
                ipcolor.erase(ipcolor.begin()+ipc);
                ipanticolor.erase(ipanticolor.begin()+ipac);
                iline++;
            }
        }
        while(!color_assigned);
        vector<particle> inParts;

        inParts.push_back(partBuf1);
        inParts.push_back(partBuf2);

        //Build the mother particle from incoming partons
        particle mother;
        mother.p4v = inParts[0].p4v + inParts[1].p4v;
        mother.mass = mother.p4v.M();

        for(int i = 0; i < confBuf.size(); i++)
        {
            masses.push_back(confBuf[i].mass);
        }
        u1.SetXYZM(0.0,0.0,x1*SQRTS/2.0,inParts[0].mass);
        u2.SetXYZM(0.0,0.0,x2*SQRTS/-2.0,inParts[1].mass);
        mom = u1 + u2;
        momM = mom.M();
        pz = mom.Pz();

        //"Decay" mother to the generated configuration
        weight = 1.0;
        Rambo ramboGeneral(confBuf.size(),mom);
        while(weight > MPW*rand.Uniform())
        {
            //gen.SetDecay(mom, confBuf.size(), &masses[0]);
            //weight = gen.Generate();
            weight = ramboGeneral.Generate();
            if(weight > MPW) cout << "The code needs to recompiled with a higher MPW" << endl;
        }

        //Extract Decay 4-vectors and assign colors to non-spectator quarks
        for(int ii = 0; ii < confBuf.size(); ii++)
        {
            //TLorentzVector prod = *gen.GetDecay(ii);
            TLorentzVector prod = *ramboGeneral.GetDecay(ii);
            confBuf[ii].p4v = prod;
            confBuf[ii].m1 = 1;
            confBuf[ii].m2 = 2;
            int kinI = confBuf[ii].pid + 6;
            if(fabs(confBuf[ii].pid) < 7) 
            {
                pt_hv[kinI].Fill(confBuf[ii].p4v.Pt());
                eta_hv[kinI].Fill(confBuf[ii].p4v.Eta());
                phi_hv[kinI].Fill(confBuf[ii].p4v.Phi());
            }
            outID_h->Fill(confBuf[ii].pid);
        }

        vector<particle> fileParts;
        fileParts.push_back(inParts[0]);
        fileParts.push_back(inParts[1]);

        for(int i = 0; i < confBuf.size(); i++)
        {
            daughtersPt.push_back(confBuf[i].p4v.Pt());
            daughtersEta.push_back(confBuf[i].p4v.Eta());
            daughtersPhi.push_back(confBuf[i].p4v.Phi());
            daughtersE.push_back(confBuf[i].p4v.E());
            daughtersID.push_back(confBuf[i].pid);
            fileParts.push_back(confBuf[i]);
        }




        if(weight > maxwt) 
        {
            maxwt = weight;
            cout << "MaxWt: " << maxwt << endl;
        }

        inQid_h->Fill(iq2,iq1);
        frac2D_h->Fill(x2,x1);
        sphM_h->Fill(momM/1000.0);
        sphPz_h->Fill(pz);
        p1x_h->Fill(x1);
        p1id_h->Fill(iq1);

        lheF.writeEvent(fileParts,sqrt(Q2));
        myT->Fill();
        NF++;
        if(NF%(Nevt/100) == 0) cout << "Produced Event " << NF << "  pdfN : " << pdfN << endl;
        //cout << "Produced Event " << NF << "  pdfN : " << pdfN << endl;
    }

    cout << "Max Weight: " << maxwt << endl;

    lheF.close();
    myF->cd();
    x1_h->Write();
    mcTot_h->Write();
    inQid_h->Write();
    frac2D_h->Write();
    p1x_h->Write();
    p1id_h->Write();
    outID_h->Write();
    sphM_h->Write();
    st_h->Write();
    sphPz_h->Write();
    for(int i = 0; i < pt_hv.size(); i++)
    {
        pt_hv[i].Write();
        eta_hv[i].Write();
        phi_hv[i].Write();
    }
    myT->Write();
    myF->Close();


    return 0;
};






