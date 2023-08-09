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
#include "../include/linearInterpolation.h"
#define ARGS 7
#define MPW 5.0e-2

using namespace LHAPDF;
using namespace std;
using namespace constants;

int main(int argc, char* argv[])
{

    //read arguements and do conversions and checks on them
    if(argc != ARGS+1 && argc != ARGS+2)
    {
        cout << "InstantonGEN sqrtS minMass maxweight Nevents Nf Filename [isWeighted (default 0)]" << endl;
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
    bool isweight = false;
    if(argc == ARGS+2) isweight = atoi(argv[7]);

    //Init structures needed during generatation
    double maxwt = 0;
    double maxMCtot = 0.0;

    int NF = 0;
    TRandom3 rand;
    int rndmSeed = atoi(argv[7]);
    rand.SetSeed(rndmSeed);
    particleBase *partBase = new particleBase();
    configBuilder confBuild;
    vector<double> daughtersPt;
    vector<double> daughtersEta;
    vector<double> daughtersPhi;
    vector<double> daughtersE;
    vector<double> daughtersID;
    TLorentzVector mom(0.0,0.0,0.0,0.0);
    TLorentzVector u1(0.0,0.0,0.0,0.0);
    TLorentzVector u2(0.0,0.0,0.0,0.0);

    double weight = 1.0;

    //TGenPhaseSpace gen;
    double x1 = 0.0;
    double x2 = 0.0;
    int iq1 = 0;
    int iq2 = 0;
    double Q2 = 0.0;
    double momM = 0.0;
    double pz = 0.0;
    double rapidity = 0.0;
    double st = 0.0;
    int nd = 0;
    //Init output files
    TFile *myF = new TFile((ofName+".root").c_str(),"RECREATE","Holds daughters from instanton decay");
    LHEWriter lheF(ofName, SQRTS, isweight);

    //Prepare for the importance sampling according to the parton-level XS.
    //Get the integral of parton level XS from the list of differential XS.
    double XSIntegral[LEN];
    getIntegral(&XSIntegral[0], ENERGYPARTON, XS, LEN);
    //Get the cumulative distribution function and the probability distribution (proportional to parton level XS) of instanton mass
    double CDF[LEN];
    getCDF(&CDF[0], ENERGYPARTON, XS, LEN);
    double XSPDF[LEN];
    getPDF(&XSPDF[0], ENERGYPARTON, XS, LEN);
    //Get the value of CDF at the lower boundary of mass and the integrated XS above this threshold.
    double CDFThr = getInterpoCDF(thr, ENERGYPARTON, CDF, XSPDF, LEN);
    double XSIntegralThr = XSIntegral[LEN-1] * (1 - CDFThr);


    string massList = "Instanton mass list: ";
    for(int i = 0; i < LEN; i++) massList += std::to_string(ENERGYPARTON[i])+"\t";

    string XSList = "Instanton XS list: ";
    for(int i = 0; i < LEN; i++) XSList += std::to_string(XS[i])+"\t";

    string CDFList = "Instanton CDF list: ";
    for(int i = 0; i < LEN; i++) CDFList += std::to_string(CDF[i])+"\t";

    string PDFList = "Instanton PDF list: ";
    for(int i = 0; i < LEN; i++) PDFList += std::to_string(XSPDF[i])+"\t";

    cout<<massList<<endl;
    cout<<XSList<<endl;
    cout<<CDFList<<endl;
    cout<<PDFList<<endl;

    cout<<"CDFThr: "<<CDFThr<<endl;
    //Init parton distribution functions
    const PDF* LHApdf = mkPDF("CT10",0);
    //cout << LHApdf->xfxQ2(2, 0.5, SQRTS*SQRTS) << endl;

    //Initialize histograms for debugging
    TH1D *x1_h = new TH1D("x1_h","x1 inclusive",1000,0.0,1.0);
    TH1D *mcTot_h = new TH1D("mcTot_h","Monte Carlo Probabilities",100,0.0,MCW);
    TH1D *instantonM_h = new TH1D("instantonM_h","Instanton Transition Energy;Invariant Mass [GeV];Events / 3 GeV",100,0.0,300.0);
    TH1D *st_h = new TH1D("st_h","S_{T};S_{T} [GeV];Events / 3 GeV", 100.0 , 0.0, 300.0);
    TH1D *instantonPz_h = new TH1D("instantonPz_h","Instanton p_{z};p_{z} [GeV];Events / 100 GeV",80,-4000.0,4000.0);
    TH1D *nd_h = new TH1D("nd_h","Daughter multiplicity;Daughter multiplicity;Events", 100.0 , 0.0, 100.0);
    TH1D *outID_h = new TH1D("outID_h","Outgoing PDG IDs;PDG ID;Entries",33,-16.5,16.5);
    TH1D *p1x_h = new TH1D("p1x_h","Parton 1 Momentum Fraction;x_{1};Events / 0.01",60,0.4,1.0);
    TH1D *p1id_h = new TH1D("p1id_h","Parton 1 PDG ID;PDG ID;Events",13,-6.5,6.5);
    TH1D *rapidity_h = new TH1D("instantonRapidity_h","Instanton rapidity",100,-8,8);
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
    myT->Branch("rapidity",&rapidity);
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
        double valXSCDF;
        double x1x2;
        double valInstMass;
        double valx1CDF;
        while(!mcPass)
        {
            //Choose x1 and x2 in proper range
            //cout<<"CDFThr: "<<CDFThr<<endl;
            //Importance sampling of the instanton mass according to parton level XS.
            valXSCDF = rand.Uniform(CDFThr,1);
            valInstMass = invertCDF(valXSCDF, ENERGYPARTON, CDF, XSPDF, LEN);
            Q2 = valInstMass * valInstMass;
            x1x2 = Q2 / (SQRTS * SQRTS);
            //cout<<"x1x2: "<<x1x2<<endl;
            //After getting x1x2 from the instanton mass
            //Importance sampling of x1 according to f(x1|x1x2) = - 1/(x1 * Ln(x1x2)), when x1~U(0,1), x2~U(0,1)
            valx1CDF = rand.Uniform(0,1);
            x1 = TMath::Power(x1x2, 1 - valx1CDF);
            //cout<<"x1: "<<x1<<endl;
            x2 = x1x2 / x1; 
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
                    //cout<<"x1: "<<x1<<", x2: "<<x2<<", x1p: "<<x1p<<", x2p: "<<x2p<<", x1p*x2p: "<<x1p*x2p<<endl;
                    mcTot += x1p*x2p;
                    if(mcTot > mcP) {mcPass = true; break;}
                }
                if(mcPass) break;
            }

            //Weight the sample by pdf values of gluons at x1, x2
            if(isweight){
                weight = mcTot * XSIntegralThr; // pdf of the partons * integrated XS(as constant)
                mcPass = true;
            }
            //For unweighted genration, mcPass is turned true with probability x1p*x2p/MCW
            //-> hit-or-miss sampling according to the parton distribution function.
            pdfN++;
            x1_h->Fill(x1);
            mcTot_h->Fill(mcTot);
            if(maxMCtot < mcTot) {maxMCtot = mcTot; cout << "Max MC Total: " << maxMCtot << endl;}
        }


        //Build incoming particles, instanton, and prepare to decay
        particle partBuf1 = partBase->getParticle(iq1);
        if(iq1 != partBuf1.pid) cout << "iq1 = " << iq1 << " != partBuf.pid = " << partBuf1.pid << endl;
        partBuf1.p4v.SetXYZM(0.0,0.0,x1*SQRTS/2.0,partBuf1.mass);
        particle partBuf2 = partBase->getParticle(iq2);
        partBuf2.p4v.SetXYZM(0.0,0.0,x2*SQRTS/-2.0,partBuf2.mass);

        //Generate vector of outgoing fermionic configuration
        vector<particle> confBuf = confBuild.build(iq1,iq2,Q2,rand.GetSeed(),Nf);

        //Assign the color lines to the daughters
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
        vector<double> masses;
        for(int i = 0; i < confBuf.size(); i++)
        {
            masses.push_back(confBuf[i].mass);
        }
        u1.SetXYZM(0.0,0.0,x1*SQRTS/2.0,inParts[0].mass);
        u2.SetXYZM(0.0,0.0,x2*SQRTS/-2.0,inParts[1].mass);
        mom = u1 + u2;
        momM = mom.M();
        pz = mom.Pz();
        rapidity = mom.Rapidity();
        nd = confBuf.size();

        //"Decay" mother to the generated configuration
        Rambo ramboGeneral(confBuf.size(),mom,masses,rand.GetSeed());
        ramboGeneral.Generate();

        //Extract Decay 4-vectors
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
        st = 0;

        for(int i = 0; i < confBuf.size(); i++)
        {
            daughtersPt.push_back(confBuf[i].p4v.Pt());
            daughtersEta.push_back(confBuf[i].p4v.Eta());
            daughtersPhi.push_back(confBuf[i].p4v.Phi());
            daughtersE.push_back(confBuf[i].p4v.E());
            daughtersID.push_back(confBuf[i].pid);
            fileParts.push_back(confBuf[i]);
            st += confBuf[i].p4v.Pt();
        }




        if(weight > maxwt) 
        {
            maxwt = weight;
            cout << "MaxWt: " << maxwt << endl;
        }

        inQid_h->Fill(iq2,iq1);
        frac2D_h->Fill(x2,x1);
        instantonM_h->Fill(momM);
        instantonPz_h->Fill(pz);
        st_h->Fill(st);
        nd_h->Fill(confBuf.size());
        p1x_h->Fill(x1);
        p1id_h->Fill(iq1);
        rapidity_h->Fill(rapidity);
        lheF.writeEvent(fileParts,sqrt(Q2), weight);
        myT->Fill();
        NF++;
        if(NF%(Nevt/10) == 0) cout << "Produced Event " << NF << "  pdfN : " << pdfN << endl;
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
    instantonM_h->Write();
    st_h->Write();
    nd_h->Write();
    instantonPz_h->Write();
    rapidity_h->Write();
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






