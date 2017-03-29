#include "l1a_latency.C"

void gen_plots () {

    TCanvas *c1 = new TCanvas("c1");
    c1->SetWindowSize(1800, 1280);
    c1->Divide(2,2);

    TFile* hfile = new TFile("tmb_model.root","READ","TMB Queue Model");

    TF1* me11_theory = l1a_latency();

    TH2F* th2f_lostevents_me11 = (TH2F*) hfile->Get("Lost Events ME11");
    TH2F* th2f_lostevents_me21 = (TH2F*) hfile -> Get ("Lost Events ME21");
    TH2F* th2f_lostevents_me31 = (TH2F*) hfile -> Get ("Lost Events ME31");
    TH2F* th2f_lostevents_me41 = (TH2F*) hfile -> Get ("Lost Events ME41");

    TH2F* th2f_lostevents_ddr_me11 = (TH2F*) hfile -> Get ("Lost Events DDR ME11");
    TH2F* th2f_lostevents_ddr_me21 = (TH2F*) hfile -> Get ("Lost Events DDR ME21");
    TH2F* th2f_lostevents_ddr_me31 = (TH2F*) hfile -> Get ("Lost Events DDR ME31");
    TH2F* th2f_lostevents_ddr_me41 = (TH2F*) hfile -> Get ("Lost events DDR ME41");

    TH2F* th2f_lostevents_gem_me11 = (TH2F*) hfile -> Get("Lost Events GEM ME11");

    TH2F* th2f_lostevents_deep_me11   = (TH2F*) hfile -> Get ("Lost Events Deep Buffer ME11");

    TH2F* th2f_lostevents_low_l1_me11 = (TH2F*) hfile -> Get ("Lost Events LOWL1 ME11");

    TH2F* th2f_occupancy = (TH2F*) hfile->Get("Mean Queue Occupancy");

    TH2F* lostevents_arr [4] = {th2f_lostevents_me11, th2f_lostevents_me21, th2f_lostevents_me31,  th2f_lostevents_me41};
    TH2F* lostevents_ddr_arr [4] = {th2f_lostevents_ddr_me11, th2f_lostevents_ddr_me21, th2f_lostevents_ddr_me31,  th2f_lostevents_ddr_me41};

    //------------------------------------------------------------------------------------------------------------------
    // Profiles
    //------------------------------------------------------------------------------------------------------------------

    //ProfileX (const char *name="_pfx", Int_t firstybin=1, Int_t lastybin=-1, Option_t *option="") const
    TProfile *me11 = th2f_lostevents_me11->ProfileX("px11",1,-1,"o");
    TProfile *me21 = th2f_lostevents_me21->ProfileX("px21",1,-1,"o");
    TProfile *me31 = th2f_lostevents_me31->ProfileX("px31",1,-1,"o");
    TProfile *me41 = th2f_lostevents_me41->ProfileX("px41",1,-1,"o");

    TProfile *me11_ddr   = th2f_lostevents_ddr_me11->ProfileX("px11_ddr",1,-1,"o");
    TProfile *me11_deep  = th2f_lostevents_deep_me11->ProfileX("px11_deep",1,-1,"o");
    TProfile *me11_lowl1 = th2f_lostevents_low_l1_me11->ProfileX("px11_lowl1",1,-1,"o");
    TProfile *me11_gem   = th2f_lostevents_gem_me11->ProfileX("px11_gem",1,-1,"o");

    me11->SetMaximum(1.0);
    me11->SetMinimum(0.0);

    me11_ddr->Rebin(4);
    me11_deep->Rebin(4);
    me11_lowl1->Rebin(4);
    me11_gem->Rebin(4);

    me11->Rebin(4);
    me21->Rebin(4);
    me31->Rebin(4);
    me41->Rebin(4);

    //------------------------------------------------------------------------------------------------------------------
    //
    //------------------------------------------------------------------------------------------------------------------

    (lostevents_arr[0])->SetMarkerColor(kBlue-9  );
    (lostevents_arr[1])->SetMarkerColor(kRed-4   );
    (lostevents_arr[2])->SetMarkerColor(kOrange+1);
    (lostevents_arr[3])->SetMarkerColor(kGreen-6 );

    (lostevents_arr[0])->SetMarkerStyle(8);
    (lostevents_arr[1])->SetMarkerStyle(8);
    (lostevents_arr[2])->SetMarkerStyle(8);
    (lostevents_arr[3])->SetMarkerStyle(8);

    me11->SetMarkerColor(kBlue-9  );
    me21->SetMarkerColor(kRed-3   );
    me31->SetMarkerColor(kOrange+1);
    me41->SetMarkerColor(kGreen-6 );

    me11->SetMarkerStyle(8);
    me21->SetMarkerStyle(8);
    me31->SetMarkerStyle(8);
    me41->SetMarkerStyle(8);

    (lostevents_ddr_arr[0])->SetMarkerColor(1);
    (lostevents_ddr_arr[1])->SetMarkerColor(1);
    (lostevents_ddr_arr[2])->SetMarkerColor(1);
    (lostevents_ddr_arr[3])->SetMarkerColor(1);

    (lostevents_ddr_arr[0])->SetMarkerStyle(3);
    (lostevents_ddr_arr[1])->SetMarkerStyle(3);
    (lostevents_ddr_arr[2])->SetMarkerStyle(3);
    (lostevents_ddr_arr[3])->SetMarkerStyle(3);


    //------------------------------------------------------------------------------------------------------------------
    //
    //------------------------------------------------------------------------------------------------------------------

    c1->cd(1);

    // th2f_occupancy->Draw();
    // th2f_occupancy->SetMarkerStyle(3);

    me11        -> SetStats(0);
    me11        -> Draw();
    me11_lowl1  -> Draw("SAME");

    me11_lowl1->SetMarkerStyle(3);
    me11_lowl1->SetMarkerColor(kBlack);

    TLegend* leg1 = new TLegend(0.1,0.65,0.48,0.85);
    leg1->AddEntry(me11,      "ME1/1 500 bx Latency","pe");
    leg1->AddEntry(me11_lowl1,"ME1/1 128 bx Latency","pe");
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->Draw();


    //------------------------------------------------------------------------------------------------------------------
    // ME1/1 Mitigation Comparison
    //------------------------------------------------------------------------------------------------------------------

    c1->cd(3);

    me11        -> SetStats(0);
    me11        -> Draw();
    me11_ddr    -> Draw("SAME");
    me11_deep   -> Draw("SAME");
    me11_gem    -> Draw("SAME");

    me11_deep->SetMarkerStyle(3);
    me11_deep->SetMarkerColor(kRed-4);

    me11_ddr->SetMarkerStyle(3);
    me11_ddr->SetMarkerColor(kGreen-3);

    me11_gem->SetMarkerStyle(3);
    me11_gem->SetMarkerColor(kMagenta-6);

    TLegend* leg3 = new TLegend(0.1,0.65,0.48,0.85);
    leg3->AddEntry(me11,"Baseline","pe");
    leg3->AddEntry(me11_deep,"Increased buffer depth","pe");
    leg3->AddEntry(me11_gem,"GEM Readout","pe");
    leg3->AddEntry(me11_ddr,"DDR Readout","pe");
    leg3->SetFillStyle(0);
    leg3->SetBorderSize(0);
    leg3->Draw();

    //------------------------------------------------------------------------------------------------------------------
    // Station Comparison
    //------------------------------------------------------------------------------------------------------------------

    c1->cd(2);

    me11        -> Draw();
    me11_theory -> Draw("SAME");
    //me11_theory -> Draw();

    TLegend* leg2 = new TLegend(0.1,0.65,0.48,0.85);
    leg2->AddEntry(me11,"ME1/1 Monte Carlo","pe");
    leg2->AddEntry(me11_theory,"ME1/1 M/D/1 Queue Model","l");
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->Draw();

    //------------------------------------------------------------------------------------------------------------------
    //
    //------------------------------------------------------------------------------------------------------------------

    c1->cd(4);

    me11->Draw();
    me21->Draw("SAME");
    me31->Draw("SAME");
    me41->Draw("SAME");

    TLegend* leg4 = new TLegend(0.1,0.65,0.48,0.85);
    leg4->AddEntry(me11,"ME1/1","pe");
    leg4->AddEntry(me21,"ME2/1","pe");
    leg4->AddEntry(me31,"ME3/1","pe");
    leg4->AddEntry(me41,"ME4/1","pe");
    leg4->SetFillStyle(0);
    leg4->SetBorderSize(0);
    leg4->Draw();


}
