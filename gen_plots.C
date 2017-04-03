#include "l1a_latency.C"

#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TProfile.h>

void gen_plots () {

    TCanvas *c1 = new TCanvas("c1");
    c1->SetWindowSize(2600, 1280);
    c1->Divide(4,2);

    TFile* hfile = new TFile("tmb_model.root","READ","TMB Queue Model");

    TF1* me11_theory = l1a_latency(0);
    TF1* me21_theory = l1a_latency(1);
    TF1* me31_theory = l1a_latency(2);
    TF1* me41_theory = l1a_latency(3);

    // TF1* l1a_latency (int station=0, bool ddr=false, bool gem_en=false, bool deep_buf=false, bool unfurled_triads=false, int n_tbins=12) {

    TF1* me11_theory_ddr          = l1a_latency(0, true);
    TF1* me11_theory_gem          = l1a_latency(0, false, true);
    TF1* me11_theory_gem_ddr      = l1a_latency(0, true,  true);
    TF1* me11_theory_deep         = l1a_latency(0, false, false, true);
    TF1* me11_theory_unfurled_ddr = l1a_latency(0, true,  false, false, true);
    TF1* me11_theory_14tbins      = l1a_latency(0, false, false, false, false, 14);
    TF1* me11_theory_14tbins_ddr  = l1a_latency(0, true,  false, false, false, 14);

    me11_theory -> SetLineColor(kBlack);
    me21_theory -> SetLineColor(kBlack);
    me31_theory -> SetLineColor(kBlack);
    me41_theory -> SetLineColor(kBlack);

    me11_theory -> SetLineStyle(2);
    me21_theory -> SetLineStyle(2);
    me31_theory -> SetLineStyle(2);
    me41_theory -> SetLineStyle(2);

    me11_theory_ddr          -> SetLineStyle(2);
    me11_theory_gem          -> SetLineStyle(2);
    me11_theory_gem_ddr      -> SetLineStyle(2);
    me11_theory_deep         -> SetLineStyle(2);
    me11_theory_unfurled_ddr -> SetLineStyle(2);

    me11_theory_ddr          -> SetLineColor (kBlack);
    me11_theory_gem          -> SetLineColor (kBlack);
    me11_theory_gem_ddr      -> SetLineColor (kBlack);
    me11_theory_deep         -> SetLineColor(kBlack);
    me11_theory_unfurled_ddr -> SetLineColor(kBlack);


    TH2F* th2f_lostevents_me11 = (TH2F*) hfile -> Get ("h2_loss_me11");
    TH2F* th2f_lostevents_me21 = (TH2F*) hfile -> Get ("h2_loss_me21");
    TH2F* th2f_lostevents_me31 = (TH2F*) hfile -> Get ("h2_loss_me31");
    TH2F* th2f_lostevents_me41 = (TH2F*) hfile -> Get ("h2_loss_me41");

    TH2F* th2f_lostevents_ddr_me11 = (TH2F*) hfile -> Get ("h2_loss_me11_ddr");
    TH2F* th2f_lostevents_ddr_me21 = (TH2F*) hfile -> Get ("h2_loss_me21_ddr");
    TH2F* th2f_lostevents_ddr_me31 = (TH2F*) hfile -> Get ("h2_loss_me31_ddr");
    TH2F* th2f_lostevents_ddr_me41 = (TH2F*) hfile -> Get ("h2_loss_me41_ddr");

    TH2F* th2f_lostevents_unfurled_ddr_me11 = (TH2F*) hfile -> Get ("h2_loss_me11_unfurled_ddr");

    TH2F* th2f_lostevents_gem_me11 = (TH2F*) hfile -> Get("h2_loss_me11_gem");
    TH2F* th2f_lostevents_gem_ddr_me11 = (TH2F*) hfile -> Get("h2_loss_me11_gem_ddr");

    TH2F* th2f_lostevents_deep_me11   = (TH2F*) hfile -> Get ("h2_loss_me11_deep");

    TH2F* th2f_lostevents_low_l1_me11 = (TH2F*) hfile -> Get ("h2_loss_me11_lowl1");


    TH2F* lostevents_arr     [4] = {th2f_lostevents_me11     , th2f_lostevents_me21     , th2f_lostevents_me31     , th2f_lostevents_me41};
    TH2F* lostevents_ddr_arr [4] = {th2f_lostevents_ddr_me11 , th2f_lostevents_ddr_me21 , th2f_lostevents_ddr_me31 , th2f_lostevents_ddr_me41};

    //------------------------------------------------------------------------------------------------------------------
    // Profiles
    //------------------------------------------------------------------------------------------------------------------

    //ProfileX (const char *name="_pfx", Int_t firstybin=1, Int_t lastybin=-1, Option_t *option="") const

    TProfile *me11 = th2f_lostevents_me11 -> ProfileX("me11" , 1 , -1 , "o");
    TProfile *me21 = th2f_lostevents_me21 -> ProfileX("me21" , 1 , -1 , "o");
    TProfile *me31 = th2f_lostevents_me31 -> ProfileX("me31" , 1 , -1 , "o");
    TProfile *me41 = th2f_lostevents_me41 -> ProfileX("me41" , 1 , -1 , "o");

    me11 -> GetYaxis() -> SetTitle("loss rate");
    me21 -> GetYaxis() -> SetTitle("loss rate");
    me31 -> GetYaxis() -> SetTitle("loss rate");
    me41 -> GetYaxis() -> SetTitle("loss rate");

    TProfile *me11_ddr          = th2f_lostevents_ddr_me11          -> ProfileX("me11_ddr"          , 1 , -1 , "o");
    TProfile *me11_deep         = th2f_lostevents_deep_me11         -> ProfileX("me11_deep"         , 1 , -1 , "o");
    TProfile *me11_lowl1        = th2f_lostevents_low_l1_me11       -> ProfileX("me11_lowl1"        , 1 , -1 , "o");
    TProfile *me11_gem          = th2f_lostevents_gem_me11          -> ProfileX("me11_gem"          , 1 , -1 , "o");
    TProfile *me11_gem_ddr      = th2f_lostevents_gem_ddr_me11      -> ProfileX("me11_gem_ddr"      , 1 , -1 , "o");
    TProfile *me11_unfurled_ddr = th2f_lostevents_unfurled_ddr_me11 -> ProfileX("me11_unfurled_ddr" , 1 , -2 , "0");

    me11->SetMaximum(1.0);
    me11->SetMinimum(0.0);

    me11_ddr          -> Rebin(4);
    me11_deep         -> Rebin(4);
    me11_lowl1        -> Rebin(4);
    me11_gem          -> Rebin(4);
    me11_gem_ddr      -> Rebin(4);
    me11_unfurled_ddr -> Rebin(4);

    me11->Rebin(50);
    me21->Rebin(50);
    me31->Rebin(50);
    me41->Rebin(50);

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

    gStyle->SetOptTitle(0);

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

    //-Baseline---------------------------------------------------------------------------------------------------------
    me11        -> Draw();

    //-Alternates-------------------------------------------------------------------------------------------------------

    me11_ddr          -> Draw("SAME");
    me11_deep         -> Draw("SAME");

    me11_deep->SetMarkerStyle(3);
    me11_deep->SetMarkerColor(kRed-4);

    me11_ddr->SetMarkerStyle(3);
    me11_ddr->SetMarkerColor(kGreen-3);

    //-Theory-----------------------------------------------------------------------------------------------------------

    me11_theory              -> Draw("SAME");
    me11_theory_ddr          -> Draw("SAME");
    me11_theory_deep         -> Draw("SAME");


    //-Legend-----------------------------------------------------------------------------------------------------------

    TLegend* leg3 = new TLegend(0.1,0.65,0.48,0.85);

    leg3->AddEntry(me11,"Baseline","pe");
    leg3->AddEntry(me11_deep,"Buffer depth x2","pe");
    leg3->AddEntry(me11_ddr,"DDR Readout","pe");
    leg3->AddEntry(me11_theory,"M/D/1 Models","l");

    leg3->SetFillStyle(0);
    leg3->SetBorderSize(0);
    leg3->Draw();

    //------------------------------------------------------------------------------------------------------------------
    // ME1/1 GEM Comparison
    //------------------------------------------------------------------------------------------------------------------

    c1->cd(7);

    me11        -> SetStats(0);

    //-Baseline---------------------------------------------------------------------------------------------------------
    me11        -> Draw();

    //-Alternates-------------------------------------------------------------------------------------------------------

    me11_gem          -> Draw("SAME");
    me11_gem_ddr      -> Draw("SAME");

    me11_gem->SetMarkerStyle(3);
    me11_gem->SetMarkerColor(kMagenta-6);

    //-Theory-----------------------------------------------------------------------------------------------------------

    me11_theory              -> Draw("SAME");
    me11_theory_gem          -> Draw("SAME");
    me11_theory_ddr          -> Draw("SAME");
    me11_theory_gem_ddr      -> Draw("SAME");

    //-Legend-----------------------------------------------------------------------------------------------------------

    TLegend* leg7 = new TLegend(0.1,0.65,0.48,0.85);

    leg7->AddEntry(me11,    "Baseline, l1=500bx","pe");
    leg7->AddEntry(me11_gem,"Baseline with GEM, l1=500bx","pe");
    leg7->AddEntry(me11_ddr,"DDR Readout, l1=500bx","pe");
    leg7->AddEntry(me11_gem_ddr,"DDR Readout with GEM, l1=500bx","pe");

    leg7->AddEntry(me11_theory,"M/D/1 Models","l");

    leg7->SetFillStyle(0);
    leg7->SetBorderSize(0);
    leg7->Draw();

    //------------------------------------------------------------------------------------------------------------------
    // ME1/1 Unfurled Triads
    //------------------------------------------------------------------------------------------------------------------

    c1->cd(8);

    me11        -> SetStats(0);

    //-Baseline---------------------------------------------------------------------------------------------------------
    me11        -> Draw();

    //-Alternates-------------------------------------------------------------------------------------------------------

    me11_unfurled_ddr -> Draw("SAME");

    me11_unfurled_ddr->SetMarkerStyle(3);
    me11_unfurled_ddr->SetMarkerColor(kGreen+4);

    //-Theory-----------------------------------------------------------------------------------------------------------

    me11_theory              -> Draw("SAME");
    me11_theory_unfurled_ddr -> Draw("SAME");


    //-Legend-----------------------------------------------------------------------------------------------------------

    TLegend* leg8 = new TLegend(0.1,0.65,0.48,0.85);

    leg8->AddEntry(me11,"Baseline","pe");
    leg8->AddEntry(me11_unfurled_ddr,"Unfurled Triads w/ DDR","pe");
    leg8->AddEntry(me11_theory,"M/D/1 Models","l");

    leg8->SetFillStyle(0);
    leg8->SetBorderSize(0);
    leg8->Draw();

    Float_t ymax = 1;
    //TLine *line = new TLine(0,ymax/2,50,ymax/2);
    TLine *line = new TLine(7.5,0,7.5,ymax/2);
    line->SetLineColor(kGray);
    line->Draw();


    //------------------------------------------------------------------------------------------------------------------
    // Monte Carlo vs. Data
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
    // Station Comparison
    //------------------------------------------------------------------------------------------------------------------

    c1->cd(4);

    me11->Draw();
    me21->Draw("SAME");
    me31->Draw("SAME");
    me41->Draw("SAME");

    me11_theory->Draw("SAME");
    me21_theory->Draw("SAME");
    me31_theory->Draw("SAME");
    me41_theory->Draw("SAME");

    TLegend* leg4 = new TLegend(0.1,0.65,0.48,0.85);
    leg4->AddEntry(me11,"ME1/1","pe");
    leg4->AddEntry(me21,"ME2/1","pe");
    leg4->AddEntry(me31,"ME3/1","pe");
    leg4->AddEntry(me41,"ME4/1","pe");
    leg4->SetFillStyle(0);
    leg4->SetBorderSize(0);
    leg4->Draw();

    //------------------------------------------------------------------------------------------------------------------
    //
    //------------------------------------------------------------------------------------------------------------------

    c1->cd(5);


    me11_theory->Draw();
    me11_theory_14tbins      ->Draw("same");
    me11_theory_14tbins_ddr  ->Draw("same");
    me11_theory_ddr          ->Draw("same");

    TLegend* leg5 = new TLegend(0.15,0.65,0.48,0.85);

    leg5->AddEntry(me11_theory             , "ME1/1 12 tbins @ 40MHz" , "l");
    leg5->AddEntry(me11_theory_14tbins     , "ME1/1 14 tbins @ 40MHz" , "l");

    leg5->AddEntry(me11_theory_ddr         , "ME1/1 12 tbins @ 80MHz" , "l");
    leg5->AddEntry(me11_theory_14tbins_ddr , "ME1/1 14 tbins @ 80MHz" , "l");

    leg5->SetFillStyle(0);
    leg5->SetBorderSize(0);

    leg5->Draw();

    //------------------------------------------------------------------------------------------------------------------
    //
    //------------------------------------------------------------------------------------------------------------------

    c1->cd(6);

    TH2F* th2f_occupancy   = (TH2F*) hfile->Get("h2_queue_occupancy");
    TH2F* th2f_event_sep   = (TH2F*) hfile->Get("h2_event_sep");
    TH2F* th2f_pretrig_sep = (TH2F*) hfile->Get("h2_pretrig_sep");

    //TProfile *event_sep = th2f_event_sep -> ProfileX("event sepaparation me1/1" , 1 , -1 , "o");
    TProfile *pretrig_sep = th2f_pretrig_sep -> ProfileX("pretrig separation me1/1" , 1 , -1 , "o");

    pretrig_sep->SetMarkerColor(kRed);
    pretrig_sep->SetMarkerStyle(3);

    pretrig_sep->Draw();

}
