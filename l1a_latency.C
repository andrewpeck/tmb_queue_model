#include <cmath>
#include "frame_count.h"

// model the TMB as a M/D/1 queue
// Poisson arrivals, constant service times, 1 server

double mean_queue_occupancy (double arrival_rate, double service_rate) {

    // rho is a traffic  parameter
    double rho = arrival_rate/service_rate;
    double n = (rho/(1.0-rho)) * 0.5*(2.0-rho);
    return n;
}


double traffic_intensity (float arrival_rate, float packet_length, float tx_rate) {
    // https://en.wikipedia.org/wiki/Traffic_intensity
    // i.e. rate_in/rate_out
    return (arrival_rate * packet_length / tx_rate);
}

// calculate the blocking probability for a finite queue length with deterministic, fixed readout time (M/D/1 queue)
// cf.
// http://dx.doi.org/10.4218/etrij.14.0113.0812
// theorem 1, formulas 18/19
double blocking_probability (float queue_size, double traffic_intensity) {

    // how to interpret this:
    // the blocking probability is the chance on a per-event basis that we will see that
    // the queue is full and just refuse service to that event (dropping the LCT)
    // and temporarily resolving the full-queue situation.

    // if you want to look at non-blocking occupancies (i.e. if we serve
    // everybody, how full will we be?), then look at the mean queue occupancy
    // or in this function, look at e_k, which is the probability of the queue spot being filled


    int k = queue_size;
    double rho = traffic_intensity;

    //printf ("calculating blocking probability for queue size = %i\n", queue_size);
    //printf ("traffic intensity (write rate/read rate)        = %f\n", traffic_intensity);

    double sigma = 0;
    for (int j=0; j<int(queue_size); j++) {
        sigma += ( pow(-1,double(j)) * pow (rho,double(j)) * pow (double(k-j),double(j)) * TMath::Exp(rho*double(k-j))) / TMath::Factorial(j);
    }

    double e_k = 1-(1-rho)*sigma;

    double pb = (1-rho)*e_k/(1-rho*e_k);

    return pb;
}


//----------------------------------------------------------------------------------------------------------------------
// Calculate the time required to readout the DMB for a single event:
// (add up all the data words, divide by the bandwidth)
//----------------------------------------------------------------------------------------------------------------------
TF1* l1a_latency (int station=0, bool ddr=false, bool gem_en=false, bool deep_buf=false, bool unfurled_triads=false, int n_tbins=12) {


    int n_cfebs = (station==0) ? 7 : 5;

    int dmb_frame_cnt = frame_count(n_cfebs, gem_en, unfurled_triads, n_tbins) / (ddr ? 2 : 1);

    int buffer_size = deep_buf ? 4096 : 2048;
    int minimum_fence = 64;

    //float maximum_occupancy = (float(buffer_size));

    float maximum_occupancy = (float(buffer_size))/(dmb_frame_cnt);
    printf ("    Maximum occupancy    = %f\n" , maximum_occupancy);
    printf ("    Maximum occupancy    = %f = %i\n" , maximum_occupancy, int(maximum_occupancy));

    // we just convert number of bx to kHz
    double dmb_readout_rate = 1000.0/((dmb_frame_cnt)*25.0); // DMB readout rate, Hz
    double service_time     = 1./dmb_readout_rate; // DMB readout time, seconds
    double service_time_bx  = 500+dmb_frame_cnt ;// DMB readout time, bx

    printf ("    Word Count           = %i\n" , dmb_frame_cnt);
    printf ("    DMB Readout Rate     = %f kHz (l1a/sec)\n" , dmb_readout_rate*1000.);


    //TCanvas *c1 = new TCanvas("c1");

    #include "lumi_scaling.h"

    TF1 *ovf = new TF1 (
            "ovf",
            "blocking_probability([0], x*[1])",
            0,
            50
    );

    ovf->SetParameter(0,maximum_occupancy);

    // lct*l1a rate  [Hz]   = lumi_rate [10^34] * 10000 * lumi_scaler
    //double intensity_coefficient = (lumi_rate_me11 * 10000.0 * dmb_frame_cnt) / (pow(10,9) * 1.0/25);
    double intensity_coefficient = (lumi_rate_arr[station] * 10000.0 * dmb_frame_cnt) / (40000000.);

    ovf->SetParameter(1,  intensity_coefficient); // service time has units seconds/packet


    ovf->SetNpx(100);
    ovf->GetHistogram()->GetXaxis()->SetTitle("luminosity (10^34 cm-2 s-1)");
    //ovf->Draw();
    //c1->SetLogy();
    ovf->SetTitle("overflow probability");


    //////////////////////////////////////////////////////////////

//    bool draw_c2 = 0;
//    if (draw_c2) {
//    TCanvas *c2 = new TCanvas("c2");
//
//    TF1 *occ = new TF1 (
//            "occ",
//            "mean_queue_occupancy(x, [0])",
//            0,
//            dmb_readout_rate*9/10
//    );
//
//    occ->SetParameter(0,dmb_readout_rate*9/10);
//    occ->SetNpx(10000);
//    occ->Draw();
//    c2->SetLogy();


//    occ->GetHistogram()->GetXaxis()->SetTitle("l1a*lct rate (kHz)");
//    occ->GetHistogram()->GetYaxis()->SetTitle("events");
//    occ->SetTitle("mean l1as in readout queue");
//    }

    return ovf;
}
