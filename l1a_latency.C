#include <cmath>

// model the TMB as a M/D/1 queue
// Poisson arrivals, constant service times, 1 server

double mean_queue_occupancy (double arrival_rate, double service_rate) {

    // rho is a traffic  parameter
    double rho = arrival_rate/service_rate;
    double n = (rho/(1.0-rho)) * 0.5*(2.0-rho);
    return n;
}

// calculate the blocking probability for a finite queue length with deterministic, fixed readout time (M/D/1 queue)
// cf. http://dx.doi.org/10.4218/etrij.14.0113.0812
// theorem 1, formulas 18/19
double blocking_probability (int queue_size, double traffic_intensity) {

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
    for (int j=0; j<queue_size; j++) {
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
void l1a_latency () {

    int n_headers    = 42;
    int n_cfebs      = 7;
    int n_rpcs       = 2;
    int n_tbins      = 6;
    int n_scope_chan = 128;

    bool rpc_en               = true;
    bool scope_en             = false;
    bool miniscope_en         = false;
    bool blocked_cfeb_readout = false;

    // header
    int wc_eob        = 1;

    // cfeb
    int wc_cfeb       = n_cfebs * 6 * n_tbins;

    // rpc
    int wc_b04        = rpc_en ? 1 : 0;
    int wc_rpcs       = rpc_en ? 2*n_rpcs*n_tbins : 0;
    int wc_e04        = rpc_en ? 1 : 0;

    // scope;
    int wc_b05        = scope_en ? 1                   : 0;
    int wc_scope      = scope_en ? n_scope_chan/16*256 : 0;
    int wc_e05        = scope_en ? 1                   : 0;

    // miniscope;
    int wc_b07        = miniscope_en ? 1  : 0;
    int wc_miniscope  = miniscope_en ? 22 : 0;

    // blocked cfeb readout;
    int wc_bcb        = blocked_cfeb_readout ? 1  : 0;
    int wc_b_cfeb     = blocked_cfeb_readout ? 22 : 0;
    int wc_ecb        = blocked_cfeb_readout ? 1  : 0;

    int wc_eoc        = 1;
    int wc_multiple   = 0; // we account for this later;
    int wc_eof        = 1;
    int wc_crc        = 2;
    int wc_wc         = 1;

    int dmb_frame_cnt;

    dmb_frame_cnt = n_headers + wc_eob + wc_cfeb + wc_b04 + wc_rpcs + wc_e04 + wc_b05 + wc_scope + wc_e05 + wc_b07 + wc_miniscope + wc_bcb + wc_b_cfeb + wc_ecb + wc_eoc + wc_multiple + wc_eof + wc_crc + wc_wc;
    dmb_frame_cnt = (dmb_frame_cnt+3) & ~(0x03); //  round up to nearest multiple of 4


    int buffer_size = 2048;
    int minimum_fence = 64;

    int maximum_occupancy = (buffer_size/dmb_frame_cnt);

    // we just convert number of bx to kHz
    double dmb_readout_rate = 1000000.0/(dmb_frame_cnt*25); // DMB readout time, kHz
    double service_time     = 1./dmb_readout_rate;



    printf ("Configured with:\n");
    printf ("    n_cfebs              = %i\n" , n_cfebs);
    printf ("    n_tbins              = %i\n" , n_tbins);
    printf ("    rpc_readout          = %i\n" , rpc_en);
    printf ("    scope_enabled        = %i\n" , scope_en);
    printf ("    miniscope_enabled    = %i\n" , miniscope_en);
    printf ("    blocked_cfeb_readout = %i\n" , blocked_cfeb_readout);
    printf ("\n");
    printf ("    Word Count           = %i\n" , dmb_frame_cnt);
    printf ("    DMB Readout Rate     = %f kHz (l1a/sec)\n" , dmb_readout_rate);


    TCanvas *c1 = new TCanvas("c1");

    TF1 *ovf = new TF1 (
            "ovf",
            "blocking_probability([0], x*[1])",
            0,
            dmb_readout_rate*2
    );

    ovf->SetParameter(0,maximum_occupancy);
    ovf->SetParameter(1,service_time);
    ovf->SetNpx(10000);
    ovf->GetHistogram()->GetXaxis()->SetTitle("l1a*lct rate (kHz)");
    ovf->Draw();
    c1->SetLogy();
    ovf->SetTitle("overflow probability");


    //////////////////////////////////////////////////////////////

    TCanvas *c2 = new TCanvas("c2");

    TF1 *occ = new TF1 (
            "occ",
            "mean_queue_occupancy(x, [0])",
            0,
            dmb_readout_rate*9/10
    );

    occ->SetParameter(0,dmb_readout_rate*9/10);
    occ->SetNpx(10000);
    occ->Draw();
    c2->SetLogy();
    occ->GetHistogram()->GetXaxis()->SetTitle("l1a*lct rate (kHz)");
    occ->GetHistogram()->GetYaxis()->SetTitle("events");
    occ->SetTitle("mean l1as in readout queue");




    // mean_queue_occupancy = mean_queue_occupancy (l1a_rate, dmb_readout_rate)

    // plt.xlabel("l1a_rate (kHz)")
    // plt.ylabel("overflow probability")
    // plt.semilogy(l1a_rate,overflow_probability)

    //printf ("    DMB Readout Rate = %f kHz\n"            ,  dmb_readout_rate);
    //printf ("    Mean queue occupancy = %f (%0.2f\%)\n"  ,  mean_occupancy, (mean_occupancy/maximum_occupancy*100));
    //printf ("    Mean latency = %f bx\n"                 ,  mean_occupancy * dmb_frame_cnt);
    //printf ("    Blocking Probability = %f\n"            ,  blocking_probability(maximum_occupancy,rho));

    // l1a_rate             = np.linspace(0,dmb_readout_rate/10-1,4068,endpoint=True)
    // overflow_probability = blocking_probability(maximum_occupancy,(l1a_rate/dmb_readout_rate))
    // mean_queue_occupancy = mean_queue_occupancy (l1a_rate, dmb_readout_rate)
    //
    // #plt.semilogy(l1a_rate,overflow_probability)
    // plt.grid(True)
    // plt.show()
    //
    // plt.xlabel("l1a_rate (kHz)")
    // plt.ylabel("mean queue occupancy (# of L1As)")
    // plt.semilogy(l1a_rate,mean_queue_occupancy)
    // #plt.plot(l1a_rate,mean_queue_occupancy)
    // plt.grid(True)
    // plt.show()
}
