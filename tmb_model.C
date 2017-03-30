#include <boost/circular_buffer.hpp>
#include <cmath>

#include "fence_distance.h"
#include "frame_count.h"
#include "ring_buffer.h"
#include "HistGetter.h"
#include <vector>


TRandom rando;

bool debug=0;
const int speedup = 1;
const bool low_lumi_only = 0;

const bool gen_once = 0; 

const int min_lost = 10;

const int seconds_cutoff = 10; 

const bool me11_only = 1;

const int loss_hunt_threshold = -1;


void tmb_model () {

	bool first_gen_done = 0; 

    HistGetter loss_hists;

    HistGetter misc_hists;

    rando.SetSeed(0);

    //------------------------------------------------------------------------------------------------------------------
    //
    //------------------------------------------------------------------------------------------------------------------

    int xbins = 500;
    int xmin  = 0;
    int xmax  = 50;

    int ybins = 10000;
    int ymin  = 0;
    int ymax  = 1;


    loss_hists.getOrMake2D ("h2_loss_me11_unfurled_ddr" , "Lost Event Rate (ME1/1)"                      , xbins , xmin , xmax , ybins , ymin , ymax);
    loss_hists.getOrMake2D ("h2_loss_me11_deep"         , "Lost Event Rate Deep Buffer (ME1/1)"          , xbins , xmin , xmax , ybins , ymin , ymax);
    loss_hists.getOrMake2D ("h2_loss_me11_lowl1"        , "Lost Event Rate 128bx Latency Buffer (ME1/1)" , xbins , xmin , xmax , ybins , ymin , ymax);
    loss_hists.getOrMake2D ("h2_loss_me11_gem"          , "Lost Event Rate with GEM (ME1/1)"             , xbins , xmin , xmax , ybins , ymin , ymax);

    TH2F* lostevents_arr [4] = {
        (TH2F*) loss_hists.getOrMake2D ("h2_loss_me11", "Lost Event Rate (ME1/1)",  xbins, xmin, xmax, ybins, ymin, ymax),
        (TH2F*) loss_hists.getOrMake2D ("h2_loss_me21", "Lost Event Rate (ME2/1)",  xbins, xmin, xmax, ybins, ymin, ymax),
        (TH2F*) loss_hists.getOrMake2D ("h2_loss_me31", "Lost Event Rate (ME3/1)",  xbins, xmin, xmax, ybins, ymin, ymax),
        (TH2F*) loss_hists.getOrMake2D ("h2_loss_me41", "Lost Event Rate (ME4/1)",  xbins, xmin, xmax, ybins, ymin, ymax)
    };

    TH2F* lostevents_ddr_arr [4] = {
        (TH2F*) loss_hists.getOrMake2D ("h2_loss_me11_ddr", "Lost Event Rate @ DDR (ME1/1)",  xbins, xmin, xmax, ybins, ymin, ymax),
        (TH2F*) loss_hists.getOrMake2D ("h2_loss_me21_ddr", "Lost Event Rate @ DDR (ME2/1)",  xbins, xmin, xmax, ybins, ymin, ymax),
        (TH2F*) loss_hists.getOrMake2D ("h2_loss_me31_ddr", "Lost Event Rate @ DDR (ME3/1)",  xbins, xmin, xmax, ybins, ymin, ymax),
        (TH2F*) loss_hists.getOrMake2D ("h2_loss_me41_ddr", "Lost Event Rate @ DDR (ME4/1)",  xbins, xmin, xmax, ybins, ymin, ymax)
    };


    for (int i=0; i<(loss_hists.getN2D() ); i++) {

        auto h2 = (TH2F*) loss_hists.get2D(i);

        h2 -> GetXaxis()->SetTitle("luminosity (10^{34} cm^{-2} s^{-1})");
        h2 -> GetYaxis()->SetTitle("fraction clcts lost");
        h2 -> SetStats(0);


		//-Logarithmic Binning--------------------------------------------------------------------------------------------------

  		TAxis *axis = h2->GetYaxis();

  		int bins = axis->GetNbins();

  		Axis_t from      = ymin; 
  		Axis_t to        = ymax; 
  		Axis_t width     = float(to - from) / ybins;

  		// Axis_t from      = axis->GetYmin();
  		// Axis_t to        = axis->GetYmax();
  		// Axis_t width     = (to - from) / bins;

		//printf("from=%f, to=%f, width=%f\n", from, to, width);

  		Axis_t *new_bins = new Axis_t[bins + 1];

  		for (int i = 0; i <= bins; i++) {
			float edge = 1/TMath::Power(10, from + 10*(bins-i) * width); 
  			new_bins[i]  = edge; // - (i==0 ? 0 : new_bins[i-1]); 
			//printf("bin=%i, edge=%f\n", i, edge); 
  		}
  		axis->Set(bins, new_bins);
  		delete new_bins;


    }

    //------------------------------------------------------------------------------------------------------------------
    //
    //------------------------------------------------------------------------------------------------------------------

    auto h2 = (TH2F*) misc_hists.getOrMake2D ("h2_buf_stall", "Buffer Stall Rate (ME1/1)",  xbins, xmin, xmax, ybins, ymin, 500000);

    h2-> GetXaxis()->SetTitle("luminosity (10^{34} cm^{-2} s^{-1})");
    h2-> GetYaxis()->SetTitle("Rate (Hz)");

    h2 = (TH2F*) misc_hists.getOrMake2D ("h2_queue_occupancy", "Mean L1A Queue Occupancy",  500, 0, 50, 100, 0, 100);
    h2 -> GetXaxis()->SetTitle("luminosity (10^{34} cm^{-2} s^{-1})");
    h2 -> GetYaxis()->SetTitle("#l1as in queue");

    h2 = (TH2F*) misc_hists.getOrMake2D ("h2_event_sep", "Event Separation",  500, 0, 50, 10000, 0, 10000);
    h2 -> GetXaxis()->SetTitle("luminosity (10^{34} cm^{-2} s^{-1})");
    h2 -> GetYaxis()->SetTitle("#l1as in queue");

    auto h1 = (TH1F*) misc_hists.getOrMake1D ("h1_event_sep", "Event Separation @ 1000Hz l1*lct",  25, 0, 100000);
    h1 -> GetXaxis()->SetTitle("luminosity (10^{34} cm^{-2} s^{-1})");
    h1 -> GetXaxis()->SetTitle("bx");



    #include "lumi_scaling.h"

    for (int istation=0; istation<4; istation++) {
    for (int ddr_readout=0; ddr_readout<2; ddr_readout++) {
    for (int gem_en=0; gem_en<2; gem_en++) {
    for (int deep_buffer=0; deep_buffer<2; deep_buffer++) {
    for (int low_l1=0; low_l1<2; low_l1++) {
    for (int unfurled_triads=0; unfurled_triads<2; unfurled_triads++) {

    //--------------------------------------------------------------------------------------------------------------
    // special cases
    //--------------------------------------------------------------------------------------------------------------


    //bool sim_deep_buffer = (deep_buffer==1) && (low_l1==0 && ddr_readout==0 && istation==0 && gem_en==0);

    if (me11_only && istation!=0)
        break;

    if (unfurled_triads==1 && !(low_l1==0 && ddr_readout==1 && istation==0 && gem_en==0)) {
        if (debug) printf("break-1\n");
        if (debug) printf("Station=%i, ddr_enabled=%i, gem_en=%i, deep_buffer=%i, low_l1_latency=%i\n", istation, ddr_readout, gem_en, deep_buffer, low_l1);
        break;
    }

    if (deep_buffer==1 && !(low_l1==0 && ddr_readout==0 && istation==0 && gem_en==0)) { // simulate deep buffer only at me11
        if (debug) printf("break0\n");
        if (debug) printf("Station=%i, ddr_enabled=%i, gem_en=%i, deep_buffer=%i, low_l1_latency=%i\n", istation, ddr_readout, gem_en, deep_buffer, low_l1);
        break;
    }

    if (low_l1==1 && !(istation==0 && ddr_readout==0 && deep_buffer==0))   { // simulate low l1a only at me11
        if (debug) printf("break1\n");
        if (debug) printf("Station=%i, ddr_enabled=%i, gem_en=%i, deep_buffer=%i, low_l1_latency=%i\n", istation, ddr_readout, gem_en, deep_buffer, low_l1);
        break;
    }

    if (ddr_readout==1 && istation!=0) { // only need ddr simulation for me11
        if (debug) printf("break1\n");
        if (debug) printf("Station=%i, ddr_enabled=%i, gem_en=%i, deep_buffer=%i, low_l1_latency=%i\n", istation, ddr_readout, gem_en, deep_buffer, low_l1);
        break;
    }

    if (gem_en==1 && !(istation==0 && ddr_readout==0 && deep_buffer==0 && low_l1==0 && loss_hunt_threshold<0)) { // only need gem simulation for me11
        if (debug) printf("break3\n");
        if (debug) printf("Station=%i, ddr_enabled=%i, gem_en=%i, deep_buffer=%i, low_l1_latency=%i\n", istation, ddr_readout, gem_en, deep_buffer, low_l1);
        break;
    }

    if (loss_hunt_threshold>0 && !(gem_en==0 && istation==0 && ddr_readout==0 && deep_buffer==0 && low_l1==1)) { // only need gem simulation for me11
        if (debug) printf("break4\n");
        if (debug) printf("Station=%i, ddr_enabled=%i, gem_en=%i, deep_buffer=%i, low_l1_latency=%i\n", istation, ddr_readout, gem_en, deep_buffer, low_l1);
        break;
    }

	if (gen_once && first_gen_done) 
		break; 

	first_gen_done=1; 

    printf("Station=%i, ddr_enabled=%i, gem_en=%i, deep_buffer=%i, low_l1_latency=%i, unfurled_triads=%i\n", istation, ddr_readout, gem_en, deep_buffer, low_l1, unfurled_triads);

    //--------------------------------------------------------------------------------------------------------------

    int n_cfebs = (istation==0) ? 7 : 5;

    float lumi_rate          = lumi_rate_arr[istation];
    float lumi_rate_pre      = lumi_rate_pre_arr[istation];
    float l1a_match_fraction = lumi_rate/lumi_rate_pre;

    int l1a_latency = low_l1 ? 128 : 500;

    int frame_cnt = frame_count (n_cfebs, gem_en, unfurled_triads);
    double dmb_readout_rate = 1000000.0/(frame_cnt*25); // DMB readout time, kHz

    printf("Frame count          : %i \n", frame_cnt);
    printf("Maximumn readout rate: %f kHz\n", dmb_readout_rate);

    //------------------------------------------------------------------------------------------------------------------

       for (int lum = 1; // luminosity 5 = (0.5*10^34)
             lum < 500 / (low_lumi_only ? 3 : 1);  // e.g. 500 = 50 * 10^{34}
             lum+=1*speedup) {

        int buf_depth = deep_buffer ? 4096 : 2048;

        ringBuffer ring_buffer(buf_depth);


        float mean_rate = lumi_rate_pre * float (lum) * 1000; // convert luminosity to preCLCT rate, based on station

        int    mean_separation = (int) ( (1./mean_rate) * ( pow(10,9) ) / (25.)); // bunches

        std::vector<int> fence_queue;
        std::vector<int> pretrig_queue;
        std::vector<bool> short_ro_queue;


        long long unsigned bx = 0;
        long long unsigned last_lost_bx = 0;
        long long unsigned pretrig_bx = 256;

        bool rd_started  = 0;
        bool rd_finished = 1;
        long long unsigned rd_finish_bx = -1;

        int distance = buf_depth;

        int preclct_lost        = 0;
        int preclct_cnt         = 0;
        int l1a_match_cnt       = 0;
        int num_stalls          = 0;

		int seconds = 0; 
        unsigned long long sum_queue_occupancy = 0;

        bool keep_stalled_events = 0;

        bool buf_stalled = 0;

        while (true) {

            //----------------------------------------------------------------------------------------------------------
            // Fence Distance
            //----------------------------------------------------------------------------------------------------------

            if (fence_queue.size() > 0) {
                // make sure that the fence buffer has space--- otherwise we should stall here
                distance = fence_distance(fence_queue[0], ring_buffer.wr_adr(), buf_depth);
            }
            else {
                buf_stalled = 0;
                distance = l1a_latency;
                //if (debug) printf("bx=%llu queue empty\n", bx);
            }

            //----------------------------------------------------------------------------------------------------------
            // Buffer write control
            //----------------------------------------------------------------------------------------------------------

            if (distance==0) {
                buf_stalled = 1;
                num_stalls++;
            }
            else if (buf_stalled)
            {
                if (distance > 64)
                    buf_stalled = 0;
            }
            else {
                // we get a readout signal l1a_delay later than the event pretrigger time
                ring_buffer.inc_wr_adr();
            }


            //----------------------------------------------------------------------------------------------------------
            // Pretrig Inserter
            //----------------------------------------------------------------------------------------------------------

            if (bx==pretrig_bx) {

                // calculate bx time for next preclct trigger

                float sep = -log(1-rando.Uniform(1)) / (mean_rate*25.0 / pow(10,9));
                float rate = mean_rate;

                misc_hists.get2D("h2_event_sep") -> Fill(rate,sep);

                if (mean_rate==1000) {
                    misc_hists.get1D("h1_event_sep") -> Fill(sep);
                }

                int separation = (int) sep;

                separation = separation==0 ? 1 : separation;

                pretrig_bx = bx + separation;

                //------------------------------------------------------------------------------------------------------
                // Event Accepted, Buffer Full
                //------------------------------------------------------------------------------------------------------
                if (buf_stalled) {
                    if (debug) printf("Event lost (buf stalled)!\n");
                    if (debug) printf("bx=%lld\n",bx);

                    if (keep_stalled_events) {
                        // keep the event but mark for short-mode readout
                        pretrig_queue.push_back (pretrig_bx);
                        short_ro_queue.push_back (1);
                    }

                    preclct_lost++;

                    if (loss_hunt_threshold>0) {
                        printf("pretrig dropped, sep=%4i, bx=%12i, delta=%12i, lost=%3i, total=%8i, lumi=%5.2f * 10^{34}\n", separation ,bx, bx-last_lost_bx, preclct_lost, preclct_cnt, float(lum)/10.);
                        last_lost_bx=bx;
                    }
                }
                //------------------------------------------------------------------------------------------------------
                // Event Accepted, Buffer Available
                //------------------------------------------------------------------------------------------------------
                else {
                    pretrig_queue.push_back (pretrig_bx);
                    short_ro_queue.push_back (0);
                    if (debug) printf("pretrig pushed   @ %4i, events in l1a queue=%i\n",ring_buffer.wr_adr(), (int) fence_queue.size());
                }

                preclct_cnt ++;
            }

            //----------------------------------------------------------------------------------------------------------
            // L1A Queue Inserter
            //----------------------------------------------------------------------------------------------------------

            if (pretrig_queue.size() > 0) {

                if (bx==pretrig_queue[0]+l1a_latency) {

                    // accepted pretrigs only have some % chance to readout (l1a match)
                    // what % is this ???

                    if (rando.Uniform(1) < l1a_match_fraction) { // fraction of preTrigs that become CLCT*L1A

                        l1a_match_cnt++;

                        // push pretrig_bx to fence queue for raw hits readout to DMB

                        int pretrig_adr = (buf_depth-1) & (ring_buffer.wr_adr()-l1a_latency);

                        fence_queue.push_back (pretrig_adr);

                        sum_queue_occupancy += fence_queue.size();

                        if (debug) printf("event pushed     @ %4i, events in l1a queue=%i\n",pretrig_adr, (int) fence_queue.size());

                    }

                    pretrig_queue.erase(pretrig_queue.begin());
                }

            }


            //----------------------------------------------------------------------------------------------------------
            // Buffer read control
            //----------------------------------------------------------------------------------------------------------

            if (fence_queue.size() > 0) { // # of l1as for readout

                // we have events in queue but they haven't been readout yet
                if (rd_started==0) {

                    // mark readout as busy (assert dmb_fifo_wr)
                    rd_started=1;

                    int cnt;
                    // wr_buffer was not available; use short readout mode
                    if (short_ro_queue[0])
                        cnt = 12;
                    // full readout
                    else {
                        cnt = frame_cnt;
                    }

                    ring_buffer.rd_adr(fence_queue[0]);
                    short_ro_queue.erase(short_ro_queue.begin());

                    //-DDR Mode-----------------------------------------------------------------------------------------
                    if (ddr_readout)
                        cnt = cnt/2;

                    rd_finish_bx = bx + cnt;

                    if (debug) printf("readout started  @ %4llu, events in l1a queue=%i (cnt=%i)\n", fence_queue[0], (int) fence_queue.size(), cnt);

                }

                // finished this readout
                else if (bx==rd_finish_bx) {

                    rd_started=0;

                    // pop from queue after readout
                    fence_queue.erase(fence_queue.begin());

                    if (debug) printf("readout finished @ %4llu, events in l1a queue=%i\n", ring_buffer.rd_adr(), (int) fence_queue.size());
                }
                else {
                    ring_buffer.inc_rd_adr();
                }

            }

            //----------------------------------------------------------------------------------------------------------
            // Bx Counter
            //----------------------------------------------------------------------------------------------------------
            bx++;

			if (bx%40000000==0) {
				seconds++; 
				printf("    seconds=%i, lost=%i (%e)\n", seconds, preclct_lost, preclct_lost/float(preclct_cnt)); 
			}

            if (debug && bx>2048)
                break;

            //----------------------------------------------------------------------------------------------------------
            // Printout & Fill Histograms
            //----------------------------------------------------------------------------------------------------------
            if (seconds>=seconds_cutoff || (preclct_cnt >= 50000/(speedup/4.0) && loss_hunt_threshold<0 && preclct_lost >= min_lost)) {

                float preclct_rate = (float(preclct_cnt)   * pow(10,9) / (25.*bx)); // Hz
                float readout_rate = (l1a_match_cnt * pow(10,9) / (25.*bx));        // Hz
                float stall_rate   = (num_stalls * pow(10,9) / (25.*bx));           // Hz

                float lost_event_frac = float (preclct_lost)/float(preclct_cnt); // fraction

                float luminosity = preclct_rate * l1a_match_fraction / lumi_rate / 10000.;
                printf("preclct_rate=%7.1f kHz, readout_rate=%5.1f kHz, preclct_cnt=%7i, pretrigs_dropped=%5i, mean_occupancy=%4.2f, luminosity=%4.1f *10^34, lhc_seconds=%i\n",
                        preclct_rate/1000.,
                        readout_rate/1000.,
                        preclct_cnt,
                        preclct_lost, // lost_event_frac * 100.0 * 1000,
                        sum_queue_occupancy/double(l1a_match_cnt),
                        luminosity,
						seconds
                      );


                // low l1a latency mode
                if (low_l1 && istation==0 && !ddr_readout && !deep_buffer && !gem_en) {
                    loss_hists.get2D("h2_loss_me11_lowl1") -> Fill(luminosity, lost_event_frac);
                }
                // deep buffers
                else if (deep_buffer==1 && !low_l1 && !ddr_readout && istation==0 && !gem_en) {
                    loss_hists.get2D("h2_loss_me11_deep") -> Fill(luminosity, lost_event_frac);
                }
                else if (ddr_readout && !low_l1 && !gem_en && !unfurled_triads) {
                        (lostevents_ddr_arr[istation])  -> Fill(luminosity, lost_event_frac);
                }
                else if (gem_en && !low_l1 && !ddr_readout) {
                        loss_hists.get2D("h2_loss_me11_gem") -> Fill(luminosity, lost_event_frac);
                }
                else if (unfurled_triads && ddr_readout && !gem_en && !low_l1) {
                        loss_hists.get2D("h2_loss_me11_unfurled_ddr") -> Fill(luminosity, lost_event_frac);
                }
                else if (!gem_en && !low_l1 && !gem_en && !ddr_readout) {

                        (lostevents_arr[istation])  -> Fill(luminosity, lost_event_frac);

                        if (istation==0) {
                            misc_hists.get2D("h2_buf_stall") -> Fill(luminosity, stall_rate);
                            misc_hists.get2D("h2_queue_occupancy") -> Fill (luminosity, sum_queue_occupancy/double(l1a_match_cnt));
                        }
                }


                break;
            }
            else if (loss_hunt_threshold>0) {
                if (preclct_lost > loss_hunt_threshold) {
                    printf("\nDONE\n");
                    return;
                }
            }

        } // while true

    } // luminosity

    }
    }
    }
    }
    }
    }

    loss_hists.write("tmb_model_loss.root");
    misc_hists.write("tmb_model_misc.root");

    if (loss_hunt_threshold<0) {
        gSystem->Exec("hadd -f tmb_model.root tmb_model_loss.root tmb_model_misc.root");
    }
}
