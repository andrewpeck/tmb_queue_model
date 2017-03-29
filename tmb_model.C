#include <boost/circular_buffer.hpp>
#include <cmath>

#include "fence_distance.h"
#include "frame_count.h"
#include "ring_buffer.h"
#include <vector>



TRandom rando;

bool debug=0;
const int speedup = 2;
const bool low_lumi_only = 0;

const int loss_hunt_threshold = -1;

void tmb_model () {

    TFile* hfile = new TFile("tmb_model_tmp.root","RECREATE","TMB Queue Model");

    rando.SetSeed(0);

    TH2F* th2f_occupancy = new TH2F ("Mean Queue Occupancy", "Mean L1A Queue Occupancy",  500, 0, 50, 100, 0, 100);
    th2f_occupancy -> GetXaxis()->SetTitle("luminosity (10^34 cm-2 s-1)");
    th2f_occupancy -> GetYaxis()->SetTitle("#l1as in queue");

    //------------------------------------------------------------------------------------------------------------------
    //
    //------------------------------------------------------------------------------------------------------------------

    int xbins = 500;
    int xmin  = 0;
    int xmax  = 50;

    int ybins = 10000;
    int ymin  = 0;
    int ymax  = 1;

    TH2F* th2f_lostevents_me11 = new TH2F ("Lost Events ME11", "Lost Event Rate (ME1/1)",  xbins, xmin, xmax, ybins, ymin, ymax);
    th2f_lostevents_me11 -> GetXaxis()->SetTitle("luminosity (10^34 cm-2 s-1)");
    th2f_lostevents_me11 -> GetYaxis()->SetTitle("fraction clcts lost");
    th2f_lostevents_me11 -> SetStats(0);

    TH2F* th2f_lostevents_ddr_me11 = new TH2F ("Lost Events DDR ME11", "Lost Event Rate (ME1/1)",  xbins, xmin, xmax, ybins, ymin, ymax);
    th2f_lostevents_ddr_me11 -> GetXaxis()->SetTitle("luminosity (10^34 cm-2 s-1)");
    th2f_lostevents_ddr_me11 -> GetYaxis()->SetTitle("fraction clcts lost");
    th2f_lostevents_ddr_me11 -> SetStats(0);

    TH2F* th2f_lostevents_deep_me11 = new TH2F ("Lost Events Deep Buffer ME11", "Lost Event Rate Deep Buffer (ME1/1)",  xbins, xmin, xmax, ybins, ymin, ymax);
    th2f_lostevents_deep_me11 -> GetXaxis()->SetTitle("luminosity (10^34 cm-2 s-1)");
    th2f_lostevents_deep_me11 -> GetYaxis()->SetTitle("fraction clcts lost");
    th2f_lostevents_deep_me11 -> SetStats(0);

    TH2F* th2f_lostevents_low_l1_me11 = new TH2F ("Lost Events LOWL1 ME11", "Lost Event Rate 128bx Latency Buffer (ME1/1)",  xbins, xmin, xmax, ybins, ymin, ymax);
    th2f_lostevents_low_l1_me11 -> GetXaxis()->SetTitle("luminosity (10^34 cm-2 s-1)");
    th2f_lostevents_low_l1_me11 -> GetYaxis()->SetTitle("fraction clcts lost");
    th2f_lostevents_low_l1_me11 -> SetStats(0);

    TH2F* th2f_lostevents_gem_me11 = new TH2F ("Lost Events GEM ME11", "Lost Event Rate with GEM (ME1/1)",  xbins, xmin, xmax, ybins, ymin, ymax);
    th2f_lostevents_gem_me11 -> GetXaxis()->SetTitle("luminosity (10^34 cm-2 s-1)");
    th2f_lostevents_gem_me11 -> GetYaxis()->SetTitle("fraction clcts lost");
    th2f_lostevents_gem_me11 -> SetStats(0);

    TH2F* th2f_lostevents_me21 = new TH2F ("Lost Events ME21", "Lost Event Rate (ME2/1)",  xbins, xmin, xmax, ybins, ymin, ymax);
    th2f_lostevents_me21 -> GetXaxis()->SetTitle("luminosity (10^34 cm-2 s-1)");
    th2f_lostevents_me21 -> GetYaxis()->SetTitle("fraction clcts lost");
    th2f_lostevents_me21 -> SetStats(0);

    TH2F* th2f_lostevents_ddr_me21 = new TH2F ("Lost Events DDR ME21", "Lost Event Rate (ME2/1)",  xbins, xmin, xmax, ybins, ymin, ymax);
    th2f_lostevents_ddr_me21 -> GetXaxis()->SetTitle("luminosity (10^34 cm-2 s-1)");
    th2f_lostevents_ddr_me21 -> GetYaxis()->SetTitle("fraction clcts lost");
    th2f_lostevents_ddr_me21 -> SetStats(0);

    TH2F* th2f_lostevents_me31 = new TH2F ("Lost Events ME31", "Lost Event Rate (ME3/1)",  xbins, xmin, xmax, ybins, ymin, ymax);
    th2f_lostevents_me31 -> GetXaxis()->SetTitle("luminosity (10^34 cm-2 s-1)");
    th2f_lostevents_me31 -> GetYaxis()->SetTitle("fraction clcts lost");
    th2f_lostevents_me31 -> SetStats(0);

    TH2F* th2f_lostevents_ddr_me31 = new TH2F ("Lost Events DDR ME31", "Lost Event Rate (ME3/1)",  xbins, xmin, xmax, ybins, ymin, ymax);
    th2f_lostevents_ddr_me31 -> GetXaxis()->SetTitle("luminosity (10^34 cm-2 s-1)");
    th2f_lostevents_ddr_me31 -> GetYaxis()->SetTitle("fraction clcts lost");
    th2f_lostevents_ddr_me31 -> SetStats(0);

    TH2F* th2f_lostevents_me41 = new TH2F ("Lost Events ME41", "Lost Event Rate (ME4/1)",  xbins, xmin, xmax, ybins, ymin, ymax);
    th2f_lostevents_me41 -> GetXaxis()->SetTitle("luminosity (10^34 cm-2 s-1)");
    th2f_lostevents_me41 -> GetYaxis()->SetTitle("fraction clcts lost");
    th2f_lostevents_me41 -> SetStats(0);

    TH2F* th2f_lostevents_ddr_me41 = new TH2F ("Lost events DDR ME41", "Lost Event Rate (ME4/1)",  xbins, xmin, xmax, ybins, ymin, ymax);
    th2f_lostevents_ddr_me41 -> GetXaxis()->SetTitle("luminosity (10^34 cm-2 s-1)");
    th2f_lostevents_ddr_me41 -> GetYaxis()->SetTitle("fraction clcts lost");
    th2f_lostevents_ddr_me41 -> SetStats(0);

    TH2F* lostevents_arr [4] = {th2f_lostevents_me11, th2f_lostevents_me21, th2f_lostevents_me31,  th2f_lostevents_me41};
    TH2F* lostevents_ddr_arr [4] = {th2f_lostevents_ddr_me11, th2f_lostevents_ddr_me21, th2f_lostevents_ddr_me31,  th2f_lostevents_ddr_me41};


    //------------------------------------------------------------------------------------------------------------------
    //
    //------------------------------------------------------------------------------------------------------------------


    TProfile* tprof_numstalls = new TProfile ("Buffer Stall Rate", "Number of Buffer Stalls",  225, 0, 0, "");
    tprof_numstalls -> GetXaxis()->SetTitle("luminosity (10^34 cm-2 s-1)");
    tprof_numstalls -> GetYaxis()->SetTitle("Buffer stall rate (Hz)");

    TProfile* th2f_meansep = new TProfile ("Mean Event Separation", "Mean Event Separation",  225, 0, 200000, "");
    TH1F* th1f_sep = new TH1F ("Event Separation", "Event Separation @ 1000Hz l1*lct",  25, 0, 100000);

    #include "lumi_scaling.h"

    for (int istation=0; istation<4; istation++) {

    for (int gem_en=0; gem_en<2; gem_en++) {

    for (int ddr_readout=0; ddr_readout<2; ddr_readout++) {

    for (int deep_buffer=0; deep_buffer<2; deep_buffer++) {

    for (int low_l1=0; low_l1<2; low_l1++) {

    for (int dummy=0; dummy<2; dummy++) { // no idea why this is necessary--- something weird is going on with these looops

    //--------------------------------------------------------------------------------------------------------------
    // special cases
    //--------------------------------------------------------------------------------------------------------------

    printf("Station=%i, ddr_enabled=%i, gem_en=%i, deep_buffer=%i, low_l1_latency=%i\n", istation, ddr_readout, gem_en, deep_buffer, low_l1);

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

    //--------------------------------------------------------------------------------------------------------------

    int n_cfebs = (istation==0) ? 7 : 5;

    float lumi_rate          = lumi_rate_arr[istation];
    float lumi_rate_pre      = lumi_rate_pre_arr[istation];
    float l1a_match_fraction = lumi_rate/lumi_rate_pre;

    int l1a_latency = low_l1 ? 128 : 500;

    printf("\nStation=%i, ddr_enabled=%i, gem_en=%i, deep_buffer=%i, low_l1_latency=%i\n", istation, ddr_readout, gem_en, deep_buffer, low_l1);

    int frame_cnt = frame_count (n_cfebs, gem_en);
    double dmb_readout_rate = 1000000.0/(frame_cnt*25); // DMB readout time, kHz

    printf("Frame count          : %i \n", frame_cnt);
    printf("Maximumn readout rate: %f kHz\n", dmb_readout_rate);

    //------------------------------------------------------------------------------------------------------------------

       for (int lum = 1; // luminosity
             lum < 500 / (low_lumi_only ? 3 : 1);  // e.g. 500 = 50 * 10^34
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

                th2f_meansep->Fill(rate, sep);

                if (mean_rate==1000) {
                    th1f_sep->Fill(sep);
                    //printf("rate=%4f, sep=%f\n",rate, sep);
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
                        printf("pretrig dropped, sep=%4i, bx=%12i, delta=%12i, lost=%3i, total=%8i, lumi=%5.2f * 10^34\n", separation ,bx, bx-last_lost_bx, preclct_lost, preclct_cnt, float(lum)/10.);
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

                    if (rando.Uniform(1) < l1a_match_fraction) {

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

            if (debug && bx>2048)
                break;

            //----------------------------------------------------------------------------------------------------------
            // Printout & Fill Histograms
            //----------------------------------------------------------------------------------------------------------
            if (preclct_cnt >= 50000/(speedup/4.0) && loss_hunt_threshold<0) {

                float preclct_rate = (float(preclct_cnt)   * pow(10,9) / (25.*bx)); // Hz
                float readout_rate = (l1a_match_cnt * pow(10,9) / (25.*bx));        // Hz
                float stall_rate   = (num_stalls * pow(10,9) / (25.*bx));           // Hz

                float lost_event_frac = float (preclct_lost)/float(preclct_cnt); // fraction

                float luminosity = preclct_rate * l1a_match_fraction / lumi_rate / 10000.;
                printf("preclct_rate=%7.1f kHz, readout_rate=%5.1f kHz, preclct_cnt=%7i, pretrigs_dropped=%5i, mean_occupancy=%4.2f, luminosity=%4.1f *10^34\n",
                        preclct_rate/1000.,
                        readout_rate/1000.,
                        preclct_cnt,
                        preclct_lost, // lost_event_frac * 100.0 * 1000,
                        sum_queue_occupancy/double(l1a_match_cnt),
                        luminosity

                      );


                // low l1a latency mode
                if (low_l1 && istation==0 && !ddr_readout && !deep_buffer && !gem_en) {
                    th2f_lostevents_low_l1_me11 -> Fill(luminosity, lost_event_frac);
                }
                // deep buffers
                else if (deep_buffer==1 && !low_l1 && !ddr_readout && istation==0 && !gem_en) {
                    th2f_lostevents_deep_me11 -> Fill(luminosity, lost_event_frac);
                }
                else if (ddr_readout && !low_l1 && !gem_en) {
                        (lostevents_ddr_arr[istation])  -> Fill(luminosity, lost_event_frac);
                }
                else if (gem_en && !low_l1 && !ddr_readout) {
                        (th2f_lostevents_gem_me11)  -> Fill(luminosity, lost_event_frac);
                }
                else if (!gem_en && !low_l1 && !gem_en && !ddr_readout) {

                        (lostevents_arr[istation])  -> Fill(luminosity, lost_event_frac);

                        if (istation==0) {
                            th2f_occupancy   -> Fill(luminosity, sum_queue_occupancy/double(l1a_match_cnt));
                            tprof_numstalls  -> Fill(luminosity, stall_rate);
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

    hfile->Write();
    hfile->Close();

    if (loss_hunt_threshold<0)
    gSystem->Exec("mv tmb_model_tmp.root tmb_model.root");
}
