#include <boost/circular_buffer.hpp>
#include <cmath>
#include <vector>
#include <TMath.h>
#include <TRandom.h>
#include <algorithm>
#include <TSystem.h>
#include <TApplication.h>

#include "fence_distance.h"
#include "frame_count.h"
#include "ring_buffer.h"
#include "HistGetter.h"


TRandom rando;

bool debug=0;
bool debug_core=0;

bool no_bx_skipping=0;

// sim options
const int speedup = 1;
const bool low_lumi_only = 0;
const bool gen_once = 0;
const int min_lost = 5;
const int seconds_cutoff = 240;
const bool me11_only = 1;

const int loss_hunt_threshold = 0;


int tmb_model () {

    bool first_gen_done = 0;

    HistGetter loss_hists;

    HistGetter misc_hists;

    //------------------------------------------------------------------------------------------------------------------
    //
    //------------------------------------------------------------------------------------------------------------------

    int xbins = 500;
    int xmin  = 0;
    int xmax  = 50;

    int ybins = 10000;
    int ymin  = 0;
    int ymax  = 1;


    loss_hists.getOrMake2D ("h2_loss_me11_unfurled_ddr" , "Lost Event Rate Unfurled Triads (ME1/1)"      , xbins , xmin , xmax , ybins , ymin , ymax);
    loss_hists.getOrMake2D ("h2_loss_me11_deep"         , "Lost Event Rate Deep Buffer (ME1/1)"          , xbins , xmin , xmax , ybins , ymin , ymax);
    loss_hists.getOrMake2D ("h2_loss_me11_deep_ddr"     , "Lost Event Rate Deep Buffer DDR (ME1/1)"      , xbins , xmin , xmax , ybins , ymin , ymax);
    loss_hists.getOrMake2D ("h2_loss_me11_lowl1"        , "Lost Event Rate 128bx Latency Buffer (ME1/1)" , xbins , xmin , xmax , ybins , ymin , ymax);
    loss_hists.getOrMake2D ("h2_loss_me11_gem"          , "Lost Event Rate with GEM (ME1/1)"             , xbins , xmin , xmax , ybins , ymin , ymax);
    loss_hists.getOrMake2D ("h2_loss_me11_gem_ddr"      , "Lost Event Rate with GEM + DDR(ME1/1)"        , xbins , xmin , xmax , ybins , ymin , ymax);

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


    for (unsigned i=0; i<(loss_hists.getN2D() ); i++) {

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
        delete [] new_bins;


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
    h2 -> GetYaxis()->SetTitle("separation between readouts, me11 baseline");

    h2 = (TH2F*) misc_hists.getOrMake2D ("h2_pretrig_sep", "Event Separation",  500, 0, 50, 10000, 0, 10000);
    h2 -> GetXaxis()->SetTitle("luminosity (10^{34} cm^{-2} s^{-1})");
    h2 -> GetYaxis()->SetTitle("separation between pretrigs, me11 baseline");

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

    if (gem_en==1 && !(istation==0 && deep_buffer==0 && low_l1==0 && loss_hunt_threshold<=0)) { // only need gem simulation for me11
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

    unsigned long long l1a_latency = (low_l1 || debug_core) ? 128 : 500;

    int frame_cnt = frame_count (n_cfebs, gem_en, unfurled_triads);
    double dmb_readout_rate = 1000000.0/(frame_cnt*25); // DMB readout time, kHz

    printf("Frame count          : %i \n", frame_cnt);
    printf("Maximumn readout rate: %f kHz\n", dmb_readout_rate);

    //------------------------------------------------------------------------------------------------------------------

       for (int lum = ((debug_core) ? 495 : 1); // luminosity 5 = (0.5*10^34)
             lum < (low_lumi_only ? 150 : 500);  // e.g. 500 = 50 * 10^{34}
             lum+=speedup) {

        //printf("lum=%f\n", float(lum));

        int buf_depth = deep_buffer ? 4096 : 2048;

        ringBuffer ring_buffer(buf_depth);


        float mean_rate = lumi_rate_pre * float (lum) * 1000; // convert luminosity to preCLCT rate, based on station

        //int    mean_separation = (int) ( (1./mean_rate) * ( pow(10,9) ) / (25.)); // bunches

        std::vector<int> fence_queue;
        std::vector<int> pretrig_queue;
        std::vector<bool> short_ro_queue;


        long long unsigned bx = 0;
        long long unsigned delta_bx = 1;
        long long unsigned last_lost_bx = 0;
        long long unsigned pretrig_bx = 256;

        bool rd_started  = 0;
        long long unsigned rd_finish_bx = -1;

        int distance = buf_depth;

        int preclct_lost        = 0;
        int preclct_cnt         = 0;
        int l1a_match_cnt       = 0;
        int num_stalls          = 0;

        double seconds = 0;

        unsigned long long sum_queue_occupancy = 0;

        bool keep_stalled_events = 0;

        bool buf_stalled = 0;

        while (true) {


            //if (bx>255)
                //if (debug_core) printf("                 @ bx=%4i, wr_adr=%4i, rd_adr=%4i, events in l1a queue=%i\n", bx, ring_buffer.wr_adr(), ring_buffer.rd_adr(), (int) fence_queue.size());

            // need to check bx
            // (1) when fence_queue[0] changes
            // in principle... we should be able to advance the bx# to the minimum of either
            //     (a)  the bx where a L1A is finished with readout
            //     (b)  the bx where a pretrigger is spawned
            // when advancing the bx number by >1 , we also need to advance the ring buffer wr_adr by the same amount


            // pretrig pushed   @ bx= 256, wr_adr= 256, rd_adr=   0, events in pretrig queue=1
            // event pushed     @ bx= 384, wr_adr= 384, rd_adr=   0, events in l1a queue=1
            // readout started  @ bx= 384, wr_adr= 384, rd_adr= 256, events in l1a queue=1
            // pretrig pushed   @ bx= 434, wr_adr= 434, rd_adr= 306, events in pretrig queue=1
            // pretrig pushed   @ bx= 481, wr_adr= 481, rd_adr= 353, events in pretrig queue=2
            // event pushed     @ bx= 562, wr_adr= 562, rd_adr= 434, events in l1a queue=2
            // pretrig pushed   @ bx= 572, wr_adr= 572, rd_adr= 444, events in pretrig queue=2
            // event pushed     @ bx= 609, wr_adr= 609, rd_adr= 481, events in l1a queue=3
            // event pushed     @ bx= 700, wr_adr= 700, rd_adr= 572, events in l1a queue=4
            // pretrig pushed   @ bx= 760, wr_adr= 760, rd_adr= 632, events in pretrig queue=1
            // event pushed     @ bx= 888, wr_adr= 888, rd_adr= 760, events in l1a queue=5
            // readout finished @ bx= 936, wr_adr= 936, rd_adr= 808, events in l1a queue=4
            // readout started  @ bx= 937, wr_adr= 937, rd_adr= 434, events in l1a queue=4
            // pretrig pushed   @ bx=1488, wr_adr=1488, rd_adr= 985, events in pretrig queue=1
            // readout finished @ bx=1489, wr_adr=1489, rd_adr= 986, events in l1a queue=3
            // readout started  @ bx=1490, wr_adr=1490, rd_adr= 481, events in l1a queue=3
            // pretrig pushed   @ bx=1534, wr_adr=1534, rd_adr= 525, events in pretrig queue=2
            // event pushed     @ bx=1616, wr_adr=1616, rd_adr= 607, events in l1a queue=4
            // event pushed     @ bx=1662, wr_adr=1662, rd_adr= 653, events in l1a queue=5
            // pretrig pushed   @ bx=1833, wr_adr=1833, rd_adr= 824, events in pretrig queue=1
            // pretrig pushed   @ bx=1896, wr_adr=1896, rd_adr= 887, events in pretrig queue=2
            // event pushed     @ bx=1961, wr_adr=1961, rd_adr= 952, events in l1a queue=6
            // event pushed     @ bx=2024, wr_adr=2024, rd_adr=1015, events in l1a queue=7
            // readout finished @ bx=2042, wr_adr=2042, rd_adr=1033, events in l1a queue=6
            // readout started  @ bx=2043, wr_adr=2043, rd_adr= 572, events in l1a queue=6


            unsigned long long delta_bx      = 0;
            unsigned long long last_event_bx = 0;
            unsigned long long next_bx_rqst  = -1;


            //----------------------------------------------------------------------------------------------------------
            // Fence Distance
            //----------------------------------------------------------------------------------------------------------

            // this can be calculated any time.

            if (fence_queue.size() > 0) {
                // make sure that the fence buffer has space--- otherwise we should stall here
                distance = fence_distance(fence_queue[0], ring_buffer.wr_adr(), buf_depth);
                //printf ("distance=%i (queue=%i - wr_adr=%i)\n", distance, fence_queue[0], ring_buffer.wr_adr());
            }
            else {
                buf_stalled = 0;
                distance = l1a_latency;
                //if (debug_core) printf("bx=%llu queue empty\n", bx);
            }

            //----------------------------------------------------------------------------------------------------------
            // Buffer write control
            //----------------------------------------------------------------------------------------------------------

            // need to increment the write address by delta BX

            // this only needs to be checked if something new pops onto the fence queue

            if (distance==0) {
                buf_stalled = 1;
                num_stalls++;
            }
            else if (buf_stalled)
            {
                if (distance > 64)
                    buf_stalled = 0;
            }


            //----------------------------------------------------------------------------------------------------------
            // Pretrig Inserter
            //----------------------------------------------------------------------------------------------------------

            // this only needs to be entered at the obvious
            //         (1)  bx==pretrig_bx

            if (bx==pretrig_bx) {

                //------------------------------------------------------------------------------------------------------
                // Calculate BX for *next* pretrigger
                //------------------------------------------------------------------------------------------------------


                rando.SetSeed(0);

                float sep = -log(1-rando.Uniform(1)) / (mean_rate*25.0 / pow(10,9));

                if (istation==0 && !ddr_readout && !gem_en && !low_l1) {
                    misc_hists.get2D("h2_pretrig_sep") -> Fill(float(lum)/10,sep);
                    //printf("filling pretrig sep Fill(%f, %i)\n", float(lum)/10, sep);
                }

                int separation = (int) sep;

                separation = separation==0 ? 1 : separation;

                pretrig_bx = bx + separation;

                // need processing at pretrig bx, to roll for the next pretrigger


                //------------------------------------------------------------------------------------------------------
                // Event Accepted, Buffer Full for *this* pretrigger
                //------------------------------------------------------------------------------------------------------

                if (buf_stalled) {
                    //if (debug_core) printf("Event lost (buf stalled)!\n");
                    //if (debug_core) printf("bx=%lld\n",bx);

                    if (keep_stalled_events) {
                        // keep the event but mark for short-mode readout
                        pretrig_queue.push_back (bx);
                        short_ro_queue.push_back (1);
                    }

                    preclct_lost++;

                    if (debug_core) printf("pretrig dropped, sep=%4i, bx=%12llu, delta=%12llu, lost=%3i, total=%8i, lumi=%5.2f * 10^{34}\n", separation ,bx, bx-last_lost_bx, preclct_lost, preclct_cnt, float(lum)/10.);

                    if (loss_hunt_threshold>0) {
                        last_lost_bx=bx;
                    }
                }

                //------------------------------------------------------------------------------------------------------
                // Event Accepted, Buffer Available
                //------------------------------------------------------------------------------------------------------

                else {
                    pretrig_queue.push_back (bx);
                    short_ro_queue.push_back (0);
                    if (debug_core) printf("pretrig pushed   @ bx=%4llu, wr_adr=%4i, rd_adr=%4i, distance=%i, events in pretrig queue=%i\n", bx, ring_buffer.wr_adr(), ring_buffer.rd_adr(), distance, (int) pretrig_queue.size());
                }

                preclct_cnt ++;

            }

            //----------------------------------------------------------------------------------------------------------
            // L1A Queue Inserter
            //----------------------------------------------------------------------------------------------------------

            if (pretrig_queue.size() > 0) {

                // this only needs to be checked  on the obvious condition
                //         (1) if (bx==pretrig_queue[0]+l1a_latency)

                if (bx==pretrig_queue[0]+l1a_latency) {

                    // accepted pretrigs only have some % chance to readout (l1a match)
                    // what % is this ???

                    if (debug_core || (rando.Uniform(1) < l1a_match_fraction)) { // fraction of preTrigs that become CLCT*L1A

                        l1a_match_cnt++;

                        // push pretrig_bx to fence queue for raw hits readout to DMB

                        int pretrig_adr = (buf_depth-1) & (ring_buffer.wr_adr()-l1a_latency);

                        if (istation==0 && !ddr_readout && !gem_en && !low_l1) {
                            misc_hists.get2D("h2_event_sep") -> Fill(float(lum)/10, bx-last_event_bx);
                            //printf("filling event sep Fill(%f, %i)\n", float(lum)/10, bx-last_event_bx);
                            last_event_bx = bx;
                        }

                        fence_queue.push_back (pretrig_adr);

                        sum_queue_occupancy += fence_queue.size();

                        if (debug_core) printf("event pushed     @ bx=%4llu, wr_adr=%4i, rd_adr=%4i, distance=%i, events in l1a queue=%i\n", bx, ring_buffer.wr_adr(), ring_buffer.rd_adr(), distance, (int) fence_queue.size());

                    }

                    pretrig_queue.erase(pretrig_queue.begin());
                }
            }


            //----------------------------------------------------------------------------------------------------------
            // Buffer read control
            //----------------------------------------------------------------------------------------------------------

            if (fence_queue.size() > 0) { // # of l1as for readout

                //  this condition will only happen on
                //    (1) on the bx after we finish a readout
                //    (2) on the bx that a pretrig gets inserted into the queue (bx==pretrig_queue[0]+l1a_latency)

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


                    if (debug_core) printf("readout started  @ bx=%4llu, wr_adr=%4i, rd_adr=%4i, distance=%i, events in l1a queue=%i\n", bx, ring_buffer.wr_adr(), ring_buffer.rd_adr(), distance, (int) fence_queue.size());

                }

                //  this condition will only happen on
                //    (1) bx==rd_finish_bx

                // finished this readout
                if (rd_started==1 && bx==rd_finish_bx) {

                    rd_started=0;

                    // pop from queue after readout
                    fence_queue.erase(fence_queue.begin());

                    if (debug_core) printf("readout finished @ bx=%4llu, wr_adr=%4i, rd_adr=%4i, distance=%i, events in l1a queue=%i\n", bx, ring_buffer.wr_adr(), ring_buffer.rd_adr(), distance, (int) fence_queue.size());
                }


            }

            //----------------------------------------------------------------------------------------------------------
            // Bx Counter
            //----------------------------------------------------------------------------------------------------------

            // request processing on the next pretrig popping out of l1a delay
            if (pretrig_queue.size() > 0)
                next_bx_rqst = std::min(next_bx_rqst, pretrig_queue[0] + l1a_latency);

            if (rd_started) { // # of l1as for readout
                next_bx_rqst = std::min(next_bx_rqst, rd_finish_bx);
            }
            else if (fence_queue.size() > 0) { // # of l1as for readout
                next_bx_rqst = bx+1;
            }

            next_bx_rqst = std::min(next_bx_rqst, bx+distance);

            next_bx_rqst = std::min(next_bx_rqst, pretrig_bx);
            //printf ("request @ bx=%i (current_bx=%i)\n", next_bx_rqst, bx);

            if (no_bx_skipping)
                next_bx_rqst=-1;

            if (next_bx_rqst!=unsigned long long(-1)) {
                delta_bx = next_bx_rqst - bx;
                bx = next_bx_rqst;
                if (delta_bx==0) {
                    delta_bx=1;
                    bx++;
                }
            }
            else {
                delta_bx = 1;
                bx++;
            }



            if (rd_started==1) {
                ring_buffer.inc_rd_adr(delta_bx);
            }

            if (!buf_stalled) {
                // we get a readout signal l1a_delay later than the event pretrigger time
                if (debug_core) printf("incrementing wr_buf by %llu", delta_bx);
                if (debug_core) printf("   from %i", ring_buffer.wr_adr());
                ring_buffer.inc_wr_adr(delta_bx);
                if (debug_core) printf("   to %i\n", ring_buffer.wr_adr());
            }


            seconds = bx / 40e6;

            if (debug_core && bx>2322)
                return 0;

            //----------------------------------------------------------------------------------------------------------
            // Printout & Fill Histograms
            //----------------------------------------------------------------------------------------------------------

            //  this only needs to be checked:
            //    (1) we had a preclct on this bx

            if (!debug_core) {

            int min_pretrigs = lum < 150 ? 1000000 : 200000;

            if (seconds>seconds_cutoff || (preclct_cnt >= min_pretrigs/(speedup) && loss_hunt_threshold<=0 && preclct_lost >= min_lost)) {

                float preclct_rate = (float(preclct_cnt) * pow(10,9) / (25. * bx)); // Hz
                float readout_rate = (l1a_match_cnt      * pow(10,9) / (25. * bx)); // Hz
                float stall_rate   = (num_stalls         * pow(10,9) / (25. * bx)); // Hz

                float lost_event_frac = float (preclct_lost)/float(preclct_cnt); // fraction

                //float luminosity = preclct_rate * l1a_match_fraction / lumi_rate / 10000.;
                float luminosity = lum/10.;
                printf("preclct_rate=%7.1f kHz, readout_rate=%5.1f kHz, preclct_cnt=%7i, pretrigs_dropped=%5i, mean_occupancy=%4.2f, luminosity=%4.1f *10^34, lhc_seconds=%f\n",
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
                    loss_hists.get2D("h2_loss_me11_lowl1") -> Fill(lum/10., lost_event_frac);
                }
                // deep buffers
                else if (deep_buffer==1 && !low_l1 && istation==0 && !gem_en) {
                    if (ddr_readout)
                        loss_hists.get2D("h2_loss_me11_deep_ddr") -> Fill(lum/10., lost_event_frac);
                    else
                        loss_hists.get2D("h2_loss_me11_deep") -> Fill(lum/10., lost_event_frac);
                }
                // ddr readout
                else if (ddr_readout && !deep_buffer && !low_l1 && !gem_en && !unfurled_triads) {
                        (lostevents_ddr_arr[istation])  -> Fill(lum/10., lost_event_frac);
                }
                // gem
                else if (gem_en && !low_l1 && !ddr_readout) {
                    if (ddr_readout)
                        loss_hists.get2D("h2_loss_me11_gem_ddr") -> Fill(lum/10., lost_event_frac);
                    else
                        loss_hists.get2D("h2_loss_me11_gem") -> Fill(lum/10., lost_event_frac);
                }
                // unfurled triads
                else if (unfurled_triads && ddr_readout && !gem_en && !low_l1) {
                        loss_hists.get2D("h2_loss_me11_unfurled_ddr") -> Fill(lum/10., lost_event_frac);
                }
                // baselines
                else if (!gem_en && !low_l1 && !gem_en && !ddr_readout) {

                        (lostevents_arr[istation])  -> Fill(lum/10., lost_event_frac);

                        if (istation==0) {
                            misc_hists.get2D("h2_buf_stall") -> Fill(lum/10., stall_rate);
                            misc_hists.get2D("h2_queue_occupancy") -> Fill (lum/10., sum_queue_occupancy/double(l1a_match_cnt));
                        }
                }


                break;
            }
            }
            else if (loss_hunt_threshold>0) {
                if (preclct_lost > loss_hunt_threshold) {
                    printf("\nDONE\n");
                    return 0;
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

    if (loss_hunt_threshold<=0) {
        gSystem->Exec("hadd -f tmb_model.root tmb_model_loss.root tmb_model_misc.root");
    }

    gApplication->Terminate();
    return 0;
}
