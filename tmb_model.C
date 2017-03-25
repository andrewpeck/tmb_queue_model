#include <boost/circular_buffer.hpp>
#include <cmath>

#include "fence_distance.h"
#include "frame_count.h"
#include "ring_buffer.h"
#include <vector>

int    mean_rate       = 150000; // Hz

ringBuffer ring_buffer(2048);

TRandom rando;

bool debug=0;


void tmb_model () {


    double dmb_readout_rate = 1000000.0/(frame_count()*25); // DMB readout time, kHz

    printf("Frame count          : %i \n", frame_count());
    printf("Maximumn readout rate: %f kHz\n", dmb_readout_rate);


    for (mean_rate = 60000; mean_rate < dmb_readout_rate*1000+10000; mean_rate+=1000) {

        int    mean_separation = (int) ( (1./mean_rate) * ( pow(10,9) ) / (25.)); // bunches

        std::vector<int> fence_queue;

        int l1a_latency = 128;


        long long unsigned bx = 0;
        long long unsigned pretrig_bx = 256;

        bool rd_started  = 0;
        bool rd_finished = 1;
        long long unsigned rd_finish_bx = -1;

        int distance = 2048;

        int num_overflows = 0;
        int num_events = 0;
        int num_flushed   = 0;
        int num_stalls    = 0;
        bool flush_on_overflow = 0;
        bool flush = 0;

        while (true) {

            //----------------------------------------------------------------------------------------------------------
            // L1A Queue Inserter
            //----------------------------------------------------------------------------------------------------------

            if (bx==pretrig_bx) {

                // calculate bx time for next readout trigger
                rando.SetSeed(0);
                int separation = -log(1-rando.Uniform(1)) / (mean_rate*25 / pow(10,9));

                if (separation < 6) {
                    // reroll
                    pretrig_bx = bx + 1;
                    continue;
                }

                //printf("separation = %i\n", separation);
                pretrig_bx = bx + separation;

                //if (debug) printf("prematched readout in bx=%lld, delta=%lld\n",pretrig_bx, separation);

                // push pretrig_bx to fence queue for raw hits readout to DMB
                // fence_queue.push_back (pretrig_bx-l1a_latency);
                int push_back_adr = 0x7ff & (ring_buffer.wr_adr()-l1a_latency);
                fence_queue.push_back (push_back_adr);

                if (debug) printf("event pushed     @ %4i, events in queue=%i\n",push_back_adr, (int) fence_queue.size());

                if (distance==0) {
                    if (debug) printf("Event lost (buf stalled)!\n");
                    if (debug) printf("bx=%lld\n",bx);
                    num_overflows ++;
                }

                if (flush) {
                    if (debug) printf("Event lost (flushing)!\n");
                    if (debug) printf("bx=%lld\n",bx);
                    num_flushed ++;
                }

                num_events ++;

            }

            //----------------------------------------------------------------------------------------------------------
            // Fence Distance
            //----------------------------------------------------------------------------------------------------------

            if (fence_queue.size() > 0) {
                // make sure that the fence buffer has space--- otherwise we should stall here
                distance = fence_distance(fence_queue[0], ring_buffer.wr_adr());
                //if (debug) printf("bx=%llu, distance=%i, fence_adr=%i, wr_adr=%i, queue_size=%i\n",bx, distance, fence_queue[0], ring_buffer.wr_adr(), (int) fence_queue.size());
            }
            else {
                flush = 0;
                distance = l1a_latency;
                //if (debug) printf("bx=%llu queue empty\n", bx);
            }

            //----------------------------------------------------------------------------------------------------------
            // Buffer write control
            //----------------------------------------------------------------------------------------------------------

            if (flush)
            {
                //                if (debug) printf("Flushing!\n");
            }
            else if (distance==0) {

                if (flush_on_overflow) {
                     flush = 1;
                }

                num_stalls++;
            }
            else {
                // we get a readout signal l1a_delay later than the event pretrigger time
                ring_buffer.inc_wr_adr();
            }


            //----------------------------------------------------------------------------------------------------------
            // Buffer read control
            //----------------------------------------------------------------------------------------------------------

            if (fence_queue.size() > 0) { // # of l1as for readout

                // we have events in queue but they haven't been readout yet
                if (rd_started==0) {

                    // mark readout as busy
                    rd_started=1;

                    ring_buffer.rd_adr(fence_queue[0]);
                    // mark the last busy_bx
                    int cnt = frame_count ();
                    rd_finish_bx = bx + cnt;

                    if (debug) printf("readout started  @ %4llu, events in queue=%i (cnt=%i)\n", fence_queue[0], (int) fence_queue.size(), cnt);

                }

                else if (bx==rd_finish_bx) {
                    rd_started=0;

                    // pop from queue after readout
                    fence_queue.erase(fence_queue.begin());

                    if (debug) printf("readout finished @ %4llu, events in queue=%i\n", ring_buffer.rd_adr(), (int) fence_queue.size());
                }
                else {
                    ring_buffer.inc_rd_adr();
                }

            }

            // increment bx
            bx++;

            if (bx>4000000) {
                printf("rate=%6i, num_overflows=%6i, num_stalls=%9i, num_events=%8i, rate=%f\n",mean_rate, num_overflows, num_stalls, num_events, float(num_overflows)/num_events * 101.0) ;
                break;
            }

        }
    }

}


// void l1a_latency () {
//
//     int dmb_frame_cnt = frame_count();
//
//     int buffer_size = 2048;
//     int minimum_fence = 64;
//
//     int maximum_occupancy = (buffer_size/dmb_frame_cnt);
//
//     // we just convert number of bx to kHz
//     double dmb_readout_rate = 1000000.0/(dmb_frame_cnt*25); // DMB readout time, kHz
//     double service_time     = 1./dmb_readout_rate;
//
//
//
//     printf ("Configured with:\n");
//     printf ("    n_cfebs              = %i\n" , n_cfebs);
//     printf ("    n_tbins              = %i\n" , n_tbins);
//     printf ("    rpc_readout          = %i\n" , rpc_en);
//     printf ("    scope_enabled        = %i\n" , scope_en);
//     printf ("    miniscope_enabled    = %i\n" , miniscope_en);
//     printf ("    blocked_cfeb_readout = %i\n" , blocked_cfeb_readout);
//     printf ("\n");
//     printf ("    Word Count           = %i\n" , dmb_frame_cnt);
//     printf ("    DMB Readout Rate     = %f kHz (l1a/sec)\n" , dmb_readout_rate);
//
//
//     TCanvas *c1 = new TCanvas("c1");
//
//     TF1 *ovf = new TF1 (
//             "ovf",
//             "blocking_probability([0], x*[1])",
//             0,
//             dmb_readout_rate*2
//             );
//
//     ovf->SetParameter(0,maximum_occupancy);
//     ovf->SetParameter(1,service_time);
//     ovf->SetNpx(10000);
//     ovf->GetHistogram()->GetXaxis()->SetTitle("l1a*lct rate (kHz)");
//     ovf->Draw();
//     c1->SetLogy();
//     ovf->SetTitle("overflow probability");
//
//
//     //////////////////////////////////////////////////////////////
//
//     TCanvas *c2 = new TCanvas("c2");
//
//     TF1 *occ = new TF1 (
//             "occ",
//             "mean_queue_occupancy(x, [0])",
//             0,
//             dmb_readout_rate*9/10
//             );
//
//     occ->SetParameter(0,dmb_readout_rate*9/10);
//     occ->SetNpx(10000);
//     occ->Draw();
//     c2->SetLogy();
//     occ->GetHistogram()->GetXaxis()->SetTitle("l1a*lct rate (kHz)");
//     occ->GetHistogram()->GetYaxis()->SetTitle("events");
//     occ->SetTitle("mean l1as in readout queue");
//
//
//
//
//     // mean_queue_occupancy = mean_queue_occupancy (l1a_rate, dmb_readout_rate)
//
//     // plt.xlabel("l1a_rate (kHz)")
//     // plt.ylabel("overflow probability")
//     // plt.semilogy(l1a_rate,overflow_probability)
//
//     //printf ("    DMB Readout Rate = %f kHz\n"            ,  dmb_readout_rate);
//     //printf ("    Mean queue occupancy = %f (%0.2f\%)\n"  ,  mean_occupancy, (mean_occupancy/maximum_occupancy*100));
//     //printf ("    Mean latency = %f bx\n"                 ,  mean_occupancy * dmb_frame_cnt);
//     //printf ("    Blocking Probability = %f\n"            ,  blocking_probability(maximum_occupancy,rho));
//
//     // l1a_rate             = np.linspace(0,dmb_readout_rate/10-1,4068,endpoint=True)
//     // overflow_probability = blocking_probability(maximum_occupancy,(l1a_rate/dmb_readout_rate))
//     // mean_queue_occupancy = mean_queue_occupancy (l1a_rate, dmb_readout_rate)
//     //
//     // #plt.semilogy(l1a_rate,overflow_probability)
//     // plt.grid(True)
//     // plt.show()
//     //
//     // plt.xlabel("l1a_rate (kHz)")
//     // plt.ylabel("mean queue occupancy (# of L1As)")
//     // plt.semilogy(l1a_rate,mean_queue_occupancy)
//     // #plt.plot(l1a_rate,mean_queue_occupancy)
//     // plt.grid(True)
//     // plt.show()
// }
