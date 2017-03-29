#ifndef FENCE_DISTANCE_H
#define FENCE_DISTANCE_H

// calculate the distance in a fence queue
int fence_distance (int buf_queue_adr, int fifo_wadr, int buffer_size) {

    int fifo_pretrig_cfeb    = 7;
    int fifo_pretrig_gem_rpc = 7;

    int READ_ADR_OFFSET = 6;
    int PRESTORE_SAFETY = 2;

    int prestore_setback =  READ_ADR_OFFSET + 1 + PRESTORE_SAFETY;

    // take maximum pretrigger setback (GEM or RPC or CSC)
    int pretrig_setback  = (fifo_pretrig_cfeb >= fifo_pretrig_gem_rpc)
        ?  fifo_pretrig_cfeb : fifo_pretrig_gem_rpc;

    // total buffer setback is the sum of (pretrig + prestore)
    int buf_setback = pretrig_setback-prestore_setback;

    int next_fence_adr = buf_queue_adr-buf_setback;  // compensate for pre-trig latency

    int distance = (buffer_size-1) & ((next_fence_adr - fifo_wadr) % buffer_size);

    return distance;
}


#endif /* FENCE_DISTANCE_H */
