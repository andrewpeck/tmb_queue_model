#ifndef FRAME_COUNT_H
#define FRAME_COUNT_H

#include <cstdio>

int frame_count (int n_cfebs=7, int gem_en=0, bool unfurled_triads=0, int n_tbins=12) {

    bool debug=0;

    int n_headers    = 42;
    //int n_cfebs      = 7;
    //int n_tbins      = 12;
    int n_gems       = 2;
    int n_rpcs       = 2;
    int n_scope_chan = 128;

    bool rpc_en               = false;
    //bool gem_en               = true;
    bool scope_en             = false;
    bool miniscope_en         = false;
    bool blocked_cfeb_readout = false;

    // header
    int wc_eob        = 1;

    // cfeb
    int wc_cfeb       = n_cfebs * 6 * n_tbins  * (unfurled_triads ? 4 : 1);


    // rpc + gem
    int wc_b04        = (rpc_en || gem_en) ? 1 : 0;
    int wc_gems       = gem_en ? n_gems  * 8 * n_tbins : 0;
    int wc_rpcs       = rpc_en ? 2*n_rpcs*n_tbins : 0;
    int wc_e04        = (rpc_en || gem_en) ? 1 : 0;

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

    dmb_frame_cnt = n_headers + wc_eob + wc_cfeb + wc_b04 + wc_gems + wc_rpcs + wc_e04 + wc_b05 + wc_scope + wc_e05 + wc_b07 + wc_miniscope + wc_bcb + wc_b_cfeb + wc_ecb + wc_eoc + wc_multiple + wc_eof + wc_crc + wc_wc;

    int wc_mod = ((dmb_frame_cnt+3) & ~(0x03) - dmb_frame_cnt);; //  round up to nearest multiple of 4

    dmb_frame_cnt = (dmb_frame_cnt+3) & ~(0x03); //  round up to nearest multiple of 4

    if (debug) printf("n_headers    = %i\n", n_headers);
    if (debug) printf("wc_eob       = %i\n", wc_eob);
    if (debug) printf("wc_cfeb      = %i\n", wc_cfeb);
    if (debug) printf("wc_b04       = %i\n", wc_b04);
    if (debug) printf("wc_gems      = %i\n", wc_gems);
    if (debug) printf("wc_rpcs      = %i\n", wc_rpcs);
    if (debug) printf("wc_e04       = %i\n", wc_e04);
    if (debug) printf("wc_b05       = %i\n", wc_b05);
    if (debug) printf("wc_scope     = %i\n", wc_scope);
    if (debug) printf("wc_e05       = %i\n", wc_e05);
    if (debug) printf("wc_b07       = %i\n", wc_b07);
    if (debug) printf("wc_miniscope = %i\n", wc_miniscope);
    if (debug) printf("wc_bcb       = %i\n", wc_bcb);
    if (debug) printf("wc_b_cfeb    = %i\n", wc_b_cfeb);
    if (debug) printf("wc_ecb       = %i\n", wc_ecb);
    if (debug) printf("wc_eoc       = %i\n", wc_eoc);
    if (debug) printf("wc_multiple  = %i\n", wc_multiple);
    if (debug) printf("wc_eof       = %i\n", wc_eof);
    if (debug) printf("wc_crc       = %i\n", wc_crc);
    if (debug) printf("wc_wc        = %i\n", wc_wc);
    if (debug) printf("wc_mod       = %i\n", wc_mod);

    if (debug) printf("word_count = %i\n", dmb_frame_cnt);

    return dmb_frame_cnt;
}

#endif /* FRAME_COUNT_H */
