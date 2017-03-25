#ifndef FRAME_COUNT_H
#define FRAME_COUNT_H 

int frame_count () {

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

    return dmb_frame_cnt;
}

#endif /* FRAME_COUNT_H */
