#!/usr/bin/python3

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams['backend'] = "Qt5Agg"

n_headers = 42
n_cfebs = 7
n_rpcs = 2
n_tbins = 4
n_chan = 128

rpc_en                =  False
scope_en              =  False
miniscope_en          =  True
blocked_cfeb_readout  =  False

wc_eob        = 1
wc_cfeb       = n_cfebs * 6 * n_tbins
wc_b04        = 1 if (rpc_en) else 0
wc_rpcs       = 2*n_rpcs*n_tbins if (rpc_en) else 0
wc_e04        = 1 if (rpc_en) else 0

#scope
wc_b05        = 1             if (scope_en) else 0
wc_scope      = n_chan/16*256 if (scope_en) else 0
wc_e05        = 1             if (scope_en) else 0

#miniscope
wc_b07        = 1  if (miniscope_en) else 0
wc_miniscope  = 22 if (miniscope_en) else 0

#blocked cfeb readout
wc_bcb        = 1  if (blocked_cfeb_readout) else 0
wc_b_cfeb     = 22 if (blocked_cfeb_readout) else 0
wc_ecb        = 1  if (blocked_cfeb_readout) else 0

wc_eoc        = 1
wc_multiple   = 0 # we account for this later
wc_eof        = 1
wc_crc        = 2
wc_wc         = 1

dmb_frame_cnt = n_headers + wc_eob + wc_cfeb + wc_b04 + wc_rpcs + wc_e04 + wc_b05 + wc_scope + wc_e05 + wc_b07 + wc_miniscope + wc_bcb + wc_b_cfeb + wc_ecb + wc_eoc + wc_multiple + wc_eof + wc_crc + wc_wc 
dmb_frame_cnt = (dmb_frame_cnt+3) & ~(0x03) # round up to nearest multiple of 4 

print ("Configured with:")
print ("    rpc_readout: %s" % rpc_en)
print ("    scope_enabled: %s" % scope_en)
print ("    miniscope_enabled: %s" % miniscope_en)
print ("    blocked_cfeb_readout: %s" % blocked_cfeb_readout)
print ("")
print ("    Word Count = %d" % dmb_frame_cnt)


buffer_size = 2048
minimum_fence = 64

maximum_occupancy = math.floor (buffer_size/dmb_frame_cnt)

# we just convert number of bx to kHz
dmb_readout_rate = 1000000.0/(dmb_frame_cnt*25) # DMB readout time, kHz
print ("    DMB Readout Rate = %d kHz (l1a/sec)" % dmb_readout_rate)

def mean_queue_occupancy (arrival, service): 
    rho = arrival/float(service)
    n = (rho/(1.0-rho)) * 0.5*(2.0-rho)
    return n

# calculate the blocking probability for a finite queue length with deterministic, fixed readout time (M/D/1 queue)
# cf. http://dx.doi.org/10.4218/etrij.14.0113.0812
def blocking_probability (queue_size, traffic_intensity): 
    k = queue_size 
    rho = traffic_intensity 

    sigma = 0
    for j in range (0,k): 
        sigma = sigma + ((-1)**j * rho**j * (k-j)**j * np.exp(rho*(k-j)))/math.factorial(j)
    e_k = 1-(1-rho)*sigma

    pb = (1-rho)*e_k/(1-rho*e_k)

    return pb

#print ("    DMB Readout Rate = %f kHz" % service_rate)
#print ("    Mean queue occupancy = %f (%0.2f%%)" % (mean_occupancy, (mean_occupancy/maximum_occupancy*100)))
#print ("    Mean latency = %f bx" % (mean_occupancy * dmb_frame_cnt))
#print ("    Blocking Probability = %f" % (blocking_probability(maximum_occupancy,rho))) 

l1a_rate             = np.linspace(0,dmb_readout_rate+1,4068,endpoint=True)
overflow_probability = blocking_probability(maximum_occupancy,(l1a_rate/dmb_readout_rate))
mean_queue_occupancy = mean_queue_occupancy (l1a_rate, dmb_readout_rate)

#plt.semilogy(l1a_rate,overflow_probability)
plt.xlabel("l1a_rate (kHz)")
plt.ylabel("overflow probability")
plt.plot(l1a_rate,overflow_probability)
#plt.semilogy(l1a_rate,overflow_probability)
plt.grid(True)
plt.show()

plt.xlabel("l1a_rate (kHz)")
plt.ylabel("mean queue occupancy (# of L1As)")
plt.semilogy(l1a_rate,mean_queue_occupancy)
#plt.plot(l1a_rate,mean_queue_occupancy)
plt.grid(True)
plt.show()
