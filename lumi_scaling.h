#ifndef LUMI_SCALING_H
#define LUMI_SCALING_H 

float lumi_rate_pre_me11 = (30.3+30.6+39.6+29.1)/4;
float lumi_rate_pre_me21 = (14.6+14.6+13.6+14.4)/4;
float lumi_rate_pre_me31 = (8.10+8.11+6.00+8.00)/4;
float lumi_rate_pre_me41 = (9.25+9.28+8.40+9.13)/4;

float lumi_rate_me11 = (0.612+0.599+0.506+0.545)/4;
float lumi_rate_me21 = (0.347+0.341+0.296+0.297)/4;
float lumi_rate_me31 = (0.211+0.197+0.196+0.192)/4;
float lumi_rate_me41 = (0.209+0.207+0.189+0.179)/4;

float lumi_rate_pre_arr [4] = { lumi_rate_pre_me11, lumi_rate_pre_me21, lumi_rate_pre_me31, lumi_rate_pre_me41};
float lumi_rate_arr     [4] = { lumi_rate_me11, lumi_rate_me21, lumi_rate_me31, lumi_rate_me41};

#endif /* LUMI_SCALING_H */
