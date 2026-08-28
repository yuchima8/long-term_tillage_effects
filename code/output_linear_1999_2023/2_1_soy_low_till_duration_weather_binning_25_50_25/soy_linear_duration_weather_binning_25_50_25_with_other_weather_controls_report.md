# 20-60-20 Binning Estimator for Conditional Low-Till Duration Effects on soybean yield (With Other Weather Controls)

## Design
- Outcome: log soybean yield in bushels per acre (`log_y`).
- Treatment: continuous low-till duration in years (`low_till_duration`).
- Control/reference duration: clean conventional-till observations with no prior low-till exposure (`till == 0` and `till_1_count == 0`).
- Low-till observations: current low-till in an uninterrupted low-till spell (`till == 1` and `till_1_count == till_1_streak`).
- Excluded years: none..
- Focal weather moderator: one of `GDD_4_5`, `PPT_4_5`, `GDD_6_9`, `PPT_6_9`, or `EDD_6_9`, depending on the panel.
- Additive weather controls: the other four weather variables enter linearly outside the bins in each panel.
- Full weather sample: no focal-weather trimming is applied before constructing 20-60-20 weather bins and fitting the estimator.
- Pooled estimator: each focal weather moderator is split into low/middle/high bins using the bottom 20%, middle 60%, and top 20% of the full focal-weather sample.
- Model specification: control-group bin intercept dummies, bin-specific low-till duration effects, bin-specific focal-weather slopes, bin-specific duration-by-focal-weather interactions, and additive controls for the other weather variables.
- Fixed effects: field and year.
- Standard errors clustered by `FIPS`.
- Weights: field size.
- Visualization note: the bottom histogram band shows only the middle 95% of the focal-weather distribution (2.5th-97.5th percentiles) to keep extreme tails from dominating the plot scale.

## Interpretation
- `estimate_log_points_per_year` is the conditional marginal effect of one additional year of continuous low-till duration at the median focal-weather value within that bin.
- `pct_effect_per_year` converts that log-point estimate to percent using `100 * (exp(estimate) - 1)`.
- `additive_control_terms` reports the linear control coefficients for the non-focal weather variables included outside the bins.

## Sample Counts
                                                sample observations  fields counties years low_till_observations mean_low_till_duration
                                                <char>        <int>   <int>    <int> <int>                 <int>                  <num>
1:                                  full_cleaned_panel     44912332 4158139      796    25                    NA                     NA
2: clean_controls_and_continuous_low_till_with_weather     23399573 3550567      791    25              10226005               8.486174
3:                                   estimation_sample     21156753 3202625      790    25               9485464               8.586234

## Conditional Low-Till Duration Effects by Weather Bin
    weather_var weather_bin weather_median observations estimate_log_points_per_year    std_error pct_effect_per_year pct_effect_low_per_year
         <char>      <char>          <num>        <int>                        <num>        <num>               <num>                   <num>
 1:     GDD_4_5         low     252.638739      5290231                 0.0012067493 0.0003163408          0.12074778             0.058689358
 2:     GDD_4_5      middle     351.737962     10577904                 0.0003977048 0.0002175306          0.03977839            -0.002865482
 3:     GDD_4_5        high     470.924409      5288618                 0.0030508338 0.0002488159          0.30554923             0.256644240
 4:     PPT_4_5         low     121.725662      5289177                 0.0017356199 0.0002798335          0.17371269             0.118785117
 5:     PPT_4_5      middle     194.432367     10576420                 0.0014036786 0.0002374156          0.14046642             0.093878437
 6:     PPT_4_5        high     280.147198      5291156                 0.0011387664 0.0002841697          0.11394151             0.058196305
 7:     GDD_6_9         low    1334.096270      5288215                -0.0001917463 0.0003999759         -0.01917279            -0.097522327
 8:     GDD_6_9      middle    1520.316629     10578844                 0.0002277854 0.0002062425          0.02278114            -0.017643427
 9:     GDD_6_9        high    1743.247590      5289694                 0.0035673302 0.0002679751          0.35737007             0.304673092
10:     PPT_6_9         low     261.118323      5288115                 0.0023745837 0.0002917713          0.23774052             0.180433782
11:     PPT_6_9      middle     365.032176     10577902                 0.0016719530 0.0002188697          0.16733515             0.124374114
12:     PPT_6_9        high     508.925419      5290736                 0.0007078746 0.0002689982          0.07081252             0.018065452
13:     EDD_6_9         low       2.319834      5289136                 0.0005430643 0.0002978347          0.05432118            -0.004069086
14:     EDD_6_9      middle      11.175295     10579104                 0.0006924669 0.0002225535          0.06927067             0.025629498
15:     EDD_6_9        high      31.033865      5288513                 0.0036405622 0.0002692783          0.36471971             0.311762637
    pct_effect_high_per_year
                       <num>
 1:               0.18284468
 2:               0.08244044
 3:               0.35447808
 4:               0.22867040
 5:               0.18707609
 6:               0.16971777
 7:               0.05923819
 8:               0.06322204
 9:               0.41009474
10:               0.29508004
11:               0.21031463
12:               0.12358741
13:               0.11274554
14:               0.11293089
15:               0.41770473

## Additive Other-Weather Control Terms
    weather_var focal_weather_source_var additive_control      estimate    std_error        ci_low       ci_high
         <char>                   <char>           <char>         <num>        <num>         <num>         <num>
 1:     GDD_4_5                  GDD_4_5          ppt_4_5 -1.278171e-04 1.553905e-05 -1.582736e-04 -9.736055e-05
 2:     GDD_4_5                  GDD_4_5          GDD_6_9  6.347902e-06 2.788934e-05 -4.831521e-05  6.101102e-05
 3:     GDD_4_5                  GDD_4_5          ppt_6_9  7.667903e-05 1.015560e-05  5.677406e-05  9.658400e-05
 4:     GDD_4_5                  GDD_4_5          EDD_6_9 -6.110452e-03 1.465272e-04 -6.397645e-03 -5.823259e-03
 5:     PPT_4_5                  ppt_4_5          GDD_4_5  6.005941e-04 3.485363e-05  5.322810e-04  6.689072e-04
 6:     PPT_4_5                  ppt_4_5          GDD_6_9 -5.501186e-06 2.776059e-05 -5.991194e-05  4.890956e-05
 7:     PPT_4_5                  ppt_4_5          ppt_6_9  7.779005e-05 1.038544e-05  5.743460e-05  9.814551e-05
 8:     PPT_4_5                  ppt_4_5          EDD_6_9 -5.832345e-03 1.478193e-04 -6.122071e-03 -5.542619e-03
 9:     GDD_6_9                  GDD_6_9          GDD_4_5  5.053079e-04 3.510011e-05  4.365117e-04  5.741041e-04
10:     GDD_6_9                  GDD_6_9          ppt_4_5 -1.200161e-04 1.374854e-05 -1.469632e-04 -9.306896e-05
11:     GDD_6_9                  GDD_6_9          ppt_6_9  7.363268e-05 1.004126e-05  5.395180e-05  9.331356e-05
12:     GDD_6_9                  GDD_6_9          EDD_6_9 -5.367606e-03 1.526988e-04 -5.666896e-03 -5.068317e-03
13:     PPT_6_9                  ppt_6_9          GDD_4_5  5.512903e-04 3.512841e-05  4.824386e-04  6.201420e-04
14:     PPT_6_9                  ppt_6_9          ppt_4_5 -1.160571e-04 1.433973e-05 -1.441630e-04 -8.795125e-05
15:     PPT_6_9                  ppt_6_9          GDD_6_9  7.037738e-05 2.947582e-05  1.260477e-05  1.281500e-04
16:     PPT_6_9                  ppt_6_9          EDD_6_9 -5.877720e-03 1.352078e-04 -6.142727e-03 -5.612712e-03
17:     EDD_6_9                  EDD_6_9          GDD_4_5  5.877259e-04 3.780245e-05  5.136331e-04  6.618187e-04
18:     EDD_6_9                  EDD_6_9          ppt_4_5 -1.187545e-04 1.398081e-05 -1.461569e-04 -9.135215e-05
19:     EDD_6_9                  EDD_6_9          GDD_6_9  6.141127e-05 2.483269e-05  1.273919e-05  1.100833e-04
20:     EDD_6_9                  EDD_6_9          ppt_6_9  7.336816e-05 1.001524e-05  5.373829e-05  9.299803e-05
    weather_var focal_weather_source_var additive_control      estimate    std_error        ci_low       ci_high

## Model Summaries
### GDD_4_5
OLS estimation, Dep. Var.: log_y
Observations: 20,753,599
Weights: size
Fixed-effects: OBJECTID: 2,799,471,  year: 25
Standard-errors: Clustered (FIPS) 
                                     Estimate Std. Error    t value   Pr(>|t|)    
bin_middle                         0.06215715 0.00411929  15.089290  < 2.2e-16 ***
bin_high                           0.12092009 0.00793572  15.237449  < 2.2e-16 ***
duration_years_low                 0.00120675 0.00031634   3.814713 1.4715e-04 ***
duration_years_middle              0.00039770 0.00021753   1.828270 6.7892e-02 .  
duration_years_high                0.00305083 0.00024882  12.261412  < 2.2e-16 ***
weather_c_low                      0.00041062 0.00006515   6.302745 4.8945e-10 ***
weather_c_middle                   0.00054916 0.00005299  10.362907  < 2.2e-16 ***
weather_c_high                     0.00041904 0.00004690   8.933799  < 2.2e-16 ***
duration_years_x_weather_c_low     0.00000911 0.00000740   1.231889 2.1836e-01    
duration_years_x_weather_c_middle  0.00000569 0.00000388   1.463732 1.4367e-01    
duration_years_x_weather_c_high    0.00004157 0.00000350  11.875070  < 2.2e-16 ***
ppt_4_5                           -0.00012782 0.00001554  -8.225542 8.1387e-16 ***
GDD_6_9                            0.00000635 0.00002789   0.227610 8.2001e-01    
ppt_6_9                            0.00007668 0.00001016   7.550422 1.2165e-13 ***
EDD_6_9                           -0.00611045 0.00014653 -41.701834  < 2.2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
RMSE: 1.30874     Adj. R2: 0.683981
                Within R2: 0.139478

### PPT_4_5
OLS estimation, Dep. Var.: log_y
Observations: 20,753,599
Weights: size
Fixed-effects: OBJECTID: 2,799,471,  year: 25
Standard-errors: Clustered (FIPS) 
                                     Estimate Std. Error    t value   Pr(>|t|)    
bin_middle                         0.00049996 0.00270294   0.184967 8.5330e-01    
bin_high                          -0.01465014 0.00345209  -4.243851 2.4622e-05 ***
duration_years_low                 0.00173562 0.00027983   6.202331 9.0293e-10 ***
duration_years_middle              0.00140368 0.00023742   5.912326 5.0471e-09 ***
duration_years_high                0.00113877 0.00028417   4.007346 6.7305e-05 ***
weather_c_low                      0.00030142 0.00006352   4.745645 2.4733e-06 ***
weather_c_middle                  -0.00012947 0.00006134  -2.110856 3.5103e-02 *  
weather_c_high                    -0.00021375 0.00004893  -4.368275 1.4219e-05 ***
duration_years_x_weather_c_low    -0.00002649 0.00000609  -4.350634 1.5383e-05 ***
duration_years_x_weather_c_middle -0.00000199 0.00000531  -0.374698 7.0799e-01    
duration_years_x_weather_c_high    0.00000819 0.00000441   1.857770 6.3579e-02 .  
GDD_4_5                            0.00060059 0.00003485  17.231896  < 2.2e-16 ***
GDD_6_9                           -0.00000550 0.00002776  -0.198165 8.4297e-01    
ppt_6_9                            0.00007779 0.00001039   7.490301 1.8677e-13 ***
EDD_6_9                           -0.00583235 0.00014782 -39.455901  < 2.2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
RMSE: 1.31001     Adj. R2: 0.683367
                Within R2: 0.137806

### GDD_6_9
OLS estimation, Dep. Var.: log_y
Observations: 20,753,599
Weights: size
Fixed-effects: OBJECTID: 2,799,471,  year: 25
Standard-errors: Clustered (FIPS) 
                                     Estimate Std. Error    t value   Pr(>|t|)    
bin_middle                         0.00029899 0.00567526   0.052684 9.5800e-01    
bin_high                          -0.04046810 0.01178165  -3.434841 6.2431e-04 ***
duration_years_low                -0.00019175 0.00039998  -0.479395 6.3179e-01    
duration_years_middle              0.00022779 0.00020624   1.104454 2.6974e-01    
duration_years_high                0.00356733 0.00026798  13.312171  < 2.2e-16 ***
weather_c_low                      0.00016806 0.00003411   4.927549 1.0178e-06 ***
weather_c_middle                  -0.00014317 0.00003611  -3.964987 8.0165e-05 ***
weather_c_high                    -0.00004091 0.00004488  -0.911521 3.6230e-01    
duration_years_x_weather_c_low    -0.00000906 0.00000446  -2.031476 4.2546e-02 *  
duration_years_x_weather_c_middle  0.00000887 0.00000189   4.695948 3.1364e-06 ***
duration_years_x_weather_c_high    0.00001323 0.00000247   5.349249 1.1623e-07 ***
GDD_4_5                            0.00050531 0.00003510  14.396190  < 2.2e-16 ***
ppt_4_5                           -0.00012002 0.00001375  -8.729373  < 2.2e-16 ***
ppt_6_9                            0.00007363 0.00001004   7.333008 5.6559e-13 ***
EDD_6_9                           -0.00536761 0.00015270 -35.151603  < 2.2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
RMSE: 1.30497     Adj. R2: 0.685802
                Within R2: 0.144436

### PPT_6_9
OLS estimation, Dep. Var.: log_y
Observations: 20,753,599
Weights: size
Fixed-effects: OBJECTID: 2,799,471,  year: 25
Standard-errors: Clustered (FIPS) 
                                     Estimate Std. Error    t value   Pr(>|t|)    
bin_middle                         0.02562178 0.00290040   8.833867  < 2.2e-16 ***
bin_high                           0.04016034 0.00378146  10.620318  < 2.2e-16 ***
duration_years_low                 0.00237458 0.00029177   8.138511 1.5828e-15 ***
duration_years_middle              0.00167195 0.00021887   7.639032 6.4339e-14 ***
duration_years_high                0.00070787 0.00026900   2.631522 8.6687e-03 ** 
weather_c_low                      0.00043993 0.00007618   5.774874 1.1129e-08 ***
weather_c_middle                   0.00018972 0.00003704   5.122375 3.8077e-07 ***
weather_c_high                    -0.00016622 0.00002382  -6.978703 6.3823e-12 ***
duration_years_x_weather_c_low     0.00000261 0.00000601   0.433890 6.6449e-01    
duration_years_x_weather_c_middle -0.00001066 0.00000317  -3.363679 8.0675e-04 ***
duration_years_x_weather_c_high    0.00000538 0.00000212   2.539635 1.1290e-02 *  
GDD_4_5                            0.00055129 0.00003513  15.693576  < 2.2e-16 ***
ppt_4_5                           -0.00011606 0.00001434  -8.093398 2.2294e-15 ***
GDD_6_9                            0.00007038 0.00002948   2.387631 1.7195e-02 *  
EDD_6_9                           -0.00587772 0.00013521 -43.471749  < 2.2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
RMSE: 1.30341     Adj. R2: 0.686552
                Within R2: 0.146478

### EDD_6_9
OLS estimation, Dep. Var.: log_y
Observations: 20,753,599
Weights: size
Fixed-effects: OBJECTID: 2,799,471,  year: 25
Standard-errors: Clustered (FIPS) 
                                   Estimate Std. Error    t value   Pr(>|t|)    
bin_middle                        -0.074451   0.005676 -13.117687  < 2.2e-16 ***
bin_high                          -0.213484   0.010244 -20.840017  < 2.2e-16 ***
duration_years_low                 0.000543   0.000298   1.823375 6.8630e-02 .  
duration_years_middle              0.000692   0.000223   3.111463 1.9296e-03 ** 
duration_years_high                0.003641   0.000269  13.519700  < 2.2e-16 ***
weather_c_low                     -0.006965   0.001515  -4.597510 4.9885e-06 ***
weather_c_middle                  -0.006017   0.000436 -13.813871  < 2.2e-16 ***
weather_c_high                    -0.006539   0.000217 -30.189283  < 2.2e-16 ***
duration_years_x_weather_c_low    -0.000762   0.000148  -5.129352 3.6737e-07 ***
duration_years_x_weather_c_middle  0.000025   0.000029   0.882584 3.7773e-01    
duration_years_x_weather_c_high    0.000098   0.000016   6.049604 2.2547e-09 ***
GDD_4_5                            0.000588   0.000038  15.547294  < 2.2e-16 ***
ppt_4_5                           -0.000119   0.000014  -8.494109  < 2.2e-16 ***
GDD_6_9                            0.000061   0.000025   2.473001 1.3611e-02 *  
ppt_6_9                            0.000073   0.000010   7.325652 5.9539e-13 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
RMSE: 1.30701     Adj. R2: 0.684817
                Within R2: 0.141754


