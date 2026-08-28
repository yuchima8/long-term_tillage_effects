# 20-60-20 Binning Estimator for Conditional Low-Till Duration Effects on corn Yield (With Other Weather Controls)

## Design
- Outcome: log corn yield in bushels per acre (`log_y`).
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
1:                                  full_cleaned_panel     49734268 4287518      798    25                    NA                     NA
2: clean_controls_and_continuous_low_till_with_weather     26014018 3741900      794    25              11697356               7.746795
3:                                   estimation_sample     23391541 3381656      794    25              10753030               7.860792

## Conditional Low-Till Duration Effects by Weather Bin
    weather_var weather_bin weather_median observations estimate_log_points_per_year    std_error pct_effect_per_year pct_effect_low_per_year
         <char>      <char>          <num>        <int>                        <num>        <num>               <num>                   <num>
 1:     GDD_4_5         low     247.536637      5848291                  0.002568866 0.0004451148           0.2572169              0.16978811
 2:     GDD_4_5      middle     338.633145     11696984                  0.001634559 0.0002928077           0.1635895              0.10612182
 3:     GDD_4_5        high     458.704294      5846266                  0.002091745 0.0003142612           0.2093934              0.14768822
 4:     PPT_4_5         low     120.280219      5848028                  0.002162536 0.0003303691           0.2164876              0.15161604
 5:     PPT_4_5      middle     188.134987     11695851                  0.001737914 0.0002919530           0.1739425              0.11663656
 6:     PPT_4_5        high     268.433457      5847662                  0.001226971 0.0003438293           0.1227724              0.05532187
 7:     GDD_6_9         low    1325.957798      5847847                  0.003003493 0.0004972808           0.3008008              0.20308822
 8:     GDD_6_9      middle    1492.703074     11695752                  0.001674035 0.0003234827           0.1675437              0.10405497
 9:     GDD_6_9        high    1712.824278      5847942                  0.002203233 0.0003067213           0.2205662              0.16033436
10:     PPT_6_9         low     261.923568      5848483                  0.003146444 0.0002933404           0.3151399              0.25748060
11:     PPT_6_9      middle     367.783878     11694562                  0.001959933 0.0002833540           0.1961855              0.14055461
12:     PPT_6_9        high     514.928912      5848496                  0.001240349 0.0003157038           0.1241119              0.06217632
13:     EDD_6_9         low       2.099847      5848373                  0.002299211 0.0004426565           0.2301857              0.14326299
14:     EDD_6_9      middle      10.024067     11694488                  0.001052467 0.0003083052           0.1053021              0.04482889
15:     EDD_6_9        high      28.373285      5848680                  0.002960121 0.0003048358           0.2964507              0.23654363
    pct_effect_high_per_year
                       <num>
 1:                0.3447219
 2:                0.2210902
 3:                0.2711366
 4:                0.2814011
 5:                0.2312812
 6:                0.1902685
 7:                0.3986087
 8:                0.2310726
 9:                0.2808343
10:                0.3728324
11:                0.2518473
12:                0.1860858
13:                0.3171838
14:                0.1658118
15:                0.3563935

## Additive Other-Weather Control Terms
    weather_var focal_weather_source_var additive_control      estimate    std_error        ci_low       ci_high
         <char>                   <char>           <char>         <num>        <num>         <num>         <num>
 1:     GDD_4_5                  GDD_4_5          ppt_4_5 -1.998422e-04 1.409585e-05 -2.274700e-04 -1.722143e-04
 2:     GDD_4_5                  GDD_4_5          GDD_6_9  9.945798e-05 3.623008e-05  2.844703e-05  1.704689e-04
 3:     GDD_4_5                  GDD_4_5          ppt_6_9 -5.110274e-06 1.144979e-05 -2.755187e-05  1.733132e-05
 4:     GDD_4_5                  GDD_4_5          EDD_6_9 -9.411864e-03 2.393710e-04 -9.881031e-03 -8.942696e-03
 5:     PPT_4_5                  ppt_4_5          GDD_4_5  2.676139e-04 4.270663e-05  1.839089e-04  3.513189e-04
 6:     PPT_4_5                  ppt_4_5          GDD_6_9  9.552616e-05 3.635740e-05  2.426566e-05  1.667867e-04
 7:     PPT_4_5                  ppt_4_5          ppt_6_9 -5.763005e-06 1.111507e-05 -2.754853e-05  1.602252e-05
 8:     PPT_4_5                  ppt_4_5          EDD_6_9 -9.177376e-03 2.350611e-04 -9.638096e-03 -8.716657e-03
 9:     GDD_6_9                  GDD_6_9          GDD_4_5  2.762242e-04 4.439298e-05  1.892140e-04  3.632345e-04
10:     GDD_6_9                  GDD_6_9          ppt_4_5 -2.068517e-04 1.413807e-05 -2.345623e-04 -1.791411e-04
11:     GDD_6_9                  GDD_6_9          ppt_6_9 -6.507938e-06 1.138594e-05 -2.882438e-05  1.580851e-05
12:     GDD_6_9                  GDD_6_9          EDD_6_9 -9.495538e-03 2.594829e-04 -1.000412e-02 -8.986951e-03
13:     PPT_6_9                  ppt_6_9          GDD_4_5  2.288623e-04 4.182639e-05  1.468826e-04  3.108420e-04
14:     PPT_6_9                  ppt_6_9          ppt_4_5 -2.162615e-04 1.424541e-05 -2.441825e-04 -1.883405e-04
15:     PPT_6_9                  ppt_6_9          GDD_6_9  1.475953e-04 3.644550e-05  7.616213e-05  2.190285e-04
16:     PPT_6_9                  ppt_6_9          EDD_6_9 -9.259853e-03 2.302874e-04 -9.711216e-03 -8.808490e-03
17:     EDD_6_9                  EDD_6_9          GDD_4_5  2.725779e-04 4.494466e-05  1.844863e-04  3.606694e-04
18:     EDD_6_9                  EDD_6_9          ppt_4_5 -2.036324e-04 1.424800e-05 -2.315585e-04 -1.757063e-04
19:     EDD_6_9                  EDD_6_9          GDD_6_9  4.033963e-05 3.705391e-05 -3.228604e-05  1.129653e-04
20:     EDD_6_9                  EDD_6_9          ppt_6_9  1.565854e-05 1.103975e-05 -5.979365e-06  3.729644e-05
    weather_var focal_weather_source_var additive_control      estimate    std_error        ci_low       ci_high

## Model Summaries
### GDD_4_5
OLS estimation, Dep. Var.: log_y
Observations: 22,970,456
Weights: size
Fixed-effects: OBJECTID: 2,960,571,  year: 25
Standard-errors: Clustered (FIPS) 
                                     Estimate Std. Error    t value   Pr(>|t|)    
bin_middle                         0.02251876 0.00442656   5.087188 4.5535e-07 ***
bin_high                           0.05553069 0.00951882   5.833779 7.9308e-09 ***
duration_years_low                 0.00256887 0.00044511   5.771244 1.1341e-08 ***
duration_years_middle              0.00163456 0.00029281   5.582362 3.2732e-08 ***
duration_years_high                0.00209174 0.00031426   6.656070 5.2915e-11 ***
weather_c_low                      0.00020558 0.00006561   3.133427 1.7922e-03 ** 
weather_c_middle                   0.00010689 0.00005962   1.792998 7.3360e-02 .  
weather_c_high                     0.00022471 0.00007266   3.092613 2.0545e-03 ** 
duration_years_x_weather_c_low    -0.00000463 0.00000736  -0.629538 5.2918e-01    
duration_years_x_weather_c_middle  0.00000258 0.00000508   0.507045 6.1227e-01    
duration_years_x_weather_c_high    0.00001273 0.00000473   2.691741 7.2599e-03 ** 
ppt_4_5                           -0.00019984 0.00001410 -14.177378  < 2.2e-16 ***
GDD_6_9                            0.00009946 0.00003623   2.745177 6.1866e-03 ** 
ppt_6_9                           -0.00000511 0.00001145  -0.446320 6.5549e-01    
EDD_6_9                           -0.00941186 0.00023937 -39.319146  < 2.2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
RMSE: 1.47357     Adj. R2: 0.688428
                Within R2: 0.186997

### PPT_4_5
OLS estimation, Dep. Var.: log_y
Observations: 22,970,456
Weights: size
Fixed-effects: OBJECTID: 2,960,571,  year: 25
Standard-errors: Clustered (FIPS) 
                                     Estimate Std. Error    t value   Pr(>|t|)    
bin_middle                        -0.00183678 0.00264246  -0.695101 4.8720e-01    
bin_high                          -0.02164763 0.00322649  -6.709340 3.7514e-11 ***
duration_years_low                 0.00216254 0.00033037   6.545818 1.0704e-10 ***
duration_years_middle              0.00173791 0.00029195   5.952718 3.9802e-09 ***
duration_years_high                0.00122697 0.00034383   3.568547 3.8085e-04 ***
weather_c_low                      0.00029626 0.00010521   2.816025 4.9846e-03 ** 
weather_c_middle                  -0.00036952 0.00006410  -5.764854 1.1761e-08 ***
weather_c_high                    -0.00036808 0.00004064  -9.056971  < 2.2e-16 ***
duration_years_x_weather_c_low    -0.00003619 0.00001163  -3.112231 1.9243e-03 ** 
duration_years_x_weather_c_middle  0.00000266 0.00000589   0.451789 6.5155e-01    
duration_years_x_weather_c_high    0.00002024 0.00000396   5.112443 4.0020e-07 ***
GDD_4_5                            0.00026761 0.00004271   6.266331 6.1022e-10 ***
GDD_6_9                            0.00009553 0.00003636   2.627420 8.7721e-03 ** 
ppt_6_9                           -0.00000576 0.00001112  -0.518486 6.0427e-01    
EDD_6_9                           -0.00917738 0.00023506 -39.042518  < 2.2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
RMSE: 1.47127     Adj. R2: 0.6894  
                Within R2: 0.189532

### GDD_6_9
OLS estimation, Dep. Var.: log_y
Observations: 22,970,456
Weights: size
Fixed-effects: OBJECTID: 2,960,571,  year: 25
Standard-errors: Clustered (FIPS) 
                                     Estimate Std. Error    t value   Pr(>|t|)    
bin_middle                         0.02068499 0.00626265   3.302911 1.0004e-03 ** 
bin_high                           0.03932311 0.01499423   2.622550 8.8973e-03 ** 
duration_years_low                 0.00300349 0.00049728   6.039833 2.3839e-09 ***
duration_years_middle              0.00167403 0.00032348   5.175037 2.8990e-07 ***
duration_years_high                0.00220323 0.00030672   7.183176 1.5890e-12 ***
weather_c_low                      0.00003026 0.00004070   0.743505 4.5740e-01    
weather_c_middle                   0.00014107 0.00004298   3.281950 1.0767e-03 ** 
weather_c_high                     0.00022116 0.00007118   3.107115 1.9575e-03 ** 
duration_years_x_weather_c_low    -0.00001322 0.00000450  -2.934961 3.4336e-03 ** 
duration_years_x_weather_c_middle -0.00000600 0.00000248  -2.421494 1.5683e-02 *  
duration_years_x_weather_c_high    0.00000359 0.00000309   1.160046 2.4638e-01    
GDD_4_5                            0.00027622 0.00004439   6.222250 7.9824e-10 ***
ppt_4_5                           -0.00020685 0.00001414 -14.630827  < 2.2e-16 ***
ppt_6_9                           -0.00000651 0.00001139  -0.571577 5.6777e-01    
EDD_6_9                           -0.00949554 0.00025948 -36.594074  < 2.2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
RMSE: 1.47299     Adj. R2: 0.688674
                Within R2: 0.187638

### PPT_6_9
OLS estimation, Dep. Var.: log_y
Observations: 22,970,456
Weights: size
Fixed-effects: OBJECTID: 2,960,571,  year: 25
Standard-errors: Clustered (FIPS) 
                                     Estimate Std. Error    t value   Pr(>|t|)    
bin_middle                         0.02253788 0.00304201   7.408883 3.3068e-13 ***
bin_high                           0.02324002 0.00344124   6.753378 2.8179e-11 ***
duration_years_low                 0.00314644 0.00029334  10.726254  < 2.2e-16 ***
duration_years_middle              0.00195993 0.00028335   6.916907 9.6043e-12 ***
duration_years_high                0.00124035 0.00031570   3.928839 9.2906e-05 ***
weather_c_low                      0.00005597 0.00007740   0.723128 4.6982e-01    
weather_c_middle                   0.00020080 0.00003496   5.743358 1.3287e-08 ***
weather_c_high                    -0.00029431 0.00002803 -10.499836  < 2.2e-16 ***
duration_years_x_weather_c_low    -0.00000578 0.00000755  -0.765584 4.4416e-01    
duration_years_x_weather_c_middle -0.00001472 0.00000357  -4.125652 4.0932e-05 ***
duration_years_x_weather_c_high    0.00000712 0.00000252   2.827082 4.8174e-03 ** 
GDD_4_5                            0.00022886 0.00004183   5.471720 6.0041e-08 ***
ppt_4_5                           -0.00021626 0.00001425 -15.181142  < 2.2e-16 ***
GDD_6_9                            0.00014760 0.00003645   4.049754 5.6380e-05 ***
EDD_6_9                           -0.00925985 0.00023029 -40.209990  < 2.2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
RMSE: 1.46463     Adj. R2: 0.692199
                Within R2: 0.196836

### EDD_6_9
OLS estimation, Dep. Var.: log_y
Observations: 22,970,456
Weights: size
Fixed-effects: OBJECTID: 2,960,571,  year: 25
Standard-errors: Clustered (FIPS) 
                                   Estimate Std. Error    t value   Pr(>|t|)    
bin_middle                        -0.067906   0.004742 -14.319521  < 2.2e-16 ***
bin_high                          -0.220353   0.010021 -21.989406  < 2.2e-16 ***
duration_years_low                 0.002299   0.000443   5.194121 2.6258e-07 ***
duration_years_middle              0.001052   0.000308   3.413717 6.7382e-04 ***
duration_years_high                0.002960   0.000305   9.710543  < 2.2e-16 ***
weather_c_low                     -0.009661   0.001951  -4.951988 9.0039e-07 ***
weather_c_middle                  -0.007335   0.000564 -13.011881  < 2.2e-16 ***
weather_c_high                    -0.010604   0.000363 -29.180106  < 2.2e-16 ***
duration_years_x_weather_c_low    -0.000534   0.000194  -2.750628 6.0856e-03 ** 
duration_years_x_weather_c_middle -0.000020   0.000036  -0.550402 5.8220e-01    
duration_years_x_weather_c_high    0.000132   0.000031   4.332924 1.6633e-05 ***
GDD_4_5                            0.000273   0.000045   6.064744 2.0565e-09 ***
ppt_4_5                           -0.000204   0.000014 -14.291996  < 2.2e-16 ***
GDD_6_9                            0.000040   0.000037   1.088674 2.7663e-01    
ppt_6_9                            0.000016   0.000011   1.418378 1.5648e-01    
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
RMSE: 1.46879     Adj. R2: 0.690447
                Within R2: 0.192265


