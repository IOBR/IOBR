# Build Prognostic Models Using LASSO and Ridge Regression

Prepares data, splits it into training and testing sets, and fits LASSO
and Ridge regression models for survival analysis. Evaluates model
performance using cross-validation and optionally generates
time-dependent ROC curves for visual assessment of predictive accuracy.

## Usage

``` r
PrognosticModel(
  x,
  y,
  scale = FALSE,
  seed = 123456,
  train_ratio = 0.7,
  nfold = 10,
  plot = TRUE,
  palette = "jama",
  cols = NULL
)
```

## Arguments

- x:

  A matrix or data frame of predictor variables (features).

- y:

  A data frame of survival outcomes with two columns: survival time and
  event status.

- scale:

  Logical indicating whether to scale predictor variables. Default is
  \`FALSE\`.

- seed:

  Integer seed for random number generation to ensure reproducibility.
  Default is \`123456\`.

- train_ratio:

  Numeric proportion of data for training (e.g., 0.7). Default is
  \`0.7\`.

- nfold:

  Integer number of folds for cross-validation. Default is \`10\`.

- plot:

  Logical indicating whether to plot ROC curves. Default is \`TRUE\`.

- palette:

  String specifying color palette. Default is \`"jama"\`.

- cols:

  Optional vector of colors for ROC curves. If \`NULL\`, uses default
  palette.

## Value

A list containing:

- lasso_result:

  Results from LASSO model including coefficients and AUC

- ridge_result:

  Results from Ridge model including coefficients and AUC

- train.x:

  Training data with sample IDs

## Author

Dongqiang Zeng

## Examples

``` r
if (requireNamespace("glmnet", quietly = TRUE) &&
  requireNamespace("survival", quietly = TRUE)) {
  library(survival)
  set.seed(123)
  # Create small example data (first column must be ID)
  x_sim <- as.data.frame(matrix(rnorm(100 * 5), 100, 5))
  colnames(x_sim) <- paste0("Sig", 1:5)
  x_sim$ID <- paste0("S", 1:100)
  x_sim <- x_sim[, c(6, 1:5)] # Move ID to first column

  y_sim <- data.frame(
    ID = paste0("S", 1:100),
    OS_days = rexp(100, 0.01),
    OS_status = rbinom(100, 1, 0.5)
  )
  prognostic_result <- PrognosticModel(
    x = x_sim, y = y_sim,
    scale = TRUE, seed = 123456,
    train_ratio = 0.7, nfold = 3, plot = FALSE
  )
  if (!is.null(prognostic_result)) head(prognostic_result)
}
#> ℹ Processing data
#> ℹ Splitting data into training and validation sets
#> ℹ Running LASSO
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> ℹ Running RIDGE REGRESSION
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> ✔ Model fitting complete
#> $lasso_result
#> $lasso_result$model
#> 
#> Call:  glmnet::cv.glmnet(x = train.x, y = as.matrix(train.y), nfolds = nfold,      family = "cox", alpha = 1) 
#> 
#> Measure: Partial Likelihood Deviance 
#> 
#>      Lambda Index Measure     SE Nonzero
#> min 0.05584    14   3.876 0.1870       3
#> 1se 0.18715     1   3.968 0.1557       0
#> 
#> $lasso_result$coefs
#>   feature  lambda.min lambda.1se
#> 1    Sig1 -0.21600149          0
#> 2    Sig2  0.31941492          0
#> 3    Sig3  0.00000000          0
#> 4    Sig4  0.08920698          0
#> 5    Sig5  0.00000000          0
#> 
#> $lasso_result$AUC
#>                    probs.3   probs.9
#> Train_lambda.min 0.6540505 0.7627230
#> Train_lambda.1se 0.5000000 0.5000000
#> Test_lambda.min  0.4159452 0.2270948
#> Test_lambda.1se  0.5000000 0.5000000
#> 
#> 
#> $ridge_result
#> $ridge_result$model
#> 
#> Call:  glmnet::cv.glmnet(x = train.x, y = as.matrix(train.y), nfolds = nfold,      family = "cox", alpha = 0) 
#> 
#> Measure: Partial Likelihood Deviance 
#> 
#>     Lambda Index Measure     SE Nonzero
#> min   0.25    72   3.887 0.1889       5
#> 1se 187.15     1   3.973 0.1599       5
#> 
#> $ridge_result$coefs
#>   feature  lambda.min    lambda.1se
#> 1    Sig1 -0.21698782 -1.695910e-37
#> 2    Sig2  0.28062189  1.935739e-37
#> 3    Sig3  0.04159439  6.573543e-38
#> 4    Sig4  0.13033033  8.702333e-38
#> 5    Sig5 -0.04593428 -2.337453e-38
#> 
#> $ridge_result$AUC
#>                    probs.3   probs.9
#> Train_lambda.min 0.6523896 0.8154881
#> Train_lambda.1se 0.6482388 0.8120806
#> Test_lambda.min  0.3221501 0.2101121
#> Test_lambda.1se  0.3795094 0.2277480
#> 
#> 
#> $train.x
#>     ID        Sig1        Sig2         Sig3        Sig4        Sig5
#> 60 S60  0.13752572 -0.27615073 -1.740608199  2.51032748 -0.40532716
#> 42 S42 -0.32681639 -0.15993054 -0.655519126  0.95807304  0.35267582
#> 71 S71 -0.63697081 -0.11151520  0.762779177  2.36141760  1.88970138
#> 54 S54  1.40027842 -0.93158458 -0.031549770  0.14715340  0.61588929
#> 74 S74 -0.87597805  2.31233678 -1.595718525  0.82943622 -0.34891996
#> 3   S3  1.60854170 -0.14389556 -0.405957160 -0.86862929 -0.74857465
#> 67 S67  0.39197814  0.76952097 -0.266998306  0.39892007  0.64188789
#> 23 S23 -1.22304003 -0.39608682  1.041588525 -1.25392272  0.02068532
#> 94 S94 -0.78691881 -0.81471300 -0.734760046 -0.76363982 -1.20872075
#> 38 S38 -0.16686565  0.87642818  0.314268771 -0.15841460  0.98946613
#> 2   S2 -0.35120270  0.37687234  1.254841839 -0.68971795 -1.28822738
#> 80 S80 -0.25119772 -0.36279565 -0.100844561 -1.68857339 -0.33712681
#> 85 S85 -0.34058618 -0.13312777 -0.002443047  0.49322819  0.61352804
#> 65 S65 -1.27319995 -0.31987081  2.287253372  0.82200403 -1.21154619
#> 56 S56  1.56226982 -0.17874967 -0.220059850  0.35630780 -2.79656893
#> 24 S24 -0.89754917 -0.15361681  0.618103217  0.67102022  1.17897691
#> 46 S46 -1.32941869 -0.43781343  1.875919753 -1.65291248  0.06011053
#> 14 S14  0.02221347  0.05375962 -0.670115841 -0.54417984 -0.23102093
#> 79 S79  0.09957931  0.56264510  0.323937751  0.02109732  0.27003069
#> 35 S35  0.80101058 -2.01212758 -0.731518252  0.39816939 -0.08752415
#> 89 S89 -0.45610238  0.89101605 -0.357430264 -0.38451076  1.53741016
#> 97 S97  2.29720706  0.73243581  1.838536067 -0.71037289  0.07681467
#> 92 S92  0.50173432 -0.22455235  0.982984094 -2.04222350  1.31217957
#> 70 S70  2.14685001  0.49277965 -0.032907231 -0.40895563  1.03012700
#> 63 S63 -0.46407310 -1.19195906 -0.077178967  0.29841644  0.41709667
#> 48 S48 -0.61026684  0.82262108 -1.558415198 -0.51563037  0.78582422
#> 29 S29 -1.34588242 -0.88347639  0.804508172  1.32052880  1.54416938
#> 84 S84  0.60688103 -0.78384339 -0.406490458 -0.98075910  0.01563402
#> 91 S91  0.98935390  0.33298507  1.103212938 -0.58801864 -0.21391418
#> 82 S82  0.32303830  1.41752943  0.648536159  0.21799896  0.21044890
#> 50 S50 -0.19037243 -1.21975182  0.384254082 -0.13831914 -1.76099216
#> 62 S62 -0.64934164 -0.97377792 -0.642540351  0.66175261 -0.45921595
#> 30 S30  1.27452758  0.03747592 -1.196002577  0.04188851 -1.56098865
#> 10 S10 -0.58726835  1.06159010 -0.183700494  1.99637506  0.12274879
#> 75 S75 -0.85276181 -0.65542717  0.767653526  0.52595815 -0.32333706
#> 4   S4 -0.02179795 -0.24818937  0.445034541 -0.97834881 -0.13614297
#> 26 S26 -1.94683206 -0.56298929 -0.063918852  0.69308957 -0.56217912
#> 18 S18 -2.25349176 -0.55136153 -0.888040594 -1.17678842  0.32792758
#> 11 S11  1.24195461 -0.48377109 -0.001284241  1.28746927  0.39053852
#> 87 S87  1.10255872  1.25903201  0.105860180  0.47407414  0.81763358
#> 44 S44  2.27707482 -1.45516065 -1.324135788  0.08879602 -0.30903326
#> 88 S88  0.37770550  0.19884876  1.600605039 -1.05832252 -0.68757100
#> 7   S7  0.40589817 -0.70048299 -0.957035562 -1.90414261 -0.46051155
#> 16 S16  1.85854263  0.42265338  1.637294730  0.71212956 -2.77860368
#> 83 S83 -0.50510289 -0.25036912  0.279573997  0.20307385  1.23552733
#> 43 S43 -1.48529653 -1.51459943  1.448179833  1.64353927  0.56001826
#> 57 S57 -1.79571669  0.69342877  1.011007110  0.43094043  1.01524277
#> 32 S32 -0.42229479  0.57813712 -0.221906908 -1.10919518  0.27550830
#> 13 S13  0.34000892 -1.56189951  1.170686804 -1.62739512  0.55330174
#> 55 S55 -0.34637532 -0.01231228  1.556033681 -0.82499052  0.82006109
#> 64 S64 -1.21490140  3.46290911  1.241983044  1.02128929 -0.50188307
#> 59 S59  0.03664303  1.12154623 -0.246457722 -2.33896727  0.12610950
#> 28 S28  0.06898128  0.19184097 -0.881884222  0.64420073 -0.09574179
#> 22 S22 -0.33783464 -0.86860336  0.523772739  1.20750307  0.18251808
#> 66 S66  0.23347834  0.41962772  1.502418732 -0.16709030  1.11605021
#> 12 S12  0.29513939  0.73993902  0.129724233  0.76339238  0.16372849
#> 34 S34  0.86296437 -0.32570258 -0.904318179  1.49737075 -1.12345672
#> 47 S47 -0.54040552 -1.40044213 -0.233124433  0.13048994  0.33655348
#> 33 S33  0.88157948  0.15385913  0.099037583 -0.65979392  0.19627881
#> 77 S77 -0.41101270  0.15029701  0.057204752 -0.92469695  1.62375385
#> 98 S98  1.57995139 -1.18277188 -1.549005559  0.84902810  0.05962520
#> 69 S69  0.91131364  0.64572646 -0.536115646  0.85980176 -0.04113526
#> 61 S61  0.31685861  1.19987005 -0.675435978 -0.16276416  0.77436426
#> 5   S5  0.04259548 -0.87288879 -0.563024425 -0.38596830  0.57092782
#> 41 S41 -0.86009994  0.83696209 -0.957055704  0.63575502 -1.36140005
#> 81 S81 -0.09272595 -0.98841011 -1.882282138  0.06813174 -0.08631989
#> 27 S27  0.81876439  0.35464128 -0.868596532 -0.02368068  2.31627966
#> 31 S31  0.36815564  1.60508703  1.931644899  1.01444031 -0.29955929
#> 76 S76  1.02448422 -1.02219561 -0.596941661 -0.53260313  0.04632311
#> 96 S96 -0.75663177  2.17661772 -0.056601147  1.48262549 -1.47969495
#> 
#> $plots
#> NULL
#> 
```
