#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Supplementary material:
#
# Paper title: Incorporating environmental covariates to explore genotype ? environment ? 
#              management (G?E?M) interactions: A one-stage predictive model
#
# R-script and data to reproduce the statistical models, predictions and figures 
# in the core manuscript
#
# Note that this code for ASReml-R version 4.1 and does not currently work for ASReml-R version 4.2
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 1.1. Set the working directory ----
# By default it will be the same location as this file when running the line of code below
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# 1.2. Load the required R-packages ----
require(asreml) # Need to have an accompanying ASReml-R license
require(lattice)
require(tidyverse)
require(plot3D)
require(plotly)
require(viridis)

# Load required source code
# This source code is available as supplementary material from Verbyla (2019).
source('icREML.R')

# 1.3. Load data file ----
Sorghum_df <- readRDS(file="Sorghum_yield_data.Rds")

## Note that each of the environmental covariates are zero centred around the mean value across all environments (See Table 2 in manuscript)

# Round density to 2 decimal places
Sorghum_df$density <- round(Sorghum_df$density, 2)

# 2.1 Run analysis for baseline model ----
Baseline_fm <- asreml( Yld ~ Genotype + density + Genotype:density,
  random =~ at(Trial):Rep  + at(Trial):MainPlot + 
    at(Trial,c(1:3)):SubPlot +
    at(Trial,'Breeza 1'):Column + 
    Trial + Env +
    spl(density, k=6) + spl(density, k=6):Genotype +
   str(~Trial:Genotype + Trial:Genotype:density,
     ~corh(2):id(48)) +
   str(~Env:Genotype + Env:Genotype:density,
     ~corh(2):id(136)),
  residual=~ dsum(~units|ResidualTOS),
  data = Sorghum_df, 
  na.action=na.method(x='include'),
  maxit=30, workspace="1Gb")


# 2.2 Obtain prediction from baseline model ----

# Create data frame of values that we want to predict for
aux_parallel_Init <- unique(Sorghum_df[,c('Trial', 'Env', 'Genotype', 'density')])
# Reorder the data frame of values that we want to predict for
aux_parallel_Init <- aux_parallel_Init[order( aux_parallel_Init$Env, aux_parallel_Init$Genotype),]


# Obtain predictions for the GxExM interaction effect from the baseline model
Baseline_GxExM_pred <- predict( Baseline_fm, classify='Trial:Env:Genotype:density', 
                        levels=list('Trial'=aux_parallel_Init$Trial,
                                    'Env'=aux_parallel_Init$Env,
                                    'Genotype'=aux_parallel_Init$Genotype,
                                    'density'=aux_parallel_Init$density),
                        parallel=T)

# Obtain predictions for the GxM interaction effect from the baseline model (averaged across all environments)
Baseline_GxM_pred <- predict( Baseline_fm, classify='Genotype:density', 
                         levels=list('Genotype'=aux_parallel_Init$Genotype,
                                     'density'=aux_parallel_Init$density),
                         parallel=T)


# 2.3 Store predicted values from the baseline model -----
# Store predicted values from the Baseline model under a unique name
Baseline_GxExM_pred$pvals$Yld_pred_baseline_GxExM <- Baseline_GxExM_pred$pvals$predicted.value
# Include yield predictions from the Baseline model with those in the original data frame
Sorghum_df <- Sorghum_df %>% left_join(Baseline_GxExM_pred$pvals[,c(1:4,8)],
                                       by=c("Trial","Env","Genotype","density"))



# 3.1 Adding Environmental Covariates ----

# Define matrix with environmental covariate information
Env_Var <- matrix(nrow=18,ncol=8)
rownames(Env_Var) <- colnames(Sorghum_df)[14:31]
colnames(Env_Var)  <- c("full_loglik","p","q","b","AIC","BIC",
                        "Cor_Yld_Pred","RMSE_Yld_Pred")


# 3.2 Create a for loop to perform the forward selection procedure ----
# Note that in order to calculate RMSEP for the EC "Irrig", need to drop spline terms from the model
# Also note that this for-loop will take around 30 minutes to run
for(i in c(1:18)) {
  
  print(c('Iteration ',i))
  # Allocate relevant EC to the 'ec' term
  Sorghum_df$ec <- Sorghum_df[,i+13]
  
  # Extend baseline model to incorporate full model for a single EC
  Temp_fm <- asreml( Yld ~ Genotype + density + Genotype:density + 
                       ec + Genotype:ec + density:ec + Genotype:density:ec,
                     random =~ at(Trial):Rep  + at(Trial):MainPlot +
                       at(Trial,c(1:3)):SubPlot +
                       at(Trial,'Breeza 1'):Column +
                       Trial + Env +
                       spl(density, k=6) + spl(density, k=6):Genotype +
                       str(~Trial:Genotype + Trial:Genotype:density,
                           ~corh(2):id(48)) +
                       str(~Env:Genotype + Env:Genotype:density,
                           ~corh(2):id(136)) +
                       spl(ec, k=6) + density:spl(ec, k=6) + spl(ec, k=6):Genotype +
                       density:spl(ec, k=6):Genotype + 
                       spl(density, k=6):ec + spl(density, k=6):ec:Genotype +
                       spl(density, k=6):spl(ec, k=6) + spl(density, k=6):spl(ec, k=6):Genotype,
                     residual=~ dsum(~units|ResidualTOS),
                     data = Sorghum_df,
                     na.action=na.method(x='include'),
                     maxit=30, workspace="2Gb")
  

  # Obtain AIC and BIC from Verbyla (2019) source code 
  ic_temp <- icREML(list(Temp_fm))
  # Change format to matrix
  ic_temp <- as.matrix(ic_temp)
  
  # Store the value of each component for each iteration in the matrix 'ic_temp'
  Env_Var[i,1] <- as.numeric(ic_temp[2])
  Env_Var[i,2] <- as.numeric(ic_temp[3])
  Env_Var[i,3] <- as.numeric(ic_temp[4])
  Env_Var[i,4] <- as.numeric(ic_temp[5])
  Env_Var[i,5] <- as.numeric(ic_temp[6])
  Env_Var[i,6] <- as.numeric(ic_temp[7])
  

  # Define a term to store the cross_validation predictions at each iteration
  CrossValid_pred <- c()
  
  # 3.3 Begin cross validation procedure to calculate RMSE ----
  # Begin for loop removing one trial at a time and rerunning the model for the candidate EC
  for (j in c(1:6)){
    
    # Subset one trial out of the data frame
    Sorghum_subset_df <- Sorghum_df[Sorghum_df$Trial!=levels(Sorghum_df$Trial)[j],]

    # Re-factorise each of the relevant terms to update the number of factor levels for each term
    Sorghum_subset_df$Trial <- factor(Sorghum_subset_df$Trial)
    Sorghum_subset_df$Env <- factor(Sorghum_subset_df$Env)
    Sorghum_subset_df$ResidualTOS <- factor(Sorghum_subset_df$ResidualTOS)
    Sorghum_subset_df$Genotype <- factor(Sorghum_subset_df$Genotype)
    # Run a slightly different version of the model based on which trial has been removed in the cross validation procedure
    if(j==1){
      Temp_subset_fm <- asreml( Yld ~ Genotype + density + Genotype:density + 
                                  ec + Genotype:ec + density:ec + Genotype:density:ec,
                                random =~ at(Trial):Rep  + at(Trial):MainPlot +
                                  at(Trial,c(1:3)):SubPlot +
                                  #at(Trial,'Breeza 1'):Column +
                                  Trial + Env +
                                  spl(density, k=6) + spl(density, k=6):Genotype +
                                  str(~Trial:Genotype + Trial:Genotype:density,
                                      ~corh(2):id(40)) +
                                  str(~Env:Genotype + Env:Genotype:density,
                                      ~corh(2):id(112)) +
                                  spl(ec, k=6) + density:spl(ec, k=6) + spl(ec, k=6):Genotype +
                                  density:spl(ec, k=6):Genotype + 
                                  spl(density, k=6):ec + spl(density, k=6):ec:Genotype +
                                  spl(density, k=6):spl(ec, k=6) + spl(density, k=6):spl(ec, k=6):Genotype,
                                residual=~ dsum(~units|ResidualTOS),
                                data = Sorghum_subset_df,
                                na.action=na.method(x='include'),
                                maxit=30, workspace="2Gb")
    } else if(j==3) {
      Temp_subset_fm <- asreml( Yld ~ Genotype + density + Genotype:density + 
                                  ec + Genotype:ec + density:ec + Genotype:density:ec,
                                random =~ at(Trial):Rep  + at(Trial):MainPlot +
                                  at(Trial,c(1:3)):SubPlot +
                                  at(Trial,'Breeza 1'):Column +
                                  Trial + Env +
                                  spl(density, k=6) + spl(density, k=6):Genotype +
                                  str(~Trial:Genotype + Trial:Genotype:density,
                                      ~corh(2):id(40)) +
                                  str(~Env:Genotype + Env:Genotype:density,
                                      ~corh(2):id(120)) +
                                  spl(ec, k=6) + density:spl(ec, k=6) + spl(ec, k=6):Genotype +
                                  density:spl(ec, k=6):Genotype + 
                                  spl(density, k=6):ec + spl(density, k=6):ec:Genotype +
                                  spl(density, k=6):spl(ec, k=6) + spl(density, k=6):spl(ec, k=6):Genotype,
                                residual=~ dsum(~units|ResidualTOS),
                                data = Sorghum_subset_df,
                                na.action=na.method(x='include'),
                                maxit=30, workspace="2Gb")
    } else {
      Temp_subset_fm <- asreml( Yld ~ Genotype + density + Genotype:density + 
                                  ec + Genotype:ec + density:ec + Genotype:density:ec,
                                random =~ at(Trial):Rep  + at(Trial):MainPlot +
                                  at(Trial,c(1:3)):SubPlot +
                                  at(Trial,'Breeza 1'):Column +
                                  Trial + Env +
                                  spl(density, k=6) + spl(density, k=6):Genotype +
                                  str(~Trial:Genotype + Trial:Genotype:density,
                                      ~corh(2):id(40)) +
                                  str(~Env:Genotype + Env:Genotype:density,
                                      ~corh(2):id(112)) +
                                  spl(ec, k=6) + density:spl(ec, k=6) + spl(ec, k=6):Genotype +
                                  density:spl(ec, k=6):Genotype + 
                                  spl(density, k=6):ec + spl(density, k=6):ec:Genotype +
                                  spl(density, k=6):spl(ec, k=6) + spl(density, k=6):spl(ec, k=6):Genotype,
                                residual=~ dsum(~units|ResidualTOS),
                                data = Sorghum_subset_df,
                                na.action=na.method(x='include'),
                                maxit=30, workspace="2Gb")
    }
    
    # Create data frame with the combination of all explanatory variables to obtain predictions for (including those in the quasi-untested environment)
    aux_parallel_temp <- unique(Sorghum_df[,c( 'Env', 'Trial', 'Genotype', 'density', 'ec')])
    # Re-order the dataframe
    aux_parallel_temp <- aux_parallel_temp[order( aux_parallel_temp$Env, aux_parallel_temp$Genotype,
                                                  aux_parallel_temp$density),]
    
    # Obtain predictions for all terms using the candidate ec as a surrogate for the environment term in the classify set
    Temp_subset_pred <- predict(  Temp_subset_fm, classify='Genotype:density:ec', 
                                  levels=list('Genotype'=aux_parallel_temp$Genotype,
                                              'ec'=aux_parallel_temp$ec,
                                              'density'=aux_parallel_temp$density),
                                  parallel=T, pworkspace="2Gb", maxit=30)
    
    # Round numeric terms for the data frame of (i) predictions and (ii) combination of variables to predict for
    # to ensure the merging of the data-frames is smooth
    Temp_subset_pred$pvals$ec <- round(Temp_subset_pred$pvals$ec,2)
    aux_parallel_temp$ec <-  round(aux_parallel_temp$ec,2)
    Temp_subset_pred$pvals$density <- round(Temp_subset_pred$pvals$density,2)
    aux_parallel_temp$density <-  round(aux_parallel_temp$density,2)
    
    # Merge data frame of (i) predictions and (ii) combination of variables to predict for
    Temp_subset_pred$pvals <- left_join(Temp_subset_pred$pvals, aux_parallel_temp, by=c("Genotype","density","ec")) 
    
    # Subset so that the table of predictions only contains predictions that were from a quasi-untested environment
    Temp_subset_pred <- Temp_subset_pred$pvals[Temp_subset_pred$pvals$Trial==levels(Sorghum_df$Trial)[j],]
    # Merge the current set of predictions from an untested environment with the previous set(s) of predictions from an untested environment
    CrossValid_pred <-  bind_rows(CrossValid_pred, Temp_subset_pred)
  }
  
  # 3.3 Merge the set of predictions for a quasi-untested environment and calculate the RMSEP ----
  
  # First, need to re-factorise Genotype since "GenoA" is not in the first trial (i.e. Breeza 1)
  CrossValid_pred$Genotype <- factor(CrossValid_pred$Genotype, levels=c("GenoA","GenoB","GenoC","GenoD",
                                                                        "GenoE","GenoF","GenoG","GenoH"))
  
  # Define a temporary dataframe consisting of the cross-validation predictions
  CV_temp <- left_join(CrossValid_pred, Sorghum_df[,c(1,2,11,12,32)], by=c("Genotype","density",
                                            "Env","Trial")) %>% drop_na(density)


  # (optional) calculate the Pearson correlation between the predicted values from the baseline model and
  # predictions from the current Candidate model using the Candidate EC as a surrogate for the Environment term
  Env_Var[i,7] <- cor(CV_temp$predicted.value, CV_temp$Yld_pred_baseline_GxExM, use="pairwise.complete.obs")
  
  # Calculate the root meet square error of prediction between the baseline model and the 
  # Candidate model using the Candidate EC as a surrogate for the Environment term
  CV_temp$e_Pred <- (CV_temp$Yld_pred_baseline_GxExM - CV_temp$predicted.value)^2
  Env_Var[i,8] <- sqrt(mean(CV_temp$e_Pred, na.rm=T))
}  

# Store the information from the first iteration of the forward selection procedure
Env_Var1 <- Env_Var
# This table matches the results provided in Table 3 in the manuscript
Env_Var1


# 3.4 Calculate yield values that are adjusted for spatial field trend for each individual plot------

#Use bind_rows to repeat the set of prediction 3 times (1 for each replicate block)
GxE_adj.pred <- bind_rows(Baseline_GxExM_pred$pvals ,  Baseline_GxExM_pred$pvals)
GxE_adj.pred <- bind_rows(GxE_adj.pred, Baseline_GxExM_pred$pvals)

# Define Rep term in the data frame of predicted values from the baseline model
GxE_adj.pred$Rep <- rep(1:3, each=1443)
# Factorise the Rep term
GxE_adj.pred$Rep <- factor(GxE_adj.pred$Rep)

# Round density values to 2 decimal place
Sorghum_df$density <- round(Sorghum_df$density, 2)
GxE_adj.pred$density <- round(GxE_adj.pred$density, 2)

# Define new data_frame that merges raw data with the data frame of predicted values from the baseline model
SorghumMET_Adj.df <- left_join(Sorghum_df, GxE_adj.pred[,c(1:5,8:9)],
                               by=c('Genotype','Trial','Env','Rep','density', 'Yld_pred_baseline_GxExM'))

# Define a column containing the fitted values for each individual plot within each trial
SorghumMET_Adj.df$fitted <- Baseline_fm$linear.predictors

# Calculate the value for the design terms (i.e. spatial adjustment for each plot based on blocking terms and random column terms) for each individual plot
SorghumMET_Adj.df$Design <- SorghumMET_Adj.df$fitted - SorghumMET_Adj.df$predicted.value
# Calculate the yield value for each individual plot which is adjusted to account for spatial field trend within each trial
SorghumMET_Adj.df$Yield_Adj <- Sorghum_df$Yld - SorghumMET_Adj.df$Design

# Check that the average value of the design terms is approximately zero
mean(SorghumMET_Adj.df$Design, na.rm=T)


# 3.5  Update baseline model to incorporate the most significant environmental covariate (PostPAW) ------

Yield_EC1_full_fm <- asreml( Yld ~ Genotype + density + Genotype:density +
                          PostPAW + Genotype:PostPAW + density:PostPAW +
                          Genotype:density:PostPAW,
                        random =~ at(Trial):Rep  + at(Trial):MainPlot + 
                          at(Trial,c(1:3)):SubPlot +
                          at(Trial,'Breeza 1'):Column + 
                          Trial + Env +
                          spl(density, k=6) + spl(density, k=6):Genotype +
                          str(~Trial:Genotype + Trial:Genotype:density,
                              ~corh(2):id(48)) +
                          str(~Env:Genotype + Env:Genotype:density,
                              ~corh(2):id(136)) +
                          spl(PostPAW, k=6) + density:spl(PostPAW, k=6) + spl(PostPAW, k=6):Genotype +
                          density:spl(PostPAW, k=6):Genotype + 
                          spl(density, k=6):PostPAW + spl(density, k=6):PostPAW:Genotype +
                          spl(density, k=6):spl(PostPAW, k=6) + spl(density, k=6):spl(PostPAW, k=6):Genotype,
                        residual=~ dsum(~units|ResidualTOS),
                        data = Sorghum_df, 
                        na.action=na.method(x='include'),
                        maxit=30, workspace=120e6)

# Final model for first EC after applying backwards selection procedure
Yield_EC1_fm <- asreml( Yld ~ Genotype + density + Genotype:density +
                        PostPAW + Genotype:PostPAW + density:PostPAW +
                        Genotype:density:PostPAW,
                      random =~ at(Trial):Rep  + at(Trial):MainPlot + 
                        at(Trial,c(1:3)):SubPlot +
                        at(Trial,'Breeza 1'):Column + 
                        Trial + Env +
                        spl(density, k=6) + spl(density, k=6):Genotype +
                        str(~Trial:Genotype + Trial:Genotype:density,
                            ~corh(2):id(48)) +
                        str(~Env:Genotype + Env:Genotype:density,
                            ~corh(2):id(136)) +
                        spl(PostPAW, k=6) + density:spl(PostPAW, k=6) + 
                        spl(PostPAW, k=6):Genotype + spl(density , k=6):PostPAW,
                      residual=~ dsum(~units|ResidualTOS),
                      data = Sorghum_df, 
                      na.action=na.method(x='include'),
                      maxit=30, workspace=120e6)

# Table of variance component values for the statistical model
summary(Yield_EC1_fm)$varcomp
# Summary of Wald table for final PostPAW model
#  Note that Genotype:density:PostPAW interaction is dropped since the term is non-significant
wald(Yield_EC1_fm, denDF='numeric', ssType='conditional', workspace=120e6)
# Summary of residual diagnostics
plot(Yield_EC1_fm, facet=F)

# AIC and BIC values for the baseline, full EC1 and final EC1 models derived using the full log-likelihood at the REML parameter estimates (Verbyla 2019)
IC_EC1_df <- icREML(list(Baseline_fm, Yield_EC1_full_fm, Yield_EC1_fm))
IC_EC1_df

# 4.1 Final model obtained after performing subset selection procedure -----

Yield_final_fm <- asreml( Yld ~ Genotype + density + Genotype:density + 
                            PostPAW + Genotype:PostPAW +
                            PrePAW + density:PrePAW +
                            ISW + density:ISW +
                            PreCumRad + density:PreCumRad +
                            ptq + Genotype:ptq +
                            PreFlwEvap + density:PreFlwEvap +
                            PostMaxT + density:PostMaxT,
                          random =~ at(Trial):Rep  + at(Trial):MainPlot +
                            at(Trial,c(1:3)):SubPlot +
                            at(Trial,'Breeza 1'):Column +
                            Trial + Env +
                            spl(density) + spl(density):Genotype +
                            str(~Trial:Genotype + Trial:Genotype:density,
                                ~corh(2):id(48)) +
                            str(~Env:Genotype + Env:Genotype:density,
                                ~corh(2):id(136)) +
                            spl(PostPAW) + 
                            spl(PostPAW):Genotype +
                            spl(density):ISW +
                            #spl(PreCumRad, k=6) + density:spl(PreCumRad, k=6) +
                            spl(ptq) + spl(ptq):Genotype +
                            spl(density):PreFlwEvap +
                            spl(PostMaxT) + density:spl(PostMaxT),
                          residual=~ dsum(~units|ResidualTOS),
                          data = Sorghum_df,
                          na.action=na.method(x='include'),
                          maxit=30, workspace="4gb")

# Table of variance component values for the final model
summary(Yield_final_fm)$varcomp
# Summary of Wald table for the final model
wald(Yield_final_fm, denDF='numeric', ssType='conditional', workspace="4Gb")
# Summary of residual diagnostics from the final model
plot(Yield_final_fm, facet=F)

# AIC and BIC values for the baseline, full EC1, final EC1 and the final overall models derived using the full log-likelihood at the REML parameter estimates (Verbyla 2019)
IC_final_df <- icREML(list(Baseline_fm, Yield_EC1_full_fm, Yield_EC1_fm, Yield_final_fm))
IC_final_df

# 4.2 Calculate raw data points adjusted for Post-flowering plant available water adjusted for all other ECs and lack-of-fit terms in the model -----

# Create data frame of values that we want to predict for
aux_parallel_final <- unique(Sorghum_df[,c('Env', 'Trial', 'Genotype', 'density', 'PostPAW','ISW', 'PreCumRad', 'PrePAW', 
                                           'ptq','PreFlwEvap','PostMaxT')])
                                            
# Reorder the data frame of values that we want to predict for
aux_parallel_final <- aux_parallel_final[order(aux_parallel_final$Env, aux_parallel_final$Genotype),]


# Obtain predictions including all terms (including Env lack-of-fit) term in the classify set
Final_model_pred <- predict( Yield_final_fm , classify='Genotype:Trial:Env:density:PostPAW:ISW:PreCumRad:PrePAW:ptq:PreFlwEvap:PostMaxT', 
                    levels=list('Genotype'=aux_parallel_final$Genotype,
                                'Trial'=aux_parallel_final$Trial,
                                'Env'=aux_parallel_final$Env,
                                'density'=aux_parallel_final$density,
                                'PostPAW'=aux_parallel_final$PostPAW,
                                'ISW'=aux_parallel_final$ISW,
                                'PreCumRad'=aux_parallel_final$PreCumRad,
                                'PrePAW'=aux_parallel_final$PrePAW,
                                'ptq'=aux_parallel_final$ptq,
                                'PreFlwEvap'=aux_parallel_final$PreFlwEvap,
                                'PostMaxT'=aux_parallel_final$PostMaxT),
                    parallel=T)

# Round numbers to 2 decimal places to ensure the data frames merge correctly
Final_model_pred$pvals$density <- round(Final_model_pred$pvals$density, 2)
Final_model_pred$pvals$PostPAW <- round(Final_model_pred$pvals$PostPAW,2)

# Obtain predictions for Post-flowering plant avaliable water for each genotype and all established plant population densities observed
PostPAW_temp_pred <- predict( Yield_final_fm, classify='Genotype:density:PostPAW:ISW:PreCumRad:PrePAW:ptq:PreFlwEvap:PostMaxT', 
                             levels=list('Genotype'=aux_parallel_final$Genotype,
                                         'PostPAW'=aux_parallel_final$PostPAW,
                                         'ISW'=rep(0,1443),
                                         'PreCumRad'=rep(0,1443),
                                         'PrePAW'=rep(0,1443),
                                         'ptq'=rep(0,1443),
                                         'PreFlwEvap'=rep(0,1443),
                                         'PostMaxT'=rep(0,1443),
                                         'density'=aux_parallel_final$density),
                             parallel=T)
# Round numbers to 2 decimal places to ensure the data frames merge correctly
PostPAW_temp_pred$pvals$density <- round(PostPAW_temp_pred$pvals$density,2)
PostPAW_temp_pred$pvals$PostPAW <- round(PostPAW_temp_pred$pvals$PostPAW,2)

PostPAW_temp_pred$pvals$pred_value_PostPAW <- PostPAW_temp_pred$pvals$predicted.value

# Now merge the two data-frames together
PostPAW_pred_df <- left_join(Final_model_pred$pvals, PostPAW_temp_pred$pvals[,c(1:3,13)],
                             by=c("Genotype","density","PostPAW"))

# For some strange reason, left_join has created duplicates of some plots, so need to ensure each row of the dataframe is unique
PostPAW_pred_df <- unique(PostPAW_pred_df)
# Also ensure there are no duplicate rows in the dataframe below
Final_model_pred$pvals <-  unique(Final_model_pred$pvals)

# Now calculate the lack of fit (LOF) effect with respect to Post-flowering plant available water
# i.e. the difference in the contribution of Post-flowering plant available water to grain yield compared to the contribution from all other environment terms
PostPAW_pred_df$LOF <- PostPAW_pred_df$pred_value_PostPAW - PostPAW_pred_df$predicted.value

Sorghum_EC1_adj_df <- left_join(Sorghum_df, PostPAW_pred_df[,c(1:4,16)],
                              by=c('Genotype','Trial','Env', 'density'))

# Finally, calculate the yield values adjusted for all ECs and LoF terms except for Post-flowering plant avaliable water
Sorghum_EC1_adj_df$Diff <- Sorghum_EC1_adj_df$Yld + Sorghum_EC1_adj_df$LOF

# Check that mean of lack of fit values is approximately equal to zero
mean(Sorghum_EC1_adj_df$LOF, na.rm=T)



# 4.3 Obtain predictions for Post-flowering plant avaliable water from the final model ----

# Obtain predictions for the Post-flowering plant avaliable water by genotype interaction effect
PostPAW_pred <- predict( Yield_final_fm, classify='PostPAW:ISW:PreCumRad:PrePAW:ptq:PreFlwEvap:PostMaxT:Genotype:density', 
         levels=list('PostPAW'=aux_parallel_final$PostPAW,
                     'Genotype'=aux_parallel_final$Genotype,
                     'ISW'=rep(0, 1443),
                     'PreCumRad'=rep(0, 1443),
                     'PrePAW'=rep(0, 1443),
                     'ptq'=rep(0, 1443),
                     'PreFlwEvap'=rep(0, 1443),
                     'PostMaxT'=rep(0, 1443),
                     'density'=rep(6.611,1443)),
         present=c("Genotype","PostPAW","density"),
         maxit=1,parallel=T)

# Remove any redundant rows from the data frame of predictions
PostPAW_pred$pvals <- unique(PostPAW_pred$pvals)

# Define a colour-blind palette for the figure
cbbPalette <- c("#0072B2" , "#E69F00", "#56B4E9" , "#CC79A7",  "#D55E00", "#009E73", "#000000")

#  Plot the yield response to Post-flowering plant avaliable water for each genotype
#  Note that this is Figure 2 in the manuscript
EnvCov_C1 <- Sorghum_EC1_adj_df %>% 
  filter(TargetPop%in%c("6","9")) %>% 
  ggplot( . , aes( PostPAW + 73 , Diff )) + # Add the PostPAW mean (see Table 2) in the x-axis value
  geom_ribbon(data=PostPAW_pred$pvals, aes(y=predicted.value, 
              ymin=predicted.value-2*std.error,ymax=predicted.value+2*std.error),
              colour=NA, alpha=0.1) +
  geom_point( size=.75, aes(colour=Trial, shape=TOS)) +
  geom_smooth(data=PostPAW_pred$pvals, aes(y=predicted.value), span=0.5, linewidth=.4) +
  facet_wrap(~Genotype, ncol=4) + 
  xlab("Post-flowering plant avaliable water (mm)") +
  ylab("Yield (t/ha)") +
  theme_bw() +
  theme(  panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"),
          legend.title = element_text(size=14),
          legend.text = element_text(size=11),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=15, angle=90, vjust=0.4),
          axis.title.y = element_text(size=16),
          axis.title.x = element_text(size=16),
          axis.text.y = element_text(size=15),
          strip.text.x = element_text(size=11),
          strip.text.y = element_text(size=11),
          strip.background = element_rect(colour="black", fill="white")) +
  scale_colour_manual(values=cbbPalette, name="Trial") + 
  scale_shape_discrete(name="TOS")
# View the ggplot object 
EnvCov_C1


# 4.4 Obtain predictions for the yield response to established plant population and Post-flowering plant avaliable water simultaneously ----

# Create a data frame consisting of all the unique combinations of established plant population and Post-flowering plant avaliable water that we want to obtain predictions for
PostPAWdensity_pred <- data.frame(PostPAW=rep(seq(-65,150,3), each=64))
PostPAWdensity_pred$density <- rep(seq(1.4,14,0.2), times=72)

# Obtain predictions for the data frame created previously
PostPAWdensity_pred <- predict( Yield_final_fm, classify='PostPAW:ISW:PreCumRad:PrePAW:ptq:PreFlwEvap:PostMaxT:density', 
                                levels=list('PostPAW'=PostPAWdensity_pred$PostPAW,
                                            'ISW'=rep(0,4608),
                                            'PreCumRad'=rep(0,4608),
                                            'PrePAW'=rep(0,4608),
                                            'ptq'=rep(0,4608),
                                            'PreFlwEvap'=rep(0,4608),
                                            'PostMaxT'=rep(0,4608),
                                            'density'=PostPAWdensity_pred$density),
                                present=c("Genotype"), pworkspace="4gb",maxit=1,
                                parallel=T)


# Plot predictions using a two-dimensional contour plot
EC_1_density <-  ggplot( PostPAWdensity_pred$pvals , aes( PostPAW , density )) +
  #geom_raster(  aes(fill=predicted.value)) +
  geom_contour_filled(  aes(z=predicted.value)) +
  xlab("Post-flowering plant avaliable water (mm)") +
  ylab(expression(~Established~plant~population~"(plants/"*m^2*")")) +
  theme_bw() +
  #guides(fill=guide_colourbar(barwidth=1, barheight=20)) +
  theme(  panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"),
          legend.title = element_text(size=12),
          legend.text = element_text(size=12),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=15),
          axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=15),
          axis.text.y = element_text(size=15),
          strip.text.x = element_text(size=11),
          strip.text.y = element_text(size=11)) +
  scale_fill_viridis(name="Yield\n(t/ha)", option="magma", discrete=T)
# Print ggplot object
EC_1_density



# Generate 3d figure using 'plotly' package
plot_3d_EC1 <- plot_ly(PostPAWdensity_pred$pvals, y = ~ density  , x = ~PostPAW,
              z = ~ predicted.value, #type='mesh3d',
              marker = list(color = ~ predicted.value,
                            colorscale = 'Viridis', showscale = TRUE)) %>%
  add_markers() %>%
  layout(scene = list(zaxis = list(title = 'Yield (t/ha)', rangemode="tozero"),
                      xaxis = list(title = 'Post-flowering rain (mm)'),
                      yaxis = list(title = 'Established plant population (plants/m<sup>2</sup>)')))
# Print plotly object
plot_3d_EC1

# Generate 3d figure using the 'plot3d' package

# Define an array to generate a grid of predicted values
yld_pred_mat <- array(PostPAWdensity_pred$pvals$predicted.value, dim=c(64,72))
plantpop_mat <- array(PostPAWdensity_pred$pvals$density, dim=c(64,72))
PostPAW_mat <- array(PostPAWdensity_pred$pvals$PostPAW, dim=c(64,72))

# Filter every 4th row and column in array
yld_pred_mat <- yld_pred_mat %>% as.data.frame() %>% 
  select(everything()[c(TRUE,FALSE, FALSE, FALSE)], last_col()) %>% 
  filter(row_number() %% 4 == 1) %>% 
  as.matrix()

plantpop_mat <- plantpop_mat %>% as.data.frame() %>% 
  select(everything()[c(TRUE,FALSE, FALSE, FALSE)], last_col()) %>% 
  filter(row_number() %% 4 == 1) %>% 
  as.matrix()

PostPAW_mat <- PostPAW_mat %>% as.data.frame() %>% 
  select(everything()[c(TRUE,FALSE, FALSE, FALSE)], last_col()) %>% 
  filter(row_number() %% 4 == 1) %>% 
  as.matrix()

# Plot the 3d image.
# This is identical to Figure 3 in the manuscript
scatter3D(x= Sorghum_EC1_adj_df$PostPAW + 73, 
          y= Sorghum_EC1_adj_df$density, 
          z= Sorghum_EC1_adj_df$Diff,
          theta=315, phi=20, pch=19, ticktype="detailed",
          clab="Yield (t/ha)",
          cex=0.5, 
          col=ramp.col(plasma(n=100)),
          xlab='Post-flowering plant avaliable water (mm)',
          ylab='Established plant population (plants/m2)',
          zlab='Yield (t/ha)',
          surf=list(x=PostPAW_mat + 73,
                    y=plantpop_mat,
                    z=yld_pred_mat, 
                    facets=NA, fit=Sorghum_EC1_adj_df$Diff))


# 5.1 Obtain predictions for the yield response to established plant population using the 7 important ECs identified as a surrogate for the environment term ----

# Calculate predicted values for the GxECxM interaction using all the ECs as a surrogate for the Env term
Final_pred <- predict( Yield_final_fm, classify='Genotype:PostPAW:PrePAW:ISW:PreCumRad:PreCumRad:ptq:PostMaxT:PreFlwEvap:density', 
                      levels=list('Genotype'=aux_parallel_final$Genotype,
                                  'PostPAW'=aux_parallel_final$PostPAW,
                                  'ISW'=aux_parallel_final$ISW,
                                  'PreCumRad'=aux_parallel_final$PreCumRad,
                                  'PostMaxT'=aux_parallel_final$PostMaxT,
                                  'PrePAW'=aux_parallel_final$PrePAW,
                                  'ptq'=aux_parallel_final$ptq,
                                  'PreFlwEvap'=aux_parallel_final$PreFlwEvap,
                                  'PreCumRad'=aux_parallel_final$PreCumRad,
                                  'density'=aux_parallel_final$density),
                      parallel=T, pworkspace="4gb", maxit=1)

# Round predicted values for each EC to two decimal places to merge dataframes smoothly
Final_pred$pvals$PostPAW <- round(Final_pred$pvals$PostPAW,2)
aux_parallel_final$PostPAW <-  round(aux_parallel_final$PostPAW,2)
# Check that the predictions of ECs in the two data frames match 
#Final_pred$pvals$PostPAW == aux_parallel_final$PostPAW

Final_pred$pvals$ISW <- round(Final_pred$pvals$ISW,2)
aux_parallel_final$ISW <-  round(aux_parallel_final$ISW,2)
#Final_pred$pvals$ISW == aux_parallel_final$ISW

Final_pred$pvals$PreCumRad <- round(Final_pred$pvals$PreCumRad,2)
aux_parallel_final$PreCumRad <-  round(aux_parallel_final$PreCumRad,2)
#Final_pred$pvals$PreCumRad == aux_parallel_final$PreCumRad

Final_pred$pvals$PrePAW <- round(Final_pred$pvals$PrePAW,2)
aux_parallel_final$PrePAW <-  round(aux_parallel_final$PrePAW,2)
#Final_pred$pvals$PrePAW == aux_parallel_final$PrePAW

Final_pred$pvals$ptq <- round(Final_pred$pvals$ptq,2)
aux_parallel_final$ptq <-  round(aux_parallel_final$ptq,2)
#Final_pred$pvals$ptq == aux_parallel_final$ptq

Final_pred$pvals$PreFlwEvap <- round(Final_pred$pvals$PreFlwEvap,2)
aux_parallel_final$PreFlwEvap <-  round(aux_parallel_final$PreFlwEvap,2)
#Final_pred$pvals$PreFlwEvap == aux_parallel_final$PreFlwEvap

Final_pred$pvals$PostMaxT <- round(Final_pred$pvals$PostMaxT,2)
aux_parallel_final$PostMaxT <-  round(aux_parallel_final$PostMaxT,2)
#Final_pred$pvals$PostMaxT == aux_parallel_final$PostMaxT

# Round predicted values for established plant population to two decimal places to merge dataframes smoothly
Final_pred$pvals$density <- round(Final_pred$pvals$density,2)
aux_parallel_final$density <-  round(aux_parallel_final$density,2)
#Final_pred$pvals$density == aux_parallel_final$density

# Merge predicted values with the auxiliary data frame
Final_pred$pvals <- left_join(Final_pred$pvals, aux_parallel_final, by=c("Genotype","density","PostPAW",
                          "ISW","PreCumRad","PrePAW","ptq","PreFlwEvap","PostMaxT"))


# Re-order genotype levels to be in chronological order
Baseline_GxExM_pred$pvals$Genotype <- factor(Baseline_GxExM_pred$pvals$Genotype, levels=c(
                                    "GenoA","GenoB","GenoC","GenoD",
                                    "GenoE","GenoF","GenoG","GenoH"))
SorghumMET_Adj.df$Genotype <- factor(SorghumMET_Adj.df$Genotype, levels=c(
                                    "GenoA","GenoB","GenoC","GenoD",
                                    "GenoE","GenoF","GenoG","GenoH"))
Final_pred$pvals$Genotype <- factor(Final_pred$pvals$Genotype, levels=c(
                              "GenoA","GenoB","GenoC","GenoD",
                              "GenoE","GenoF","GenoG","GenoH"))

# Plot the predicted values (tested environment) for the yield response to plant population for each GxE combination
# Red = Baseline GxM predictions
# Blue = Baseline GxExM predictions
# Green = Final model GxECxM predictions
print(ggplot( Baseline_GxExM_pred$pvals , aes( density , predicted.value )) +
        #geom_ribbon(data=Final_pred$pvals, aes(ymin=predicted.value-2*std.error, #95% prediction interval
        #               ymax=predicted.value+2*std.error, group=Env),colour="Green4", alpha=0.2) +
        geom_line(colour="Blue", size=0.5) +
        geom_line(data=Baseline_GxM_pred$pvals, colour="Red", size=0.4) +
        geom_point(data = SorghumMET_Adj.df , aes(x = density , y = Yield_Adj) ,
                   stat = "identity" , size=0.8) +
        geom_line(data=Final_pred$pvals, colour="Green4", size=0.55) +
        facet_grid(Genotype~Env) + 
        xlab(expression(~Established~plant~population~"(plants/"*m^2*")")) +
        ylab("Yield (t/ha)") +
        xlim(c(2,10.5)) +
        theme_bw() +
        theme(  panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"),
                legend.title = element_text(size=17),
                legend.text = element_text(size=15),
                axis.line.y = element_line(colour = "black"),
                axis.text.x = element_text(size=18, angle=90, vjust=0.4), 
                axis.title.y = element_text(size=20),
                axis.title.x = element_text(size=20),
                axis.text.y = element_text(size=18),
                strip.text.x = element_text(size=14), 
                strip.text.y = element_text(size=14),
                strip.background = element_rect(colour="black", fill="white")) +
        scale_colour_manual(values=cbbPalette))


# 5.2 Perform leave-one-trial-out cross-validation for the final model ----

# Define a term to store the cross_validation predictions at each iteration
CrossValid_final_pred <- c()
for (j in c(1:6)){
  
  # Subset one trial out of the data frame
  Sorghum_subset_df <- Sorghum_df[Sorghum_df$Trial!=levels(Sorghum_df$Trial)[j],]
  
  # Re-factorise each of the relevant terms to update the number of factor levels for each term
  Sorghum_subset_df$Trial <- factor(Sorghum_subset_df$Trial)
  Sorghum_subset_df$Env <- factor(Sorghum_subset_df$Env)
  Sorghum_subset_df$ResidualTOS <- factor(Sorghum_subset_df$ResidualTOS)
  Sorghum_subset_df$Genotype <- factor(Sorghum_subset_df$Genotype)
  # Run a slightly different version of the model based on which trial has been removed in the cross validation procedure
  if(j==1){
    # Note that the knot points are reduced back to 6 for cross validation to significantly speed up the cross validation process
    Temp_subset_fm <- asreml( Yld ~ Genotype + density + Genotype:density +
                                PostPAW + Genotype:PostPAW +
                                PrePAW + density:PrePAW +
                                ISW + density:ISW +
                                PreCumRad + density:PreCumRad +
                                ptq + Genotype:ptq +
                                PreFlwEvap + density:PreFlwEvap +
                                PostMaxT + density:PostMaxT,
                              random =~ at(Trial):Rep  + at(Trial):MainPlot +
                                at(Trial,c(1:3)):SubPlot +
                                #at(Trial,'Breeza 1'):Column +
                                Trial + Env +
                                spl(density, k=6) + spl(density, k=6):Genotype +
                                str(~Trial:Genotype + Trial:Genotype:density,
                                    ~corh(2):id(40)) +
                                str(~Env:Genotype + Env:Genotype:density,
                                    ~corh(2):id(112)) +
                                spl(PostPAW, k=6) +
                                spl(PostPAW, k=6):Genotype +
                                spl(density , k=6):ISW +
                                spl(ptq, k=6) + spl(ptq, k=6):Genotype +
                                spl(density , k=6):PreFlwEvap +
                                spl(PostMaxT, k=6) + density:spl(PostMaxT, k=6),
                              residual=~ dsum(~units|ResidualTOS),
                              data = Sorghum_subset_df,
                              na.action=na.method(x='include'),
                              maxit=100, workspace="1gb")
  } else if(j==3) {
    Temp_subset_fm <- asreml( Yld ~ Genotype + density + Genotype:density +
                                PostPAW + Genotype:PostPAW +
                                PrePAW + density:PrePAW +
                                ISW + density:ISW +
                                PreCumRad + density:PreCumRad +
                                ptq + Genotype:ptq +
                                PreFlwEvap + density:PreFlwEvap +
                                PostMaxT + density:PostMaxT,
                              random =~ at(Trial):Rep  + at(Trial):MainPlot +
                                at(Trial,c(1:3)):SubPlot +
                                at(Trial,'Breeza 1'):Column +
                                Trial + Env +
                                spl(density, k=6) + spl(density, k=6):Genotype +
                                str(~Trial:Genotype + Trial:Genotype:density,
                                    ~corh(2):id(40)) +
                                str(~Env:Genotype + Env:Genotype:density,
                                    ~corh(2):id(120)) +
                                spl(PostPAW, k=6) +
                                spl(PostPAW, k=6):Genotype +
                                spl(density , k=6):ISW +
                                spl(ptq, k=6) + spl(ptq, k=6):Genotype +
                                spl(density , k=6):PreFlwEvap +
                                spl(PostMaxT, k=6) + density:spl(PostMaxT, k=6),
                              residual=~ dsum(~units|ResidualTOS),
                              data = Sorghum_subset_df,
                              na.action=na.method(x='include'),
                              maxit=100, workspace="1gb")
  } else {
    Temp_subset_fm <- asreml( Yld ~ Genotype + density + Genotype:density +
                                PostPAW + Genotype:PostPAW +
                                PrePAW + density:PrePAW +
                                ISW + density:ISW +
                                PreCumRad + density:PreCumRad +
                                ptq + Genotype:ptq +
                                PreFlwEvap + density:PreFlwEvap +
                                PostMaxT + density:PostMaxT,
                              random =~ at(Trial):Rep  + at(Trial):MainPlot +
                                at(Trial,c(1:3)):SubPlot +
                                at(Trial,'Breeza 1'):Column +
                                Trial + Env +
                                spl(density, k=6) + spl(density, k=6):Genotype +
                                str(~Trial:Genotype + Trial:Genotype:density,
                                    ~corh(2):id(40)) +
                                str(~Env:Genotype + Env:Genotype:density,
                                    ~corh(2):id(112)) +
                                spl(PostPAW, k=6) +
                                spl(PostPAW, k=6):Genotype +
                                spl(density , k=6):ISW +
                                spl(ptq, k=6) + spl(ptq, k=6):Genotype +
                                spl(density , k=6):PreFlwEvap +
                                spl(PostMaxT, k=6) + density:spl(PostMaxT, k=6),
                              residual=~ dsum(~units|ResidualTOS),
                              data = Sorghum_subset_df,
                              na.action=na.method(x='include'),
                              maxit=100, workspace="1gb")
  }
  
  # Create data frame with the combination of all explanatory variables to obtain predictions for (including those in the quasi-untested environment)
  aux_parallel_temp <- unique(Sorghum_df[,c( 'Env', 'Trial', 'Genotype', 'density', 'PostPAW', 'ISW', 'PreCumRad',
                                             'PrePAW','ptq','PreFlwEvap','PostMaxT')])
  # Re-order the dataframe
  aux_parallel_temp <- aux_parallel_temp[order( aux_parallel_temp$Env, aux_parallel_temp$Genotype,
                                                aux_parallel_temp$density),]
  
  # Obtain predictions for all terms using the candidate ec as a surrogate for the environment term in the classify set
  Temp_subset_final_pred <- predict(Temp_subset_fm, classify='Genotype:density:PostPAW:ISW:PreCumRad:PrePAW:ptq:PreFlwEvap:PostMaxT', 
                                levels=list('Genotype'=aux_parallel_temp$Genotype,
                                            'PostPAW'=aux_parallel_temp$PostPAW,
                                            'ISW'=aux_parallel_temp$ISW,
                                            'PreCumRad'=aux_parallel_temp$PreCumRad,
                                            'PrePAW'=aux_parallel_temp$PrePAW,
                                            'ptq'=aux_parallel_temp$ptq,
                                            'PreFlwEvap'=aux_parallel_temp$PreFlwEvap,
                                            'PostMaxT'=aux_parallel_temp$PostMaxT,
                                            'density'=aux_parallel_temp$density),
                                parallel=T, pworkspace="2Gb", maxit=30)
  
  # Round numeric terms for the data frame of (i) predictions and (ii) combination of variables to predict for
  # to ensure the merging of the data-frames is smooth
  Temp_subset_final_pred$pvals$PostPAW <- round(Temp_subset_final_pred$pvals$PostPAW,2)
  aux_parallel_temp$PostPAW <-  round(aux_parallel_temp$PostPAW,2)
  Temp_subset_final_pred$pvals$ISW <- round(Temp_subset_final_pred$pvals$ISW,2)
  aux_parallel_temp$ISW <-  round(aux_parallel_temp$ISW,2)
  Temp_subset_final_pred$pvals$PreCumRad <- round(Temp_subset_final_pred$pvals$PreCumRad,2)
  aux_parallel_temp$PreCumRad <-  round(aux_parallel_temp$PreCumRad,2)
  Temp_subset_final_pred$pvals$PrePAW<- round(Temp_subset_final_pred$pvals$PrePAW,2)
  aux_parallel_temp$PrePAW <-  round(aux_parallel_temp$PrePAW,2)
  Temp_subset_final_pred$pvals$ptq <- round(Temp_subset_final_pred$pvals$ptq,2)
  aux_parallel_temp$ptq <-  round(aux_parallel_temp$ptq,2)
  Temp_subset_final_pred$pvals$PreFlwEvap <- round(Temp_subset_final_pred$pvals$PreFlwEvap,2)
  aux_parallel_temp$PreFlwEvap <-  round(aux_parallel_temp$PreFlwEvap,2)
  Temp_subset_final_pred$pvals$PostMaxT <- round(Temp_subset_final_pred$pvals$PostMaxT,2)
  aux_parallel_temp$PostMaxT <-  round(aux_parallel_temp$PostMaxT,2)
  Temp_subset_final_pred$pvals$density <- round(Temp_subset_final_pred$pvals$density,2)
  aux_parallel_temp$density <-  round(aux_parallel_temp$density,2)
  
  # Merge data frame of (i) predictions and (ii) combination of variables to predict for
  Temp_subset_final_pred$pvals <- left_join(Temp_subset_final_pred$pvals, aux_parallel_temp, 
                                            by=c("Genotype","density","PostPAW","ISW","PreCumRad","PrePAW","ptq","PreFlwEvap","PostMaxT")) 
  
  # Subset so that the table of predictions only contains predictions that were from a quasi-untested environment
  Temp_subset_final_pred <- Temp_subset_final_pred$pvals[Temp_subset_final_pred$pvals$Trial==levels(Sorghum_df$Trial)[j],]
  # Merge the current set of predictions from an untested environment with the previous set(s) of predictions from an untested environment
  CrossValid_final_pred <-  bind_rows(CrossValid_final_pred, Temp_subset_final_pred)
}


# 5.3 Output the predictions for the yield response to established plant population ----
# for a subset of 2 genotypes and 4 environments

# Subset for the genotypes and environments included in Figure 4 in the manuscript
SorghumMET_Adj.df <- SorghumMET_Adj.df[SorghumMET_Adj.df$Genotype%in%c("GenoC","GenoG"),]
SorghumMET_Adj.df <- SorghumMET_Adj.df[SorghumMET_Adj.df$Env%in%c("Moree 08/08/2018","Moree 27/09/2018", 
                                                                  "Surat 08/08/2018","Surat 24/01/2019"),]
# Re-factorise to update levels of Genotype and Env terms
SorghumMET_Adj.df$Genotype <- factor(SorghumMET_Adj.df$Genotype)
SorghumMET_Adj.df$Env <- factor(SorghumMET_Adj.df$Env)
# Repeat for final table of predicted values (tested environment)
Final_pred$pvals <- Final_pred$pvals[Final_pred$pvals$Genotype%in%c("GenoC","GenoG"),]
Final_pred$pvals <- Final_pred$pvals[Final_pred$pvals$Env%in%c("Moree 08/08/2018","Moree 27/09/2018", 
                                                               "Surat 08/08/2018","Surat 24/01/2019"),]
Final_pred$pvals$Genotype <- factor(Final_pred$pvals$Genotype)
Final_pred$pvals$Env <- factor(Final_pred$pvals$Env)

# Repeat for final table of predicted values (untested environment)
CrossValid_final_subset_pred <- CrossValid_final_pred[CrossValid_final_pred$Genotype%in%c("GenoC","GenoG"),]
CrossValid_final_subset_pred <- CrossValid_final_subset_pred[CrossValid_final_subset_pred$Env%in%c("Moree 08/08/2018","Moree 27/09/2018", 
                                                                                                "Surat 08/08/2018","Surat 24/01/2019"),]
CrossValid_final_subset_pred$Genotype <- factor(CrossValid_final_subset_pred$Genotype)
CrossValid_final_subset_pred$Env <- factor(CrossValid_final_subset_pred$Env)

# Plot below is Figure 4 in the manuscript 
GxExM_final <- ggplot( Final_pred$pvals , aes( density , predicted.value )) +
  geom_ribbon(aes(ymin=predicted.value-2*std.error,
                  ymax=predicted.value+2*std.error, group=Env),fill="Red",
              alpha=0.1, linetype=0) +
  geom_ribbon(data=CrossValid_final_subset_pred, aes(ymin=predicted.value-2*std.error,
                              ymax=predicted.value+2*std.error, group=Env),fill="Blue",
              alpha=0.1, linetype=0) +
  geom_point(data = SorghumMET_Adj.df , aes(x = density , y = Yield_Adj) ,
             stat = "identity" , size=0.8) +
  geom_smooth(colour="Red", size=0.4) +
  geom_smooth(data=CrossValid_final_subset_pred, colour="Blue", size=0.4) +
  facet_grid(Genotype~Env) + 
  xlab(expression(~Established~plant~population~"(plants/"*m^2*")")) +
  ylab("Yield (t/ha)") +
  xlim(c(2,13)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"),
        plot.title = element_text(size=16, face="bold"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=15),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=15, angle=90, vjust=0.4), 
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=15),
        strip.text.x = element_text(size=10), 
        strip.text.y = element_text(size=10),
        strip.background = element_rect(colour="black", fill="white")) +
  scale_colour_manual(values=cbbPalette)
# Display the figure
GxExM_final












