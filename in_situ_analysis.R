# Random Forest
# Script by Lucas Ustick
# Compare genes and satellite measurements using random forest

#load data
sat_data <- readRDS("sat_matchup_data.RData")
all_data <- read.csv("in_situ_data.csv")

#Find best fit for RF model ####
library(randomForest)

nut_fits <- matrix(ncol = 15, nrow = 6)
samps <- matrix(ncol = 15, nrow = 6)

#Create RF_model for each fit of the data
for (time in 1:15) {
  for (space in 1:6) {
    nut_lim_sat <- sat_data[[space]][,time]
    samp <- !is.na(nut_lim_sat) 
    nut_sat <- nut_lim_sat[samp]
    
    #extract the data
    p_lim <- all_data$Omega_P_High[samp]
    n_lim <- all_data$Omega_N_High[samp]
    nitri <- all_data$nutricline[samp]
    omega_nitri <- cbind(p_lim,n_lim,nitri)
    
    model_rf <- randomForest(omega_nitri,nut_sat)
    
    #calculate the variance explained (r2)
    n <- length(model_rf$y)
    rsq <- 1 - model_rf$mse/(var(model_rf$y) * (n - 1)/n)
    r <- round(100 * rsq[length(rsq)], digits = 2)
    
    nut_fits[space,time] <- r
  }
}


library(reshape)
#df has all the R2 values for each model
df <- melt(nut_fits)
colnames(df) <- c("x", "y", "R2")
  
# Df_2 contains all p(t,s) scores
df_2 <- df

#check previous best fit, calculate p(t,s) score (eq 1&2 in the methods)
for (i in 1:90) { #improvement from smaller fits
  tm = df[i,2] #1-15
  sp = df[i,1] #1-6
  if (tm == 1 & sp == 1) {#corner 1
    df_2[i,4] <- 0
  } else if (tm == 1) {#side 1
    df_2[i,4] <- (df_2[i,3]-max(nut_fits[1:sp-1,tm]))/max(nut_fits[1:sp-1,tm])
  } else if (sp == 1) {#side 2
    df_2[i,4] <- (df_2[i,3]-max(nut_fits[sp,1:tm-1]))/max(nut_fits[sp,1:tm-1])
  } else  {
    df_2[i,4] <- (df_2[i,3]-max(c(nut_fits[1:sp-1,1:tm-1]),nut_fits[sp,1:tm-1],nut_fits[1:sp-1,tm]))/max(c(nut_fits[1:sp-1,1:tm-1]),nut_fits[sp,1:tm-1],nut_fits[1:sp-1,tm])
  }
  
} 



#Run the best model####
library(randomForest)
library(pdp)

space <- 5 
time <- 6

nut_lim_sat <- sat_data[[space]][,time]
samp <- !is.na(nut_lim_sat) 
nut_sat <- nut_lim_sat[samp]
#extract the data
p_lim <- all_data$Omega_P_High[samp]
n_lim <- all_data$Omega_N_High[samp]
nitri <- all_data$nutricline[samp]
omega_nitri <- cbind(p_lim,n_lim,nitri)

# run the random forest
best_rf <- randomForest(omega_nitri,nut_sat,ntree=500)
  
# calculate the R2 for the random forest
n <- length(best_rf$y)
rsq <- 1 - best_rf$mse/(var(best_rf$y) * (n - 1)/n)
r_squared_final <- rsq[length(rsq)] #get R2

#get model prediction
rf_pred <- predict(best_rf,omega_nitri) 
  
#calculate partials
partial_n_range <- seq(-1, 1.8, by = 0.1)
partial_p_range <- seq(-.6, 3, by = 0.1)
partial_nitri_range <- seq(5,215, by = 5)
ni_n_range <- expand.grid(partial_nitri_range, partial_n_range)
colnames(ni_n_range) <- c("nitri","n_lim")
n_p_range <- expand.grid(partial_n_range, partial_p_range)
colnames(n_p_range) <- c("n_lim","p_lim")

ni_n_partial_final <- partial(best_rf,pred.var = c("nitri","n_lim"),pred.grid=ni_n_range, chull = TRUE)
n_p_partial_final <- partial(best_rf,pred.var = c("n_lim","p_lim"),pred.grid=n_p_range, chull = TRUE)

