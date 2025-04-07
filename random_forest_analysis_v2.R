# Script by Lucas Ustick
# Compare genes and satellite measurements
# 01.00 Find best fit for RF model ####
library(randomForest)

all_data_mod <- read.csv("data/data_insitu_v2_25-04.csv",header = TRUE)

sat_data_1p <- read.csv("data/sat_data_25-04/genomic_OCmatchups_v3_medians_1p.csv")
sat_data_3p <- read.csv("data/sat_data_25-04/genomic_OCmatchups_v3_medians_3p.csv")
sat_data_5p <- read.csv("data/sat_data_25-04/genomic_OCmatchups_v3_medians_5p.csv")
sat_data_1d <- read.csv("data/sat_data_25-04/genomic_OCmatchups_v3_medians_1d.csv")
sat_data_2d <- read.csv("data/sat_data_25-04/genomic_OCmatchups_v3_medians_2d.csv")
sat_data_5d <- read.csv("data/sat_data_25-04/genomic_OCmatchups_v3_medians_5d.csv")
sat_data <- list(sat_data_1d,sat_data_3p,sat_data_5p,sat_data_1d,sat_data_2d,sat_data_5d)

nut_fits <- matrix(ncol = 15, nrow = 6)
samps <- matrix(ncol = 15, nrow = 6)

#Create RF_model for each fit of the data
for (time in 1:15) {
  for (space in 1:6) {
    nut_lim_sat <- sat_data[[space]][,time]
    samp <- !is.na(nut_lim_sat)
    samp_count <- sum(samp)
    nut_sat <- nut_lim_sat[samp]
    #extract the data
    p_lim_h <- all_data_mod$Omega_P_High[samp]
    n_lim_h <- all_data_mod$Omega_N_High[samp]
    n_lim_m <- all_data_mod$Omega_N_Med[samp]
    fe_lim_h <- all_data_mod$Omega_Fe_High[samp]
    fe_lim_m <- all_data_mod$Omega_Fe_Med[samp]
    nitri <- all_data_mod$insitu_and_woa_nutricline[samp]
    
    omegas_nitri <- cbind(nitri,p_lim_h,n_lim_h,n_lim_m,fe_lim_h,fe_lim_m)
    
    set.seed(123)
    model_rf <- randomForest(omegas_nitri,nut_sat)
    
    #calculate the variance explained (r2)
    n <- length(model_rf$y)
    rsq <- 1 - model_rf$mse/(var(model_rf$y) * (n - 1)/n)
    r <- round(100 * rsq[length(rsq)], digits = 2)
    print("hi")
    nut_fits[space,time] <- r
    samps[space,time] <- samp_count
    #break
  }
  #break
}

#make a figure of the best fits
library(reshape)
#make table into form ggplot likes
df <- melt(nut_fits)
count_df <- melt(samps)
colnames(df) <- c("x", "y", "R2")
colnames(count_df) <- c("x1", "y1", "counts")
count_r_df <- cbind(df,count_df)

model_1 <- lm(R2 ~ counts, data=count_r_df)
summary(model_1)
r_squared <- summary(model_1)$r.squared
p_value <- summary(model_1)$coefficients[2, 4]  # p-value for the slope

library(ggplot2)
plot1 <- ggplot(df, aes(x = x, y = y, fill = R2)) +
  geom_tile() +
  geom_text(aes(label = round(R2, digits=2)), color = "black") +
  scale_x_discrete(limits=1:6,labels=c("1P","3x3P","5x5P","1x1D","2x2D","5x5D")) +
  scale_y_discrete(limits=1:15,labels=c("±2D","±4D","±6D","±8D","±10D","±12D"
                                        ,"±14D","±16D","±18D","±20D","±22D"
                                        ,"±24D","±26D","±28D","±30D")) +
  xlab("Space") +
  ylab("Time") +
  ggtitle("RF Nutricline+Genes")+
  scale_fill_distiller(palette = "RdYlBu")

plot1

# 01.01 Nutricline vs. Full model ####
library(randomForest)
library(pdp)

space <- 5 
time <- 10

nut_lim_sat <- sat_data[[space]][,time]
samp <- !is.na(nut_lim_sat)
samp_count <- sum(samp)
nut_sat <- nut_lim_sat[samp]
#extract the data
p_lim_h <- all_data_mod$Omega_P_High[samp]
p_lim_m <- all_data_mod$Omega_P_Med[samp]
n_lim_h <- all_data_mod$Omega_N_High[samp]
n_lim_m <- all_data_mod$Omega_N_Med[samp]
fe_lim_h <- all_data_mod$Omega_Fe_High[samp]
fe_lim_m <- all_data_mod$Omega_Fe_Med[samp]
nitri <- all_data_mod$insitu_and_woa_nutricline[samp]

omegas_nitri <- cbind(nitri,p_lim_h,n_lim_h,n_lim_m,fe_lim_h,fe_lim_m)

set.seed(123)
nutricline_rf <- randomForest(as.data.frame(nitri),nut_sat,ntree=500)
nutricline_rf

set.seed(123)
nut_gen_rf <- randomForest(omegas_nitri,nut_sat)
nut_gen_rf

# Extract predicted values
all_data_mod_save <- all_data_mod
all_data_mod_save$nut_lim_sat[samp] <- nut_sat
all_data_mod_save$nutricline_rf[samp] <- predict(nutricline_rf)
all_data_mod_save$nut_gene_rf[samp] <- predict(nut_gen_rf)

all_data_mod_save$nutricline_rf_residual[samp] <- all_data_mod_save$nut_lim_sat[samp]-all_data_mod_save$nutricline_rf[samp]
all_data_mod_save$nut_gene_rf_residual[samp] <- all_data_mod_save$nut_lim_sat[samp]-all_data_mod_save$nut_gene_rf[samp]

# calculate the R2 for the random forest
n <- length(nutricline_rf$y)
rsq <- 1 - nutricline_rf$mse/(var(nutricline_rf$y) * (n - 1)/n)
r2_nutricline_rf <- rsq[length(rsq)] #get R2 for each bootstrap
# calculate the R2 for the random forest
n <- length(nut_gen_rf$y)
rsq <- 1 - nut_gen_rf$mse/(var(nut_gen_rf$y) * (n - 1)/n)
r2_nut_gen_rf <- rsq[length(rsq)] #get R2 for each bootstrap

rf_df <- data.frame(r2 = c(r2_nutricline_rf, r2_nut_gen_rf))
rf_df$model <- c("nutricline","nutricline+genes")

library(ggplot2)

ggplot(rf_df,aes(y=r2,x=model))+
  geom_bar(stat="identity")+
  xlab("Model")+
  ylab("RandomForest R^2")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))