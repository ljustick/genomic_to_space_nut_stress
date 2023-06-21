# Script by Lucas Ustick
# Compare genes and satellite measurements using random forest

#load data
sat_data <- readRDS("sat_matchup_data.RData")
all_data <- read.csv("in_situ_data.csv")

# Find best fit for RF model ####
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

# Vsualize
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
  ggtitle("R (t,s)")+
  scale_fill_distiller(palette = "RdYlBu")

library(ggplot2)
plot2 <- ggplot(df_2, aes(x = x, y = y, fill = V4)) +
  geom_tile() +
  geom_text(aes(label = round(V4, digits=2)), color = "black") +
  scale_x_discrete(limits=1:6,labels=c("1P","3x3P","5x5P","1x1D","2x2D","5x5D")) +
  scale_y_discrete(limits=1:15,labels=c("±2D","±4D","±6D","±8D","±10D","±12D"
                                        ,"±14D","±16D","±18D","±20D","±22D"
                                        ,"±24D","±26D","±28D","±30D")) +
  xlab("Space") +
  ylab("Time") +
  ggtitle("P (t,s)")+
  scale_fill_distiller(palette = "RdYlBu")

library(ggpubr) #package for combining ggplots
ggarrange(plot1,plot2,
          labels = c("A)", "B)"),
          ncol = 2, nrow = 1)

# Run the best model####
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

# comparison of in situ with remote sensing ####

# prep data
lat <- all_data$LatitudeDegN[samp]
lon <- all_data$LongitudeDegE[samp]
samples <- all_data$Sample[samp]
p_conc <- all_data$phosphate[samp]
n_conc <- all_data$nitrate[samp]
SST <- all_data$SST_underway[samp]
diff <- (nut_sat-rf_pred)
model_data <- as.data.frame(cbind(omega_nitri,nut_sat,rf_pred,n_conc,p_conc,SST,diff,lat,lon,samples))

model_data <- transform(model_data,
                        p_lim=as.numeric(p_lim),
                        n_lim=as.numeric(n_lim),
                        nitri=as.numeric(nitri),
                        nut_sat=as.numeric(nut_sat),
                        rf_pred=as.numeric(rf_pred),
                        n_conc=as.numeric(n_conc),
                        p_conc=as.numeric(p_conc),
                        SST=as.numeric(SST),
                        diff=as.numeric(diff),
                        lat=as.numeric(lat),
                        lon=as.numeric(lon),
                        samples=samples
                        )

# pearson correlations
plim_test <- cor.test(model_data$nut_sat,model_data$p_lim ,method = "pearson")
nlim_test <- cor.test(model_data$nut_sat,model_data$n_lim ,method = "pearson")
sst_test <- cor.test(model_data$nut_sat,model_data$SST ,method = "pearson")
nitri_test <- cor.test(model_data$nut_sat,model_data$nitri ,method = "pearson")
pconc_test <- cor.test(model_data$nut_sat,model_data$p_conc ,method = "pearson")
nconc_test <- cor.test(model_data$nut_sat,model_data$n_conc ,method = "pearson")

Category<-c("ΩP","ΩN","SST","Nutricline","PO4","NO3")
Pearson_corr<-c(plim_test[4][[1]][[1]],nlim_test[4][[1]][[1]],sst_test[4][[1]][[1]],nitri_test[4][[1]][[1]],pconc_test[4][[1]][[1]],nconc_test[4][[1]][[1]])
pearson_pvals<-c(plim_test[3][[1]][[1]],nlim_test[3][[1]][[1]],sst_test[3][[1]][[1]],nitri_test[3][[1]][[1]],pconc_test[3][[1]][[1]],nconc_test[3][[1]][[1]])
pearson_df<-data.frame(cbind(Category,Pearson_corr,c(1:6),pearson_pvals))
pearson_df <- transform(pearson_df,
                        Pearson_corr=as.numeric(Pearson_corr)
)

# Visualize
library(ggplot2)
ggplot(pearson_df,aes(y=Pearson_corr,x=V3))+
  geom_bar(stat="identity")+
  scale_x_discrete(labels=Category)+
  ylab("Pearson Correlation Coefficient")+
  xlab("")+
  theme_bw()

# Transect plots ####
library(ggplot2)

#AMT
samps_fig=12:70
name="AMT"

theme2 <-theme(panel.background = element_blank() #remove background
               ,panel.grid.major = element_blank()
               ,panel.grid.minor = element_blank()
               ,strip.background=element_blank()
               ,axis.text.x=element_text(colour="black") #make all text black
               ,axis.text.y=element_text(colour="black") #make all text black
               ,axis.ticks=element_line(colour="black") #make all text black
               ,axis.line.x=element_line(size=.5)
               ,axis.line.y=element_line(size=.5)
)

point_size=1
line_size=.5

cols1 <- c("SAT"="#FE6100","RF"="#785EF0")
fig1 <- ggplot(model_data[samps_fig,], aes(x=lat)) +
  labs(y="Θ'",x="",color="data")+ #remove name from labs
  #RF prediction
  geom_point( aes(y=rf_pred,color="RF"),size=point_size,shape=15)+
  geom_smooth(aes(y=rf_pred), span = 0.3, se = FALSE, color="#785EF0",size=line_size)+
  #satelite data
  geom_point( aes(y=nut_sat,color="SAT"),size=point_size,shape=15)+
  geom_smooth(aes(y=nut_sat), span = 0.3, se = FALSE, color="#FE6100",size=line_size)+
  #add in custom legend
  scale_colour_manual(name="data",values=cols1)+
  theme2+
  theme(legend.position = "none")


cols2 <- c("P"="#648FFF","N"="#FFB000", "Nutricline"="#000000")
coeff2 <- 75
fig2 <- ggplot(model_data[samps_fig,], aes(x=lat)) +
  labs(x="")+
  #P lim
  geom_point( aes(y=p_lim,color="P"),size=point_size)+
  geom_smooth(aes(y=p_lim), span = 0.3, se = FALSE, color="#648FFF",size=line_size)+
  #N Lim
  geom_point( aes(y=n_lim,color="N"),size=point_size)+
  geom_smooth(aes(y=n_lim), span = 0.3, se = FALSE, color="#FFB000",size=line_size)+
  #Nutricline
  geom_point( aes(y=(nitri/coeff2)-1,color="Nutricline"),size=point_size,shape=6)+
  geom_smooth(aes(y=(nitri/coeff2)-1), span = 0.3, se = FALSE, color="#000000",size=line_size,linetype="dashed")+
  # Add in second axis for Nitricline
  scale_y_continuous(name = "Ω",sec.axis = sec_axis(~.*coeff2+coeff2, name="Nutricline Depth (M)"),limits=c(-1,3) )+
  #add in custom legend
  scale_colour_manual(name="Ω",values=cols2)+
  theme2+
  theme(legend.position = "none")

cols3 <- c("P"="#648FFF","N"="#FFB000", "Temp"="#000000")
coeff3 <- 3.2
fig3 <- ggplot(model_data[samps_fig,], aes(x=lat)) +
  labs(x="Latitude Degrees N",y="Concentration")+
  #Phosphate
  geom_point( aes(y=p_conc,color="P"),size=point_size,shape=5)+
  geom_smooth(aes(y=p_conc), span = 0.3, se = FALSE, color="#648FFF",size=line_size)+
  #Nitrate
  geom_point( aes(y=n_conc,color="N"),size=point_size,shape=5)+
  geom_smooth(aes(y=n_conc), span = 0.3, se = FALSE, color="#FFB000",size=line_size)+
  #Temperature
  geom_point( aes(y=(SST/coeff3),color="Temp"),size=point_size,shape=8)+
  geom_smooth(aes(y=(SST/coeff3)), span = 0.3, se = FALSE, color="#000000",size=line_size,linetype="dashed")+
  # Add in second axis for Temperature
  scale_y_continuous(name = "Concentration",sec.axis = sec_axis(~.*coeff3, name="Temperature (Deg C)"),limits=c(0,10) )+
  #add in custom legend
  scale_colour_manual(name="Concentration",values=cols3)+
  theme2+
  theme(legend.position = "none")

#https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
library(grid)
library(gtable)
library(gridExtra)
g1 <- ggplotGrob(fig1)
g2 <- ggplotGrob(fig2)
g3 <- ggplotGrob(fig3)
g <- rbind(g1, g2, g3, size = "first")
g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths)

grid.newpage()
grid.draw(g)
