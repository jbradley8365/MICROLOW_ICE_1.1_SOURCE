
## MicroLow_ICE_1.1

## 



# James Bradley
# jbradley.earth@gmail.com
# Queen Mary University of London, London, UK
# GFZ German Research Centre for Geosciences, Potsdam, Germany


# EXECUTE THIS CODE TO RUN MODEL
# SWITCH DEPENDS ON TEMPERATURE. TEMPERATURE IS VARIABLE - INPUT VIA "icetemp.csv".


# 1.0 SET UP WORKSPACE

# Clear workspace

rm(list = ls())

# Set working directory

setwd("~/Desktop/RFolder/MICROLOW_ICE_1.1_SOURCE")

# Download, install and run Packages


library("rootSolve", lib="/Users/jamesbradley/Desktop/RFolder/New_Packages")
library("shape", lib="/Users/jamesbradley/Desktop/RFolder/New_Packages")
library("deSolve", lib="/Users/jamesbradley/Desktop/RFolder/New_Packages")
library("ReacTran", lib="/Users/jamesbradley/Desktop/RFolder/New_Packages")
library("bvpSolve", lib="/Users/jamesbradley/Desktop/RFolder/New_Packages")
library("scatterplot3d", lib="/Users/jamesbradley/Desktop/RFolder/New_Packages")


######
######

# 2.0 SPECIFY PATHS

path    <-"~/Desktop/RFolder/MICROLOW_ICE_1.1_SOURCE/"
pathte  <-"~/Desktop/RFolder/MICROLOW_ICE_1.1_SOURCE/"


######
######


# 3.0 LOAD VALIDATION DATA

######
######

# 4.0 DEFINE TIMEPERIOD

# nhours is in hours

nhours<-240

times<-seq(0,nhours-1,by=1)      # Times at which model provides output

######
######


# 5.0 SET PARAMETER VALUES

V_max_B1 <-  0.048       # Maximum growth rate of B1
K_V_B1 <- 8000            # Half-saturation constant (DOC) for growth
YG <- 0.2                  # True growth yield
mq_B1 <-  0.0005           # Maintenance demand B1
mq_B2 <-  0.00005          # Maintenance demand B2
alpha_B1 <- 0.004       # Mortality rate B1
alpha_B2 <- 0.0004       # Mortality rate B2

st_S <- 2                # Steepness of state-change function
K_IceTemp_S <- 0.1         # Threshold ice temp state-change
R_S_D  <- 0.04           # Rate constant, deactivation  
R_S_A  <- R_S_D            # Rate constant, activation  
#st_M <- 0.1                # Steepness of maintenance function



parms<-c(
  V_max_B1,
  K_V_B1,
  YG,
  mq_B1,
  mq_B2,
  alpha_B1,
  alpha_B2,
  st_S,
  K_IceTemp_S,
  R_S_D, 
  R_S_A
)



######
######

# 6.0 LOAD DRIVERS

#ice_temp=-2.0;

ice_temp_driver <-paste(path,"icetemp.csv",sep="")
driverfile1<-read.csv(ice_temp_driver,header=TRUE)
driver_icetemp=data.frame(1,1)
driver_icetemp[1:dim(driverfile1)[1],1] = 1:(1*dim(driverfile1)[1])
driver_icetemp[1:dim(driverfile1)[1],2] = driverfile1[,2]



# 7.0 INITIAL VALUES

start<-c(B1=0.00,
         B2=0.55,
         #         B3=0.0,
         #         B4=0.0,
         Corg=1200,
         c_Cons_Corg_Growth_B1 = 0,
         c_Death_total = 0,
         c_M_total_Corg = 0,
         #         c_M_total_bio = 0,
         c_M_total = 0,
         #         c_M_B1_bio = 0,
         c_M_B1_Corg = 0,
         c_Cons_Corg_total = 0,
         c_Growth_B1 = 0,
         c_Death_B1 = 0,
         c_B1_D = 0,
         c_B2_A = 0,
         c_Death_B2=0
         #         c_Death_B3=0,
         #         c_Death_B4=0,
         #         c_B2_D=0,
         #         c_B3_D=0,
         #         c_B3_A=0,
         #         c_B4_A=0,
         #         c_M_B2_bio=0,
         #         c_M_B3_bio=0,
         #         c_M_B4_bio=0
)


#.................................................................................


# 8.0 CONSTRUCT ARRAYS FOR OUTPUT

#out_list = list()


# 9.0 MODEL DEFINITION AND EXECUTION WITHIN COUNTER

#for(counter in 1:1) {
  
  # 9.1 BEGINNING OF MODEL DEFINITION
  
  model<-function(t,xx,parms){
    
    B1<-xx[1]
    B2<-xx[2]
    #    B3<-xx[3]
    #    B4<-xx[4]
    Corg<-xx[3]
    c_Cons_Corg_Growth_B1<-xx[4]
    c_Death_total<-xx[5]
    c_M_total_Corg<-xx[6]
    #    c_M_total_bio<-xx[7]
    c_M_total<-xx[7]
    #    c_M_B1_bio<-xx[9]
    c_M_B1_Corg<-xx[8]
    c_Cons_Corg_total <-xx[9]
    c_Growth_B1<-xx[10]
    c_Death_B1<-xx[11]
    B1_D<-xx[12]
    B2_A<-xx[13]
    Death_B2<-xx[14]

    
    with(as.list(parms),{
      
      # Switch function for state change
      Theta_S <- 1/(exp((-driver_icetemp[t+1,2]+K_IceTemp_S)/(st_S*K_IceTemp_S))+1)        
      
 
      
      # Growth
      Growth_B1 <- V_max_B1*B1*(Corg/(Corg+K_V_B1))
      
      # Death
      Death_B1 <- alpha_B1*B1
      Death_B2 <- alpha_B2*B2

      
      Death_total = Death_B1 + Death_B2 # + Death_B3 + Death_B4  
      
      # Activation and deactivation 
      
      B2_A <- Theta_S*R_S_A*B2
      B1_D <- (1-Theta_S)*R_S_D*B1

      
      # Consumption of substrate due to active biomass growth
      
      Cons_Corg_Growth_B1 <- Growth_B1*(1/YG) 
      

      
      # Exogenous maintenance
      
      M_B1_Corg = mq_B1*B1#*Theta_M
      M_B2_Corg = mq_B2*B2#*Theta_M  
 
      
      M_total_Corg = M_B1_Corg + M_B2_Corg #+ M_B3_Corg + M_B4_Corg 
      
      # Total maintenance (from biomass and from organic carbon)
      
      M_total = M_total_Corg
      
      # Total substrate consumption
      
      Cons_Corg_total = Cons_Corg_Growth_B1 + M_total_Corg
      
      
      # BALANCE EQUATIONS
      
      dB1 <-  Growth_B1 - B1_D - Death_B1 + B2_A #
      
      dB2 <-  B1_D - Death_B2 - B2_A #
      

      
      dCorg <- Death_total - Cons_Corg_total
      
      
      #Derived variables
      
      dc_Cons_Corg_Growth_B1 <- Cons_Corg_Growth_B1
      
      dc_Death_total <- Death_total
      
      dc_M_total_Corg <- M_total_Corg
      
  
      
      dc_M_total <- M_total
      

      dc_M_B1_Corg <- M_B1_Corg
      
      dc_Cons_Corg_total <- Cons_Corg_total
      
      dc_Growth_B1 <- Growth_B1
      dc_Death_B1 <- Death_B1
      dc_B1_D <- B1_D
      dc_B2_A <- B2_A
      
      dc_Death_B2<- Death_B2


      
      
      # List the state variables and derived variables for which you want output
      
      list(c(dB1, dB2, dCorg, dc_Cons_Corg_Growth_B1, dc_Death_total, dc_M_total_Corg, dc_M_total, dc_M_B1_Corg,dc_Cons_Corg_total,dc_Growth_B1, dc_Death_B1, dc_B1_D, dc_B2_A,dc_Death_B2))
      
    })
  }
  
  
  # 9.2 END OF MODEL DEFINITION 
  
  
  out<-as.data.frame(lsoda(start,times,model,parms))

  # COMPUTE TOTALS AND ASSIGN TO NEW VARIABLES
  
              
  
  out$Btotal<- out$B1 + out$B2 #+ out$B3  + out$B4 
  
  

# END OF MODEL RUN

#####
#####

# 10.0 PLOTTING AND RESULTS

par(mfrow = c(2, 2))


plot(out$time,out$B1,type='l', xlab='time hours', main='B1 biomass', ylab='ugC/cm3')
plot(out$time,out$B2,type='l', xlab='time hours', main='B2 biomass', ylab='ugC/cm3')


plot(out$time,out$Btotal,type='l', xlab='time hours',main='Total biomass', ylab='ugC/cm3')


plot(out$time,out$Corg,type='l', xlab='time hours', main='Organic C', ylab='ugC/cm3')


