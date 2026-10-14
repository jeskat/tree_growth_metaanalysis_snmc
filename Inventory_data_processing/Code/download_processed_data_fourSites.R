# IMPORTANT: Due to new authentication requirements for direct access to EDI data you must provide an 
# authentication key. You can find instructions for obtaining a group-based access key at:
# https://youtu.be/fieZSmHk2H4?si=Wo9a5GsAOYp3dnWS 
# and the IAM web site at: 
# https://auth.edirepository.org
#
# Once you have the access key, just paste it between the quotes below and the code should run correctly and
# automatically download data from EDI for the code to use. 

myEDIAccessKey="q0oXZihyUx_Xre0jD-pmmy8yPac"   # put your access key between the quotes 

#
# Package ID: edi.2475.1 Cataloging System:https://pasta.edirepository.org.
# Data set title: Data for estimating the effects of fuel reduction treatments on tree growth at four study sites in Sierra Nevada mixed-conifer forests, California.
# Data set creator:  Jessica Katz -  
# Data set creator:  Eric Knapp -  
# Data set creator:  Kathryn Low -  
# Data set creator:  Patricia Manley -  
# Data set creator:  Malcolm North -  
# Data set creator:  Angela White -  
# Data set creator:  Lara Kueppers -  
# Metadata Provider:  Jessica Katz -  
# Contact:  Jessica Katz -    - jessica_katz@berkeley.edu
# Contact:  Lara Kueppers -    - lmkueppers@berkeley.edu
# Stylesheet v2.17 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu      
# Uncomment the following lines to have R clear previous work, or set a working directory
# rm(list=ls())      

# setwd("C:/users/my_name/my_dir")       

# check for required access key
if (myEDIAccessKey == ""){
  stop(paste0("A VALID EDI ACCESS KEY IS REQUIRED FOR THIS CODE TO DOWNLOAD THE NEEDED DATA \n",
              "Please set the myEDIAccessKey variable to a valid key and rerun this script"))
}

options(HTTPUserAgent="EDI_CodeGen")


inUrl1  <- paste0("https://pasta.lternet.edu/package/data/eml/edi/2475/1/0719fee56469d42f31439d893a901bae","?key=",myEDIAccessKey) 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl",extra=paste0(' -A "',getOption("HTTPUserAgent"),'"')))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "Site",     
                 "TreeID",     
                 "start_year",     
                 "end_year",     
                 "PlotID",     
                 "n_obs",     
                 "PFT"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$Site)!="factor") dt1$Site<- as.factor(dt1$Site)
if (class(dt1$TreeID)!="factor") dt1$TreeID<- as.factor(dt1$TreeID)
if (class(dt1$start_year)=="factor") dt1$start_year <-as.numeric(levels(dt1$start_year))[as.integer(dt1$start_year) ]               
if (class(dt1$start_year)=="character") dt1$start_year <-as.numeric(dt1$start_year)
if (class(dt1$end_year)=="factor") dt1$end_year <-as.numeric(levels(dt1$end_year))[as.integer(dt1$end_year) ]               
if (class(dt1$end_year)=="character") dt1$end_year <-as.numeric(dt1$end_year)
if (class(dt1$PlotID)!="factor") dt1$PlotID<- as.factor(dt1$PlotID)
if (class(dt1$n_obs)=="factor") dt1$n_obs <-as.numeric(levels(dt1$n_obs))[as.integer(dt1$n_obs) ]               
if (class(dt1$n_obs)=="character") dt1$n_obs <-as.numeric(dt1$n_obs)
if (class(dt1$PFT)!="factor") dt1$PFT<- as.factor(dt1$PFT)

# Convert Missing Values to NA for non-dates



# # Here is the structure of the input data frame:
# print("dt1) Structure")		    
# str(dt1)                            
# attach(dt1)                            
# # The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 
# 

# print(" ")
# print("Summary of Site")
# print(summary(Site))
# print(" ")
# print("Summary of TreeID")
# print(summary(TreeID))
# print(" ")
# print("Summary of start_year")
# print(summary(start_year))
# print(" ")
# print("Summary of end_year")
# print(summary(end_year))
# print(" ")
# print("Summary of PlotID")
# print(summary(PlotID))
# print(" ")
# print("Summary of n_obs")
# print(summary(n_obs))
# print(" ")
# print("Summary of PFT")
# print(summary(PFT)) 
# # Get more details on character variables
# 
# 
# print(" ")
# print("Summary of Site")
# print(summary(as.factor(dt1$Site))) 
# 
# print(" ")
# print("Summary of TreeID")
# print(summary(as.factor(dt1$TreeID))) 
# 
# print(" ")
# print("Summary of PlotID")
# print(summary(as.factor(dt1$PlotID))) 
# 
# print(" ")
# print("Summary of PFT")
# print(summary(as.factor(dt1$PFT)))
# detach(dt1)               



inUrl2  <- paste0("https://pasta.lternet.edu/package/data/eml/edi/2475/1/3d9180eb6a98fdd0dcba47e6e716f19b","?key=",myEDIAccessKey) 
infile2 <- tempfile()
try(download.file(inUrl2,infile2,method="curl",extra=paste0(' -A "',getOption("HTTPUserAgent"),'"')))
if (is.na(file.size(infile2))) download.file(inUrl2,infile2,method="auto")


dt2 <-read.csv(infile2,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "Site",     
                 "PlotID",     
                 "treatment",     
                 "trt_yr_1",     
                 "trt_yr_2"    ), check.names=TRUE)

unlink(infile2)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt2$Site)!="factor") dt2$Site<- as.factor(dt2$Site)
if (class(dt2$PlotID)!="factor") dt2$PlotID<- as.factor(dt2$PlotID)
if (class(dt2$treatment)!="factor") dt2$treatment<- as.factor(dt2$treatment)
if (class(dt2$trt_yr_1)=="factor") dt2$trt_yr_1 <-as.numeric(levels(dt2$trt_yr_1))[as.integer(dt2$trt_yr_1) ]               
if (class(dt2$trt_yr_1)=="character") dt2$trt_yr_1 <-as.numeric(dt2$trt_yr_1)
if (class(dt2$trt_yr_2)=="factor") dt2$trt_yr_2 <-as.numeric(levels(dt2$trt_yr_2))[as.integer(dt2$trt_yr_2) ]               
if (class(dt2$trt_yr_2)=="character") dt2$trt_yr_2 <-as.numeric(dt2$trt_yr_2)

# Convert Missing Values to NA for non-dates

dt2$trt_yr_1 <- ifelse((trimws(as.character(dt2$trt_yr_1))==trimws("NA")),NA,dt2$trt_yr_1)               
suppressWarnings(dt2$trt_yr_1 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt2$trt_yr_1))==as.character(as.numeric("NA"))),NA,dt2$trt_yr_1))
dt2$trt_yr_2 <- ifelse((trimws(as.character(dt2$trt_yr_2))==trimws("NA")),NA,dt2$trt_yr_2)               
suppressWarnings(dt2$trt_yr_2 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt2$trt_yr_2))==as.character(as.numeric("NA"))),NA,dt2$trt_yr_2))


# # Here is the structure of the input data frame:
# print("dt2) Structure")		    
# str(dt2)                            
# attach(dt2)                            
# # The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 
# 
# 
# print(" ")
# print("Summary of Site")
# print(summary(Site))
# print(" ")
# print("Summary of PlotID")
# print(summary(PlotID))
# print(" ")
# print("Summary of treatment")
# print(summary(treatment))
# print(" ")
# print("Summary of trt_yr_1")
# print(summary(trt_yr_1))
# print(" ")
# print("Summary of trt_yr_2")
# print(summary(trt_yr_2)) 
# # Get more details on character variables
# 
# 
# print(" ")
# print("Summary of Site")
# print(summary(as.factor(dt2$Site))) 
# 
# print(" ")
# print("Summary of PlotID")
# print(summary(as.factor(dt2$PlotID))) 
# 
# print(" ")
# print("Summary of treatment")
# print(summary(as.factor(dt2$treatment)))
# detach(dt2)               



inUrl3  <- paste0("https://pasta.lternet.edu/package/data/eml/edi/2475/1/92ad2ef013cf435a323836f57c6752eb","?key=",myEDIAccessKey) 
infile3 <- tempfile()
try(download.file(inUrl3,infile3,method="curl",extra=paste0(' -A "',getOption("HTTPUserAgent"),'"')))
if (is.na(file.size(infile3))) download.file(inUrl3,infile3,method="auto")


dt3 <-read.csv(infile3,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "Site",     
                 "TreeID",     
                 "1",     
                 "2",     
                 "3",     
                 "4",     
                 "5",     
                 "6"    ), check.names=TRUE)

unlink(infile3)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt3$Site)!="factor") dt3$Site<- as.factor(dt3$Site)
if (class(dt3$TreeID)!="factor") dt3$TreeID<- as.factor(dt3$TreeID)
if (class(dt3$X1)=="factor") dt3$X1 <-as.numeric(levels(dt3$X1))[as.integer(dt3$X1) ]               
if (class(dt3$X1)=="character") dt3$X1 <-as.numeric(dt3$X1)
if (class(dt3$X2)=="factor") dt3$X2 <-as.numeric(levels(dt3$X2))[as.integer(dt3$X2) ]               
if (class(dt3$X2)=="character") dt3$X2 <-as.numeric(dt3$X2)
if (class(dt3$X3)=="factor") dt3$X3 <-as.numeric(levels(dt3$X3))[as.integer(dt3$X3) ]               
if (class(dt3$X3)=="character") dt3$X3 <-as.numeric(dt3$X3)
if (class(dt3$X4)=="factor") dt3$X4 <-as.numeric(levels(dt3$X4))[as.integer(dt3$X4) ]               
if (class(dt3$X4)=="character") dt3$X4 <-as.numeric(dt3$X4)
if (class(dt3$X5)=="factor") dt3$X5 <-as.numeric(levels(dt3$X5))[as.integer(dt3$X5) ]               
if (class(dt3$X5)=="character") dt3$X5 <-as.numeric(dt3$X5)
if (class(dt3$X6)=="factor") dt3$X6 <-as.numeric(levels(dt3$X6))[as.integer(dt3$X6) ]               
if (class(dt3$X6)=="character") dt3$X6 <-as.numeric(dt3$X6)

# Convert Missing Values to NA for non-dates

dt3$X3 <- ifelse((trimws(as.character(dt3$X3))==trimws("NA")),NA,dt3$X3)               
suppressWarnings(dt3$X3 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt3$X3))==as.character(as.numeric("NA"))),NA,dt3$X3))
dt3$X4 <- ifelse((trimws(as.character(dt3$X4))==trimws("NA")),NA,dt3$X4)               
suppressWarnings(dt3$X4 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt3$X4))==as.character(as.numeric("NA"))),NA,dt3$X4))
dt3$X5 <- ifelse((trimws(as.character(dt3$X5))==trimws("NA")),NA,dt3$X5)               
suppressWarnings(dt3$X5 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt3$X5))==as.character(as.numeric("NA"))),NA,dt3$X5))
dt3$X6 <- ifelse((trimws(as.character(dt3$X6))==trimws("NA")),NA,dt3$X6)               
suppressWarnings(dt3$X6 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt3$X6))==as.character(as.numeric("NA"))),NA,dt3$X6))


# # Here is the structure of the input data frame:
# print("dt3) Structure")		    
# str(dt3)                            
# attach(dt3)                            
# # The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 
# 
# 
# print(" ")
# print("Summary of Site")
# print(summary(Site))
# print(" ")
# print("Summary of TreeID")
# print(summary(TreeID))
# print(" ")
# print("Summary of X1")
# print(summary(X1))
# print(" ")
# print("Summary of X2")
# print(summary(X2))
# print(" ")
# print("Summary of X3")
# print(summary(X3))
# print(" ")
# print("Summary of X4")
# print(summary(X4))
# print(" ")
# print("Summary of X5")
# print(summary(X5))
# print(" ")
# print("Summary of X6")
# print(summary(X6)) 
# # Get more details on character variables
# 
# 
# print(" ")
# print("Summary of Site")
# print(summary(as.factor(dt3$Site))) 
# 
# print(" ")
# print("Summary of TreeID")
# print(summary(as.factor(dt3$TreeID)))
# detach(dt3)               



inUrl4  <- paste0("https://pasta.lternet.edu/package/data/eml/edi/2475/1/403f4b208f5a878dfa3fd8203d1e17a3","?key=",myEDIAccessKey) 
infile4 <- tempfile()
try(download.file(inUrl4,infile4,method="curl",extra=paste0(' -A "',getOption("HTTPUserAgent"),'"')))
if (is.na(file.size(infile4))) download.file(inUrl4,infile4,method="auto")


dt4 <-read.csv(infile4,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "Site",     
                 "TreeID",     
                 "X1",     
                 "X2",     
                 "X3",     
                 "X4",     
                 "X5",     
                 "X6"    ), check.names=TRUE)

unlink(infile4)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt4$Site)!="factor") dt4$Site<- as.factor(dt4$Site)
if (class(dt4$TreeID)!="factor") dt4$TreeID<- as.factor(dt4$TreeID)
if (class(dt4$X1)=="factor") dt4$X1 <-as.numeric(levels(dt4$X1))[as.integer(dt4$X1) ]               
if (class(dt4$X1)=="character") dt4$X1 <-as.numeric(dt4$X1)
if (class(dt4$X2)=="factor") dt4$X2 <-as.numeric(levels(dt4$X2))[as.integer(dt4$X2) ]               
if (class(dt4$X2)=="character") dt4$X2 <-as.numeric(dt4$X2)
if (class(dt4$X3)=="factor") dt4$X3 <-as.numeric(levels(dt4$X3))[as.integer(dt4$X3) ]               
if (class(dt4$X3)=="character") dt4$X3 <-as.numeric(dt4$X3)
if (class(dt4$X4)=="factor") dt4$X4 <-as.numeric(levels(dt4$X4))[as.integer(dt4$X4) ]               
if (class(dt4$X4)=="character") dt4$X4 <-as.numeric(dt4$X4)
if (class(dt4$X5)=="factor") dt4$X5 <-as.numeric(levels(dt4$X5))[as.integer(dt4$X5) ]               
if (class(dt4$X5)=="character") dt4$X5 <-as.numeric(dt4$X5)
if (class(dt4$X6)=="factor") dt4$X6 <-as.numeric(levels(dt4$X6))[as.integer(dt4$X6) ]               
if (class(dt4$X6)=="character") dt4$X6 <-as.numeric(dt4$X6)

# Convert Missing Values to NA for non-dates

dt4$X3 <- ifelse((trimws(as.character(dt4$X3))==trimws("NA")),NA,dt4$X3)               
suppressWarnings(dt4$X3 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt4$X3))==as.character(as.numeric("NA"))),NA,dt4$X3))
dt4$X4 <- ifelse((trimws(as.character(dt4$X4))==trimws("NA")),NA,dt4$X4)               
suppressWarnings(dt4$X4 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt4$X4))==as.character(as.numeric("NA"))),NA,dt4$X4))
dt4$X5 <- ifelse((trimws(as.character(dt4$X5))==trimws("NA")),NA,dt4$X5)               
suppressWarnings(dt4$X5 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt4$X5))==as.character(as.numeric("NA"))),NA,dt4$X5))
dt4$X6 <- ifelse((trimws(as.character(dt4$X6))==trimws("NA")),NA,dt4$X6)               
suppressWarnings(dt4$X6 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt4$X6))==as.character(as.numeric("NA"))),NA,dt4$X6))


# # Here is the structure of the input data frame:
# print("dt4) Structure")		    
# str(dt4)                            
# attach(dt4)                            
# # The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 
# 
# 
# print(" ")
# print("Summary of Site")
# print(summary(Site))
# print(" ")
# print("Summary of TreeID")
# print(summary(TreeID))
# print(" ")
# print("Summary of X1")
# print(summary(X1))
# print(" ")
# print("Summary of X2")
# print(summary(X2))
# print(" ")
# print("Summary of X3")
# print(summary(X3))
# print(" ")
# print("Summary of X4")
# print(summary(X4))
# print(" ")
# print("Summary of X5")
# print(summary(X5))
# print(" ")
# print("Summary of X6")
# print(summary(X6)) 
# # Get more details on character variables
# 
# 
# print(" ")
# print("Summary of Site")
# print(summary(as.factor(dt4$Site))) 
# 
# print(" ")
# print("Summary of TreeID")
# print(summary(as.factor(dt4$TreeID)))
# detach(dt4)               



inUrl5  <- paste0("https://pasta.lternet.edu/package/data/eml/edi/2475/1/f47682cf4fd1ba03ee785c3df98477a7","?key=",myEDIAccessKey) 
infile5 <- tempfile()
try(download.file(inUrl5,infile5,method="curl",extra=paste0(' -A "',getOption("HTTPUserAgent"),'"')))
if (is.na(file.size(infile5))) download.file(inUrl5,infile5,method="auto")


dt5 <-read.csv(infile5,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "Site",     
                 "TreeID",     
                 "X1",     
                 "X2",     
                 "X3",     
                 "X4",     
                 "X5",     
                 "X6",     
                 "X7",     
                 "X8",     
                 "X9",     
                 "X10",     
                 "X11",     
                 "X12",     
                 "X13",     
                 "X14",     
                 "X15",     
                 "X16",     
                 "X17",     
                 "X18"    ), check.names=TRUE)

unlink(infile5)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt5$Site)!="factor") dt5$Site<- as.factor(dt5$Site)
if (class(dt5$TreeID)!="factor") dt5$TreeID<- as.factor(dt5$TreeID)
if (class(dt5$X1)=="factor") dt5$X1 <-as.numeric(levels(dt5$X1))[as.integer(dt5$X1) ]               
if (class(dt5$X1)=="character") dt5$X1 <-as.numeric(dt5$X1)
if (class(dt5$X2)=="factor") dt5$X2 <-as.numeric(levels(dt5$X2))[as.integer(dt5$X2) ]               
if (class(dt5$X2)=="character") dt5$X2 <-as.numeric(dt5$X2)
if (class(dt5$X3)=="factor") dt5$X3 <-as.numeric(levels(dt5$X3))[as.integer(dt5$X3) ]               
if (class(dt5$X3)=="character") dt5$X3 <-as.numeric(dt5$X3)
if (class(dt5$X4)=="factor") dt5$X4 <-as.numeric(levels(dt5$X4))[as.integer(dt5$X4) ]               
if (class(dt5$X4)=="character") dt5$X4 <-as.numeric(dt5$X4)
if (class(dt5$X5)=="factor") dt5$X5 <-as.numeric(levels(dt5$X5))[as.integer(dt5$X5) ]               
if (class(dt5$X5)=="character") dt5$X5 <-as.numeric(dt5$X5)
if (class(dt5$X6)=="factor") dt5$X6 <-as.numeric(levels(dt5$X6))[as.integer(dt5$X6) ]               
if (class(dt5$X6)=="character") dt5$X6 <-as.numeric(dt5$X6)
if (class(dt5$X7)=="factor") dt5$X7 <-as.numeric(levels(dt5$X7))[as.integer(dt5$X7) ]               
if (class(dt5$X7)=="character") dt5$X7 <-as.numeric(dt5$X7)
if (class(dt5$X8)=="factor") dt5$X8 <-as.numeric(levels(dt5$X8))[as.integer(dt5$X8) ]               
if (class(dt5$X8)=="character") dt5$X8 <-as.numeric(dt5$X8)
if (class(dt5$X9)=="factor") dt5$X9 <-as.numeric(levels(dt5$X9))[as.integer(dt5$X9) ]               
if (class(dt5$X9)=="character") dt5$X9 <-as.numeric(dt5$X9)
if (class(dt5$X10)=="factor") dt5$X10 <-as.numeric(levels(dt5$X10))[as.integer(dt5$X10) ]               
if (class(dt5$X10)=="character") dt5$X10 <-as.numeric(dt5$X10)
if (class(dt5$X11)=="factor") dt5$X11 <-as.numeric(levels(dt5$X11))[as.integer(dt5$X11) ]               
if (class(dt5$X11)=="character") dt5$X11 <-as.numeric(dt5$X11)
if (class(dt5$X12)=="factor") dt5$X12 <-as.numeric(levels(dt5$X12))[as.integer(dt5$X12) ]               
if (class(dt5$X12)=="character") dt5$X12 <-as.numeric(dt5$X12)
if (class(dt5$X13)=="factor") dt5$X13 <-as.numeric(levels(dt5$X13))[as.integer(dt5$X13) ]               
if (class(dt5$X13)=="character") dt5$X13 <-as.numeric(dt5$X13)
if (class(dt5$X14)=="factor") dt5$X14 <-as.numeric(levels(dt5$X14))[as.integer(dt5$X14) ]               
if (class(dt5$X14)=="character") dt5$X14 <-as.numeric(dt5$X14)
if (class(dt5$X15)=="factor") dt5$X15 <-as.numeric(levels(dt5$X15))[as.integer(dt5$X15) ]               
if (class(dt5$X15)=="character") dt5$X15 <-as.numeric(dt5$X15)
if (class(dt5$X16)=="factor") dt5$X16 <-as.numeric(levels(dt5$X16))[as.integer(dt5$X16) ]               
if (class(dt5$X16)=="character") dt5$X16 <-as.numeric(dt5$X16)
if (class(dt5$X17)=="factor") dt5$X17 <-as.numeric(levels(dt5$X17))[as.integer(dt5$X17) ]               
if (class(dt5$X17)=="character") dt5$X17 <-as.numeric(dt5$X17)
if (class(dt5$X18)=="factor") dt5$X18 <-as.numeric(levels(dt5$X18))[as.integer(dt5$X18) ]               
if (class(dt5$X18)=="character") dt5$X18 <-as.numeric(dt5$X18)

# Convert Missing Values to NA for non-dates

dt5$X1 <- ifelse((trimws(as.character(dt5$X1))==trimws("NA")),NA,dt5$X1)               
suppressWarnings(dt5$X1 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt5$X1))==as.character(as.numeric("NA"))),NA,dt5$X1))
dt5$X2 <- ifelse((trimws(as.character(dt5$X2))==trimws("NA")),NA,dt5$X2)               
suppressWarnings(dt5$X2 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt5$X2))==as.character(as.numeric("NA"))),NA,dt5$X2))
dt5$X3 <- ifelse((trimws(as.character(dt5$X3))==trimws("NA")),NA,dt5$X3)               
suppressWarnings(dt5$X3 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt5$X3))==as.character(as.numeric("NA"))),NA,dt5$X3))
dt5$X4 <- ifelse((trimws(as.character(dt5$X4))==trimws("NA")),NA,dt5$X4)               
suppressWarnings(dt5$X4 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt5$X4))==as.character(as.numeric("NA"))),NA,dt5$X4))
dt5$X5 <- ifelse((trimws(as.character(dt5$X5))==trimws("NA")),NA,dt5$X5)               
suppressWarnings(dt5$X5 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt5$X5))==as.character(as.numeric("NA"))),NA,dt5$X5))
dt5$X6 <- ifelse((trimws(as.character(dt5$X6))==trimws("NA")),NA,dt5$X6)               
suppressWarnings(dt5$X6 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt5$X6))==as.character(as.numeric("NA"))),NA,dt5$X6))
dt5$X7 <- ifelse((trimws(as.character(dt5$X7))==trimws("NA")),NA,dt5$X7)               
suppressWarnings(dt5$X7 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt5$X7))==as.character(as.numeric("NA"))),NA,dt5$X7))
dt5$X8 <- ifelse((trimws(as.character(dt5$X8))==trimws("NA")),NA,dt5$X8)               
suppressWarnings(dt5$X8 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt5$X8))==as.character(as.numeric("NA"))),NA,dt5$X8))
dt5$X9 <- ifelse((trimws(as.character(dt5$X9))==trimws("NA")),NA,dt5$X9)               
suppressWarnings(dt5$X9 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt5$X9))==as.character(as.numeric("NA"))),NA,dt5$X9))
dt5$X10 <- ifelse((trimws(as.character(dt5$X10))==trimws("NA")),NA,dt5$X10)               
suppressWarnings(dt5$X10 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt5$X10))==as.character(as.numeric("NA"))),NA,dt5$X10))
dt5$X11 <- ifelse((trimws(as.character(dt5$X11))==trimws("NA")),NA,dt5$X11)               
suppressWarnings(dt5$X11 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt5$X11))==as.character(as.numeric("NA"))),NA,dt5$X11))
dt5$X12 <- ifelse((trimws(as.character(dt5$X12))==trimws("NA")),NA,dt5$X12)               
suppressWarnings(dt5$X12 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt5$X12))==as.character(as.numeric("NA"))),NA,dt5$X12))
dt5$X13 <- ifelse((trimws(as.character(dt5$X13))==trimws("NA")),NA,dt5$X13)               
suppressWarnings(dt5$X13 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt5$X13))==as.character(as.numeric("NA"))),NA,dt5$X13))
dt5$X14 <- ifelse((trimws(as.character(dt5$X14))==trimws("NA")),NA,dt5$X14)               
suppressWarnings(dt5$X14 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt5$X14))==as.character(as.numeric("NA"))),NA,dt5$X14))
dt5$X15 <- ifelse((trimws(as.character(dt5$X15))==trimws("NA")),NA,dt5$X15)               
suppressWarnings(dt5$X15 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt5$X15))==as.character(as.numeric("NA"))),NA,dt5$X15))
dt5$X16 <- ifelse((trimws(as.character(dt5$X16))==trimws("NA")),NA,dt5$X16)               
suppressWarnings(dt5$X16 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt5$X16))==as.character(as.numeric("NA"))),NA,dt5$X16))
dt5$X17 <- ifelse((trimws(as.character(dt5$X17))==trimws("NA")),NA,dt5$X17)               
suppressWarnings(dt5$X17 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt5$X17))==as.character(as.numeric("NA"))),NA,dt5$X17))
dt5$X18 <- ifelse((trimws(as.character(dt5$X18))==trimws("NA")),NA,dt5$X18)               
suppressWarnings(dt5$X18 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt5$X18))==as.character(as.numeric("NA"))),NA,dt5$X18))


# # Here is the structure of the input data frame:
# print("dt5) Structure")		    
# str(dt5)                            
# attach(dt5)                            
# # The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 
# 
# 
# print(" ")
# print("Summary of Site")
# print(summary(Site))
# print(" ")
# print("Summary of TreeID")
# print(summary(TreeID))
# print(" ")
# print("Summary of X1")
# print(summary(X1))
# print(" ")
# print("Summary of X2")
# print(summary(X2))
# print(" ")
# print("Summary of X3")
# print(summary(X3))
# print(" ")
# print("Summary of X4")
# print(summary(X4))
# print(" ")
# print("Summary of X5")
# print(summary(X5))
# print(" ")
# print("Summary of X6")
# print(summary(X6))
# print(" ")
# print("Summary of X7")
# print(summary(X7))
# print(" ")
# print("Summary of X8")
# print(summary(X8))
# print(" ")
# print("Summary of X9")
# print(summary(X9))
# print(" ")
# print("Summary of X10")
# print(summary(X10))
# print(" ")
# print("Summary of X11")
# print(summary(X11))
# print(" ")
# print("Summary of X12")
# print(summary(X12))
# print(" ")
# print("Summary of X13")
# print(summary(X13))
# print(" ")
# print("Summary of X14")
# print(summary(X14))
# print(" ")
# print("Summary of X15")
# print(summary(X15))
# print(" ")
# print("Summary of X16")
# print(summary(X16))
# print(" ")
# print("Summary of X17")
# print(summary(X17))
# print(" ")
# print("Summary of X18")
# print(summary(X18)) 
# # Get more details on character variables
# 
# 
# print(" ")
# print("Summary of Site")
# print(summary(as.factor(dt5$Site))) 
# 
# print(" ")
# print("Summary of TreeID")
# print(summary(as.factor(dt5$TreeID)))
# detach(dt5)               



inUrl6  <- paste0("https://pasta.lternet.edu/package/data/eml/edi/2475/1/660adacbac4c1925723ed5af94ba9e1d","?key=",myEDIAccessKey) 
infile6 <- tempfile()
try(download.file(inUrl6,infile6,method="curl",extra=paste0(' -A "',getOption("HTTPUserAgent"),'"')))
if (is.na(file.size(infile6))) download.file(inUrl6,infile6,method="auto")


dt6 <-read.csv(infile6,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "Site",     
                 "TreeID",     
                 "X1",     
                 "X2",     
                 "X3",     
                 "X4",     
                 "X5",     
                 "X6",     
                 "X7",     
                 "X8",     
                 "X9",     
                 "X10",     
                 "X11",     
                 "X12",     
                 "X13",     
                 "X14",     
                 "X15",     
                 "X16",     
                 "X17",     
                 "X18"    ), check.names=TRUE)

unlink(infile6)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt6$Site)!="factor") dt6$Site<- as.factor(dt6$Site)
if (class(dt6$TreeID)!="factor") dt6$TreeID<- as.factor(dt6$TreeID)
if (class(dt6$X1)=="factor") dt6$X1 <-as.numeric(levels(dt6$X1))[as.integer(dt6$X1) ]               
if (class(dt6$X1)=="character") dt6$X1 <-as.numeric(dt6$X1)
if (class(dt6$X2)=="factor") dt6$X2 <-as.numeric(levels(dt6$X2))[as.integer(dt6$X2) ]               
if (class(dt6$X2)=="character") dt6$X2 <-as.numeric(dt6$X2)
if (class(dt6$X3)=="factor") dt6$X3 <-as.numeric(levels(dt6$X3))[as.integer(dt6$X3) ]               
if (class(dt6$X3)=="character") dt6$X3 <-as.numeric(dt6$X3)
if (class(dt6$X4)=="factor") dt6$X4 <-as.numeric(levels(dt6$X4))[as.integer(dt6$X4) ]               
if (class(dt6$X4)=="character") dt6$X4 <-as.numeric(dt6$X4)
if (class(dt6$X5)=="factor") dt6$X5 <-as.numeric(levels(dt6$X5))[as.integer(dt6$X5) ]               
if (class(dt6$X5)=="character") dt6$X5 <-as.numeric(dt6$X5)
if (class(dt6$X6)=="factor") dt6$X6 <-as.numeric(levels(dt6$X6))[as.integer(dt6$X6) ]               
if (class(dt6$X6)=="character") dt6$X6 <-as.numeric(dt6$X6)
if (class(dt6$X7)=="factor") dt6$X7 <-as.numeric(levels(dt6$X7))[as.integer(dt6$X7) ]               
if (class(dt6$X7)=="character") dt6$X7 <-as.numeric(dt6$X7)
if (class(dt6$X8)=="factor") dt6$X8 <-as.numeric(levels(dt6$X8))[as.integer(dt6$X8) ]               
if (class(dt6$X8)=="character") dt6$X8 <-as.numeric(dt6$X8)
if (class(dt6$X9)=="factor") dt6$X9 <-as.numeric(levels(dt6$X9))[as.integer(dt6$X9) ]               
if (class(dt6$X9)=="character") dt6$X9 <-as.numeric(dt6$X9)
if (class(dt6$X10)=="factor") dt6$X10 <-as.numeric(levels(dt6$X10))[as.integer(dt6$X10) ]               
if (class(dt6$X10)=="character") dt6$X10 <-as.numeric(dt6$X10)
if (class(dt6$X11)=="factor") dt6$X11 <-as.numeric(levels(dt6$X11))[as.integer(dt6$X11) ]               
if (class(dt6$X11)=="character") dt6$X11 <-as.numeric(dt6$X11)
if (class(dt6$X12)=="factor") dt6$X12 <-as.numeric(levels(dt6$X12))[as.integer(dt6$X12) ]               
if (class(dt6$X12)=="character") dt6$X12 <-as.numeric(dt6$X12)
if (class(dt6$X13)=="factor") dt6$X13 <-as.numeric(levels(dt6$X13))[as.integer(dt6$X13) ]               
if (class(dt6$X13)=="character") dt6$X13 <-as.numeric(dt6$X13)
if (class(dt6$X14)=="factor") dt6$X14 <-as.numeric(levels(dt6$X14))[as.integer(dt6$X14) ]               
if (class(dt6$X14)=="character") dt6$X14 <-as.numeric(dt6$X14)
if (class(dt6$X15)=="factor") dt6$X15 <-as.numeric(levels(dt6$X15))[as.integer(dt6$X15) ]               
if (class(dt6$X15)=="character") dt6$X15 <-as.numeric(dt6$X15)
if (class(dt6$X16)=="factor") dt6$X16 <-as.numeric(levels(dt6$X16))[as.integer(dt6$X16) ]               
if (class(dt6$X16)=="character") dt6$X16 <-as.numeric(dt6$X16)
if (class(dt6$X17)=="factor") dt6$X17 <-as.numeric(levels(dt6$X17))[as.integer(dt6$X17) ]               
if (class(dt6$X17)=="character") dt6$X17 <-as.numeric(dt6$X17)
if (class(dt6$X18)=="factor") dt6$X18 <-as.numeric(levels(dt6$X18))[as.integer(dt6$X18) ]               
if (class(dt6$X18)=="character") dt6$X18 <-as.numeric(dt6$X18)

# Convert Missing Values to NA for non-dates

dt6$X1 <- ifelse((trimws(as.character(dt6$X1))==trimws("NA")),NA,dt6$X1)               
suppressWarnings(dt6$X1 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt6$X1))==as.character(as.numeric("NA"))),NA,dt6$X1))
dt6$X2 <- ifelse((trimws(as.character(dt6$X2))==trimws("NA")),NA,dt6$X2)               
suppressWarnings(dt6$X2 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt6$X2))==as.character(as.numeric("NA"))),NA,dt6$X2))
dt6$X3 <- ifelse((trimws(as.character(dt6$X3))==trimws("NA")),NA,dt6$X3)               
suppressWarnings(dt6$X3 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt6$X3))==as.character(as.numeric("NA"))),NA,dt6$X3))
dt6$X4 <- ifelse((trimws(as.character(dt6$X4))==trimws("NA")),NA,dt6$X4)               
suppressWarnings(dt6$X4 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt6$X4))==as.character(as.numeric("NA"))),NA,dt6$X4))
dt6$X5 <- ifelse((trimws(as.character(dt6$X5))==trimws("NA")),NA,dt6$X5)               
suppressWarnings(dt6$X5 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt6$X5))==as.character(as.numeric("NA"))),NA,dt6$X5))
dt6$X6 <- ifelse((trimws(as.character(dt6$X6))==trimws("NA")),NA,dt6$X6)               
suppressWarnings(dt6$X6 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt6$X6))==as.character(as.numeric("NA"))),NA,dt6$X6))
dt6$X7 <- ifelse((trimws(as.character(dt6$X7))==trimws("NA")),NA,dt6$X7)               
suppressWarnings(dt6$X7 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt6$X7))==as.character(as.numeric("NA"))),NA,dt6$X7))
dt6$X8 <- ifelse((trimws(as.character(dt6$X8))==trimws("NA")),NA,dt6$X8)               
suppressWarnings(dt6$X8 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt6$X8))==as.character(as.numeric("NA"))),NA,dt6$X8))
dt6$X9 <- ifelse((trimws(as.character(dt6$X9))==trimws("NA")),NA,dt6$X9)               
suppressWarnings(dt6$X9 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt6$X9))==as.character(as.numeric("NA"))),NA,dt6$X9))
dt6$X10 <- ifelse((trimws(as.character(dt6$X10))==trimws("NA")),NA,dt6$X10)               
suppressWarnings(dt6$X10 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt6$X10))==as.character(as.numeric("NA"))),NA,dt6$X10))
dt6$X11 <- ifelse((trimws(as.character(dt6$X11))==trimws("NA")),NA,dt6$X11)               
suppressWarnings(dt6$X11 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt6$X11))==as.character(as.numeric("NA"))),NA,dt6$X11))
dt6$X12 <- ifelse((trimws(as.character(dt6$X12))==trimws("NA")),NA,dt6$X12)               
suppressWarnings(dt6$X12 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt6$X12))==as.character(as.numeric("NA"))),NA,dt6$X12))
dt6$X13 <- ifelse((trimws(as.character(dt6$X13))==trimws("NA")),NA,dt6$X13)               
suppressWarnings(dt6$X13 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt6$X13))==as.character(as.numeric("NA"))),NA,dt6$X13))
dt6$X14 <- ifelse((trimws(as.character(dt6$X14))==trimws("NA")),NA,dt6$X14)               
suppressWarnings(dt6$X14 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt6$X14))==as.character(as.numeric("NA"))),NA,dt6$X14))
dt6$X15 <- ifelse((trimws(as.character(dt6$X15))==trimws("NA")),NA,dt6$X15)               
suppressWarnings(dt6$X15 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt6$X15))==as.character(as.numeric("NA"))),NA,dt6$X15))
dt6$X16 <- ifelse((trimws(as.character(dt6$X16))==trimws("NA")),NA,dt6$X16)               
suppressWarnings(dt6$X16 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt6$X16))==as.character(as.numeric("NA"))),NA,dt6$X16))
dt6$X17 <- ifelse((trimws(as.character(dt6$X17))==trimws("NA")),NA,dt6$X17)               
suppressWarnings(dt6$X17 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt6$X17))==as.character(as.numeric("NA"))),NA,dt6$X17))
dt6$X18 <- ifelse((trimws(as.character(dt6$X18))==trimws("NA")),NA,dt6$X18)               
suppressWarnings(dt6$X18 <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt6$X18))==as.character(as.numeric("NA"))),NA,dt6$X18))


# # Here is the structure of the input data frame:
# print("dt6) Structure")		    
# str(dt6)                            
# attach(dt6)                            
# # The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 
# 
# 
# print(" ")
# print("Summary of Site")
# print(summary(Site))
# print(" ")
# print("Summary of TreeID")
# print(summary(TreeID))
# print(" ")
# print("Summary of X1")
# print(summary(X1))
# print(" ")
# print("Summary of X2")
# print(summary(X2))
# print(" ")
# print("Summary of X3")
# print(summary(X3))
# print(" ")
# print("Summary of X4")
# print(summary(X4))
# print(" ")
# print("Summary of X5")
# print(summary(X5))
# print(" ")
# print("Summary of X6")
# print(summary(X6))
# print(" ")
# print("Summary of X7")
# print(summary(X7))
# print(" ")
# print("Summary of X8")
# print(summary(X8))
# print(" ")
# print("Summary of X9")
# print(summary(X9))
# print(" ")
# print("Summary of X10")
# print(summary(X10))
# print(" ")
# print("Summary of X11")
# print(summary(X11))
# print(" ")
# print("Summary of X12")
# print(summary(X12))
# print(" ")
# print("Summary of X13")
# print(summary(X13))
# print(" ")
# print("Summary of X14")
# print(summary(X14))
# print(" ")
# print("Summary of X15")
# print(summary(X15))
# print(" ")
# print("Summary of X16")
# print(summary(X16))
# print(" ")
# print("Summary of X17")
# print(summary(X17))
# print(" ")
# print("Summary of X18")
# print(summary(X18)) 
# # Get more details on character variables
# 
# 
# print(" ")
# print("Summary of Site")
# print(summary(as.factor(dt6$Site))) 
# 
# print(" ")
# print("Summary of TreeID")
# print(summary(as.factor(dt6$TreeID)))
# detach(dt6)               
