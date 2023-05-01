library(MASS)
library(R.matlab)
library(rmatio)
# library(GeodRegr)
library(optmatch)
library(senstrat)
library(abind)
library(ggplot2)


####################################
## data load and preprocessing
####################################
setwd('C:/Users/qsoon/Desktop/AATE-AMTE')

# covariates
# treated (AD=1)
ADF_x_data <- readMat('data/ADinfoF.mat')$ADinfoF # x data for women with Alzheimer's / 88 x 10
ADM_x_data <- readMat('data/ADinfoM.mat')$ADinfoM # x data for men with Alzheimer's / 98 x 10
AD_x_data <- rbind(ADF_x_data, ADM_x_data) # dim: 186 x 10

# controlled (AD=0)
NCF_x_data <- readMat('data/NCinfoF.mat')$NCinfoF # x data for women w/o Alzheimer's / 107 x 10
NCM_x_data <- readMat('data/NCinfoM.mat')$NCinfoM # x data for men w/o Alzheimer's / 116 x 10
NC_x_data <- rbind(NCF_x_data, NCM_x_data) # dim: 223 x 10

colnames(AD_x_data) <- c("rid", "gender", "handedness", "marital1", "marital2", "marital3", 
                          "educationlength", "retirement", "age", "AD")
colnames(NC_x_data) <- c("rid", "gender", "handedness", "marital1", "marital2", "marital3", 
                         "educationlength", "retirement", "age", "AD")

x_data <- data.frame(rbind(AD_x_data, NC_x_data))
x_data$AD <- ifelse(x_data$AD==3, 1, 0) # 1 = with Alzheimer, 0 = w/o Alzheimer
x_data$gender <- ifelse(x_data$gender==2, 1, 0) # 1 = Female, 0 = Male
sbj_id <- x_data[,1] # subject ID

# outcome
# treated (AD=1)
ADF_y_data <- readMat('data/ADLdataF.mat')$ADLdataF # y data for women with Alzheimer's / dim: 50 x 2 x 88
ADM_y_data <- readMat('data/ADLdataM.mat')$ADLdataM # y data for men with Alzheimer's / dim: 50 x 2 x 98

# controlled (AD=0)
NCF_y_data <- readMat('data/NCLdataF.mat')$NCLdataF # y data for women w/o Alzheimer's / dim: 50 x 2 x 107
NCM_y_data <- readMat('data/NCLdataM.mat')$NCLdataM # y data for men w/o Alzheimer's / dim: 50 x 2 x 116


y_data <- abind(ADF_y_data, ADM_y_data, NCF_y_data, NCM_y_data, along=3)
y_data <- aperm(y_data, c(2,1,3))


sum(is.na(x_data))
sum(is.na(y_data))


#################
## matching
#################
propscore.model <- glm(AD ~ gender + handedness + marital1 + marital2 + marital3 + 
                         educationlength + retirement + age,
                       family=binomial, x=TRUE, data=x_data)

x_data.trt <- x_data[x_data$AD==1,]
x_data.ctrl <- x_data[x_data$AD==0,]

prop.trt <- predict(propscore.model, newdata=x_data.trt[,2:9], type="response")
prop.ctrl <- predict(propscore.model, newdata=x_data.ctrl[,2:9], type="response")

ggplot(data.frame(Z = c(rep("treatment", length(prop.trt)), rep("control", length(prop.ctrl))), 
                  propensity_score = c(prop.trt, prop.ctrl)), aes(x = Z, y = propensity_score)) + geom_boxplot()



remove.trt.index <- as.numeric(names(which(prop.trt > max(prop.ctrl)))) # Indices of treated subjects to be removed 
remove.ctrl.index <- as.numeric(names(which(prop.ctrl < min(prop.trt)))) # Indices of controlled subjects to be removed 

remove.trt.id <- sbj_id[remove.trt.index] # ID of treated subjects to be removed 
remove.ctrl.id <- sbj_id[remove.ctrl.index] # ID of controlled subjects to be removed 


rankmahal <- function(z,X){
  X <- as.matrix(X)
  n <- dim(X)[1]
  rownames(X) <- 1:n
  k <- dim(X)[2]
  m <- sum(z)
  
  for(j in 1:k){
    X[,j] <- rank(X[,j])
  }
  cv <- cov(X)
  vuntied <- var(1:n)
  rat <- sqrt(vuntied/diag(cv))
  cv <- diag(rat)%*%cv%*%diag(rat)
  out <- matrix(NA,m,n-m)
  Xc <- X[z==0,]
  Xt <- X[z==1,]
  rownames(out) <- rownames(X)[z==1]
  colnames(out) <- rownames(X)[z==0]
  
  icov <- ginv(cv)
  if(m==1){
    out[1,] <- mahalanobis(Xc, Xt, icov, inverted=T)
  } else{
    for(i in 1:m){
      out[i,] <- mahalanobis(Xc, Xt[i,], icov, inverted=T)
    }
  }
  
  
  return(out)
}




addcaliper <- function(dmat, z, logitp, calipersd=.2, penalty=1000){
  sd.logitp <- sd(logitp)
  adif <- abs(outer(logitp[z==1],logitp[z==0],"-"))
  adif <- (adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
  dmat <- dmat+adif*penalty
  return(dmat)
}

Xmat <- (propscore.model$x[,-1])[-c(remove.trt.index, remove.ctrl.index),]
z <- x_data$AD[-c(remove.trt.index, remove.ctrl.index)]
  
# Rank based Mahalanobis distance
distmat <- rankmahal(z,Xmat)
# Add caliper
logit.propscore <- (predict(propscore.model))[-c(remove.trt.index, remove.ctrl.index)]
distmat2 <- addcaliper(distmat, z, logit.propscore) # 184 x 222


'''
#######################
## paired matching
#######################
matchvec.paired <- pairmatch(distmat2)

summary(matchvec.paired)
stratumStructure(matchvec.paired)

# assess covariate balance
treated.subject.index.paired <- rep(0,sum(z==1))
matched.control.subject.index.paired <- rep(0,length(treated.subject.index.paired))
matchedset.index.paired <- substr(matchvec.paired, start=3, stop=10) # group name
matchedset.index.paired.numeric <- as.numeric(matchedset.index.paired) # group name (numeric)
subjects.match.order.paired <- as.numeric(names(matchvec.paired)) # The subject indices 

# the order of matchvec
for(i in 1:length(treated.subject.index.paired)){
  matched.set.temp <- which(matchedset.index.paired.numeric==i) # ith matched group
  matched.set.temp.indices <- subjects.match.order.paired[matched.set.temp] # matched data indices within the entire data
  if(z[matched.set.temp.indices[1]]==1){ # if the first element is treated
    treated.subject.index.paired[i] <- matched.set.temp.indices[1] # treated subject index for the ith matched group
    matched.control.subject.index.paired[i] <- matched.set.temp.indices[2] # controlled subject index for the ith matched group
  }
  if(z[matched.set.temp.indices[2]]==1){ # if the second element is treated
    treated.subject.index.paired[i] <- matched.set.temp.indices[2]
    matched.control.subject.index.paired[i] <- matched.set.temp.indices[1]
  }
}


# Standardized differences before paired matching
controlmean.before <- apply(Xmat[z==0,], 2, mean)
treatmean <- apply(Xmat[z==1,], 2, mean)
treatvar <- apply(Xmat[z==1,], 2, var)
controlvar <- apply(Xmat[z==0,], 2, var)

stand.diff.before <- (treatmean - controlmean.before) / sqrt((treatvar + controlvar)/2)


# Standardized differences after paired matching
controlmean.after <- apply(Xmat[matched.control.subject.index.paired,], 2 ,mean)

# Standardized differences after matching
stand.diff.after.paired <- (treatmean - controlmean.after) / sqrt((treatvar + controlvar)/2)

pairedmatch.std.diff <- cbind(stand.diff.before, stand.diff.after.paired)
pairedmatch.std.diff


####################################
# full matching w/o restriction
####################################
matchvec.full <- fullmatch(distmat2)

treated.subject.index.full <- rep(0,sum(z==1))
# The subject indices in the order of matchvec
matchedset.index.full <- substr(matchvec.full, start=3, stop=10) # group name
matchedset.index.full.numeric <- as.numeric(matchedset.index.full) # group name (numeric)
subjects.match.order.full <- as.numeric(names(matchvec.full)) # subject indices

# Create a numeric variable for which stratum each unit belongs to
# 0 denotes that the unit was not matched
# there is not units not matched b/c we use full matching
stratum.short.full <- substr(matchvec.full, start=3, stop=10)
stratum.full.numeric <- as.numeric(stratum.short.full) # i-th element : group name of unit i

# Reassign numbers to each stratum that go from 1,..., no. of straum
sort.unique.stratum.full <- sort(unique(stratum.full.numeric)) # 136 groups
stratum.myindex.matchvecorder.full <- rep(0,length(stratum.full.numeric))
for(i in 1:length(sort.unique.stratum.full)){
  stratum.myindex.matchvecorder.full[stratum.full.numeric==sort.unique.stratum.full[i]] <- i  # 136 matched group sense indices / i-th element : group name (1~136) of unit i
}

stratum.myindex.full <- rep(0, length(stratum.myindex.matchvecorder.full)) # length 406
stratum.myindex.full[subjects.match.order.full] <- stratum.myindex.matchvecorder.full # i-th element : matched group index for unit i (1~136)

stratumStructure(matchvec.full)
# 8:1 6:1 5:1 4:1 3:1 2:1 1:1 1:2 1:3 1:4 1:5 1:6 1:7 
# 1   1   1   3   6  11  78  10  11   7   3   3   1 


# Calculate standardized difference before and after a full match
# Drop observations with missing values from the calculations
# stratum.myindex should contain strata for each subject with numbers 1:I (=
# number of strata), 0 means a unit was not matched
standardized.diff.func <- function(x, treatment, stratum.myindex, missing=rep(0,length(x))){
  xtreated <- x[treatment==1 & missing==0]
  xcontrol <- x[treatment==0 & missing==0]
  var.xtreated <- var(xtreated)
  var.xcontrol <- var(xcontrol)
  combinedsd <- sqrt(.5*(var.xtreated + var.xcontrol))
  std.diff.before.matching <- (mean(xtreated) - mean(xcontrol))/combinedsd
  nostratum <- length(unique(stratum.myindex)) - 1*min(stratum.myindex==0)
  diff.in.stratum <- rep(0, nostratum)
  treated.in.stratum <- rep(0,nostratum)
  
  for(i in 1:nostratum){
    diff.in.stratum[i] <- mean(x[stratum.myindex==i & treatment==1 & missing==0]) - mean(x[stratum.myindex==i & treatment==0 & missing==0])
    treated.in.stratum[i] <- sum(stratum.myindex==i & treatment==1 & missing==0)
    if(sum(stratum.myindex==i & treatment==0 & missing==0)==0){
      treated.in.stratum[i] <- 0
      diff.in.stratum[i] <- 0
    }
  }
  std.diff.after.matching <- (sum(treated.in.stratum*diff.in.stratum)/sum(treated.in.stratum))/combinedsd
  return(list(std.diff.before.matching = std.diff.before.matching, std.diff.after.matching = std.diff.after.matching))
}


# Calculate the standardized differences
Xmat.matchorder.full <- Xmat[subjects.match.order.full,]

std.diff.before <- rep(0, ncol(Xmat.matchorder.full))
std.diff.after.full <- rep(0, ncol(Xmat.matchorder.full))
names(std.diff.before) <- names(Xmat.matchorder.full[1,])
names(std.diff.after.full) <- names(Xmat.matchorder.full[1,])

for(i in 1:ncol(Xmat.matchorder.full)){
  temp.stand.diff <- standardized.diff.func(Xmat[,i], z, stratum.myindex.full)
  std.diff.before[i] <- temp.stand.diff$std.diff.before.matching
  std.diff.after.full[i] <- temp.stand.diff$std.diff.after.matching
}

fullmatch.std.diff <- cbind(std.diff.before, std.diff.after.full)
fullmatch.std.diff


#############################################################
## full matching with max.control = 2, min.control = 1/2
#############################################################
matchvec.full2 <- fullmatch(distmat2, max.controls = 2, min.controls = 1/2)

treated.subject.index.full2 <- rep(0,sum(z==1))
# The subject indices in the order of matchvec
matchedset.index.full2 <- substr(matchvec.full2, start=3, stop=10) # group name
matchedset.index.full2.numeric <- as.numeric(matchedset.index.full2) # group name (numeric)
subjects.match.order.full2 <- as.numeric(names(matchvec.full2)) # subject indices

# Create a numeric variable for which stratum each unit belongs to
# 0 denotes that the unit was not matched
# there is not units not matched b/c we use full matching
stratum.short.full2 <- substr(matchvec.full2, start=3, stop=10)
stratum.full2.numeric <- as.numeric(stratum.short.full2) # i-th element : group name of unit i

# Reassign numbers to each stratum that go from 1,..., no. of straum
sort.unique.stratum.full2 <- sort(unique(stratum.full2.numeric)) # 148 groups
stratum.myindex.matchvecorder.full2 <- rep(0,length(stratum.full2.numeric))
for(i in 1:length(sort.unique.stratum.full2)){
  stratum.myindex.matchvecorder.full2[stratum.full2.numeric==sort.unique.stratum.full2[i]] <- i  # 148 matched group sense indices / i-th element : group name (1~148) of unit i
}

stratum.myindex.full2 <- rep(0, length(stratum.myindex.matchvecorder.full2)) # length 406
stratum.myindex.full2[subjects.match.order.full2] <- stratum.myindex.matchvecorder.full2 # i-th element : matched group index for unit i (1~148)

stratumStructure(matchvec.full2)
# 2:1 1:1 1:2 
# 36  38  74 


# Calculate the standardized differences
Xmat.matchorder.full2 <- Xmat[subjects.match.order.full2,]

std.diff.before <- rep(0, ncol(Xmat.matchorder.full2))
std.diff.after.full2 <- rep(0, ncol(Xmat.matchorder.full2))
names(std.diff.before) <- names(Xmat.matchorder.full2[1,])
names(std.diff.after.full2) <- names(Xmat.matchorder.full2[1,])

for(i in 1:ncol(Xmat.matchorder.full2)){
  temp.stand.diff <- standardized.diff.func(Xmat[,i], z, stratum.myindex.full2)
  std.diff.before[i] <- temp.stand.diff$std.diff.before.matching
  std.diff.after.full2[i] <- temp.stand.diff$std.diff.after.matching
}

fullmatch.std.diff2 <- cbind(std.diff.before, std.diff.after.full2)
fullmatch.std.diff2


#############################################################
## full matching with max.control = 3, min.control = 1/3
#############################################################
matchvec.full3 <- fullmatch(distmat2, max.controls = 3, min.controls = 1/3)

treated.subject.index.full3 <- rep(0,sum(z==1))
# The subject indices in the order of matchvec
matchedset.index.full3 <- substr(matchvec.full3, start=3, stop=10) # group name
matchedset.index.full3.numeric <- as.numeric(matchedset.index.full3) # group name (numeric)
subjects.match.order.full3 <- as.numeric(names(matchvec.full3)) # subject indices

# Create a numeric variable for which stratum each unit belongs to
# 0 denotes that the unit was not matched
# there is not units not matched b/c we use full matching
stratum.short.full3 <- substr(matchvec.full3, start=3, stop=10)
stratum.full3.numeric <- as.numeric(stratum.short.full3) # i-th element : group name of unit i

# Reassign numbers to each stratum that go from 1,..., no. of straum
sort.unique.stratum.full3 <- sort(unique(stratum.full3.numeric)) # 141 groups
stratum.myindex.matchvecorder.full3 <- rep(0,length(stratum.full3.numeric))
for(i in 1:length(sort.unique.stratum.full3)){
  stratum.myindex.matchvecorder.full3[stratum.full3.numeric==sort.unique.stratum.full3[i]] <- i  # 141 matched group sense indices / i-th element : group name (1~141) of unit i
}

stratum.myindex.full3 <- rep(0, length(stratum.myindex.matchvecorder.full3)) # length 406
stratum.myindex.full3[subjects.match.order.full3] <- stratum.myindex.matchvecorder.full3 # i-th element : matched group index for unit i (1~141)

stratumStructure(matchvec.full3)
# 3:1 2:1 1:1 1:2 1:3 
# 15  13  66  13  34 


# Calculate the standardized differences
Xmat.matchorder.full3 <- Xmat[subjects.match.order.full3,]

std.diff.before <- rep(0, ncol(Xmat.matchorder.full3))
std.diff.after.full3 <- rep(0, ncol(Xmat.matchorder.full3))
names(std.diff.before) <- names(Xmat.matchorder.full3[1,])
names(std.diff.after.full3) <- names(Xmat.matchorder.full3[1,])

for(i in 1:ncol(Xmat.matchorder.full3)){
  temp.stand.diff <- standardized.diff.func(Xmat[,i], z, stratum.myindex.full3)
  std.diff.before[i] <- temp.stand.diff$std.diff.before.matching
  std.diff.after.full3[i] <- temp.stand.diff$std.diff.after.matching
}

fullmatch.std.diff3 <- cbind(std.diff.before, std.diff.after.full3)
fullmatch.std.diff3


#############################################################
## full matching with max.control = 4, min.control = 1/4
#############################################################
matchvec.full4 <- fullmatch(distmat2, max.controls = 4, min.controls = 1/4)

treated.subject.index.full4 <- rep(0,sum(z==1))
# The subject indices in the order of matchvec
matchedset.index.full4 <- substr(matchvec.full4, start=3, stop=10) # group name
matchedset.index.full4.numeric <- as.numeric(matchedset.index.full4) # group name (numeric)
subjects.match.order.full4 <- as.numeric(names(matchvec.full4)) # subject indices

# Create a numeric variable for which stratum each unit belongs to
# 0 denotes that the unit was not matched
# there is not units not matched b/c we use full matching
stratum.short.full4 <- substr(matchvec.full4, start=3, stop=10)
stratum.full4.numeric <- as.numeric(stratum.short.full4) # i-th element : group name of unit i

# Reassign numbers to each stratum that go from 1,..., no. of straum
sort.unique.stratum.full4 <- sort(unique(stratum.full4.numeric)) # 138 groups
stratum.myindex.matchvecorder.full4 <- rep(0,length(stratum.full4.numeric))
for(i in 1:length(sort.unique.stratum.full4)){
  stratum.myindex.matchvecorder.full4[stratum.full4.numeric==sort.unique.stratum.full4[i]] <- i  # 138 matched group sense indices / i-th element : group name (1~138) of unit i
}

stratum.myindex.full4 <- rep(0, length(stratum.myindex.matchvecorder.full4)) # length 406
stratum.myindex.full4[subjects.match.order.full4] <- stratum.myindex.matchvecorder.full4 # i-th element : matched group index for unit i (1~138)

stratumStructure(matchvec.full4)
# 4:1 3:1 2:1 1:1 1:2 1:3 1:4 
# 8   6  10  75  10  13  16 


# Calculate the standardized differences
Xmat.matchorder.full4 <- Xmat[subjects.match.order.full4,]

std.diff.before <- rep(0, ncol(Xmat.matchorder.full4))
std.diff.after.full4 <- rep(0, ncol(Xmat.matchorder.full4))
names(std.diff.before) <- names(Xmat.matchorder.full4[1,])
names(std.diff.after.full4) <- names(Xmat.matchorder.full4[1,])

for(i in 1:ncol(Xmat.matchorder.full4)){
  temp.stand.diff <- standardized.diff.func(Xmat[,i], z, stratum.myindex.full4)
  std.diff.before[i] <- temp.stand.diff$std.diff.before.matching
  std.diff.after.full4[i] <- temp.stand.diff$std.diff.after.matching
}

fullmatch.std.diff4 <- cbind(std.diff.before, std.diff.after.full4)
fullmatch.std.diff4


#########################
## almost exact matching
#########################
max.distmat2 <- max(distmat2)
diff <- abs(outer(Xmat[z==1,8], Xmat[z==0,8], "-")) # 184 x 222
distmat3 <- distmat2 + diff*(max.distmat2 + 10)


# paired matching
matchvec.paired.almostexact <- pairmatch(distmat3)

summary(matchvec.paired.almostexact)


# assess covariate balance
treated.subject.index.paired.almostexact <- rep(0,sum(z==1))
matched.control.subject.index.paired.almostexact <- rep(0,length(treated.subject.index.paired.almostexact))
matchedset.index.paired.almostexact <- substr(matchvec.paired.almostexact, start=3, stop=10) # group name
matchedset.index.paired.almostexact.numeric <- as.numeric(matchedset.index.paired.almostexact) # group name (numeric)
subjects.match.order.paired.almostexact <- as.numeric(names(matchvec.paired.almostexact)) # The subject indices 

# the order of matchvec
for(i in 1:length(treated.subject.index.paired.almostexact)){
  matched.set.temp <- which(matchedset.index.paired.almostexact.numeric==i) # ith matched group
  matched.set.temp.indices <- subjects.match.order.paired.almostexact[matched.set.temp] # matched data indices within the entire data
  if(z[matched.set.temp.indices[1]]==1){ # if the first element is treated
    treated.subject.index.paired.almostexact[i] <- matched.set.temp.indices[1] # treated subject index for the ith matched group
    matched.control.subject.index.paired.almostexact[i] <- matched.set.temp.indices[2] # controlled subject index for the ith matched group
  }
  if(z[matched.set.temp.indices[2]]==1){ # if the second element is treated
    treated.subject.index.paired.almostexact[i] <- matched.set.temp.indices[2]
    matched.control.subject.index.paired.almostexact[i] <- matched.set.temp.indices[1]
  }
}


# Standardized differences before paired.almostexact matching
controlmean.before <- apply(Xmat[z==0,], 2, mean)
treatmean <- apply(Xmat[z==1,], 2, mean)
treatvar <- apply(Xmat[z==1,], 2, var)
controlvar <- apply(Xmat[z==0,], 2, var)

stand.diff.before <- (treatmean - controlmean.before) / sqrt((treatvar + controlvar)/2)


# Standardized differences after paired.almostexact matching
controlmean.after <- apply(Xmat[matched.control.subject.index.paired.almostexact,], 2 ,mean)

# Standardized differences after matching
stand.diff.after.paired.almostexact <- (treatmean - controlmean.after) / sqrt((treatvar + controlvar)/2)

paired.almostexactmatch.std.diff <- cbind(stand.diff.before, stand.diff.after.paired.almostexact)
paired.almostexactmatch.std.diff




# full matching w/o restriction
matchvec.full.almostexact <- fullmatch(distmat3)

treated.subject.index.full.almostexact <- rep(0,sum(z==1))
# The subject indices in the order of matchvec
matchedset.index.full.almostexact <- substr(matchvec.full.almostexact, start=3, stop=10) # group name
matchedset.index.full.almostexact.numeric <- as.numeric(matchedset.index.full.almostexact) # group name (numeric)
subjects.match.order.full.almostexact <- as.numeric(names(matchvec.full.almostexact)) # subject indices

# Create a numeric variable for which stratum each unit belongs to
# 0 denotes that the unit was not matched
# there is not units not matched b/c we use full.almostexact matching
stratum.short.full.almostexact <- substr(matchvec.full.almostexact, start=3, stop=10)
stratum.full.almostexact.numeric <- as.numeric(stratum.short.full.almostexact) # i-th element : group name of unit i

# Reassign numbers to each stratum that go from 1,..., no. of straum
sort.unique.stratum.full.almostexact <- sort(unique(stratum.full.almostexact.numeric)) # 136 groups
stratum.myindex.matchvecorder.full.almostexact <- rep(0,length(stratum.full.almostexact.numeric))
for(i in 1:length(sort.unique.stratum.full.almostexact)){
  stratum.myindex.matchvecorder.full.almostexact[stratum.full.almostexact.numeric==sort.unique.stratum.full.almostexact[i]] <- i  # 136 matched group sense indices / i-th element : group name (1~136) of unit i
}

stratum.myindex.full.almostexact <- rep(0, length(stratum.myindex.matchvecorder.full.almostexact)) # length 406
stratum.myindex.full.almostexact[subjects.match.order.full.almostexact] <- stratum.myindex.matchvecorder.full.almostexact # i-th element : matched group index for unit i (1~136)

stratumStructure(matchvec.full.almostexact)



# Calculate the standardized differences
Xmat.matchorder.full.almostexact <- Xmat[subjects.match.order.full.almostexact,]

std.diff.before <- rep(0, ncol(Xmat.matchorder.full.almostexact))
std.diff.after.full.almostexact <- rep(0, ncol(Xmat.matchorder.full.almostexact))
names(std.diff.before) <- names(Xmat.matchorder.full.almostexact[1,])
names(std.diff.after.full.almostexact) <- names(Xmat.matchorder.full.almostexact[1,])

for(i in 1:ncol(Xmat.matchorder.full.almostexact)){
  temp.stand.diff <- standardized.diff.func(Xmat[,i], z, stratum.myindex.full.almostexact)
  std.diff.before[i] <- temp.stand.diff$std.diff.before.matching
  std.diff.after.full.almostexact[i] <- temp.stand.diff$std.diff.after.matching
}

full.almostexactmatch.std.diff <- cbind(std.diff.before, std.diff.after.full.almostexact)
full.almostexactmatch.std.diff



#########################
## exact matching
#########################
max.distmat2 <- max(distmat2)
diff2 <- abs(outer(Xmat[z==1,8], Xmat[z==0,8], "-")) # 184 x 222
diff2 <- ifelse(diff2!=0, 1, 0)
distmat4 <- distmat2 + diff*(max.distmat2 + 10000)

# paired matching
matchvec.paired.exact <- pairmatch(distmat4)

summary(matchvec.paired.exact)


# assess covariate balance
treated.subject.index.paired.exact <- rep(0,sum(z==1))
matched.control.subject.index.paired.exact <- rep(0,length(treated.subject.index.paired.exact))
matchedset.index.paired.exact <- substr(matchvec.paired.exact, start=3, stop=10) # group name
matchedset.index.paired.exact.numeric <- as.numeric(matchedset.index.paired.exact) # group name (numeric)
subjects.match.order.paired.exact <- as.numeric(names(matchvec.paired.exact)) # The subject indices 

# the order of matchvec
for(i in 1:length(treated.subject.index.paired.exact)){
  matched.set.temp <- which(matchedset.index.paired.exact.numeric==i) # ith matched group
  matched.set.temp.indices <- subjects.match.order.paired.exact[matched.set.temp] # matched data indices within the entire data
  if(z[matched.set.temp.indices[1]]==1){ # if the first element is treated
    treated.subject.index.paired.exact[i] <- matched.set.temp.indices[1] # treated subject index for the ith matched group
    matched.control.subject.index.paired.exact[i] <- matched.set.temp.indices[2] # controlled subject index for the ith matched group
  }
  if(z[matched.set.temp.indices[2]]==1){ # if the second element is treated
    treated.subject.index.paired.exact[i] <- matched.set.temp.indices[2]
    matched.control.subject.index.paired.exact[i] <- matched.set.temp.indices[1]
  }
}


# Standardized differences before paired.exact matching
controlmean.before <- apply(Xmat[z==0,], 2, mean)
treatmean <- apply(Xmat[z==1,], 2, mean)
treatvar <- apply(Xmat[z==1,], 2, var)
controlvar <- apply(Xmat[z==0,], 2, var)

stand.diff.before <- (treatmean - controlmean.before) / sqrt((treatvar + controlvar)/2)


# Standardized differences after paired.exact matching
controlmean.after <- apply(Xmat[matched.control.subject.index.paired.exact,], 2 ,mean)

# Standardized differences after matching
stand.diff.after.paired.exact <- (treatmean - controlmean.after) / sqrt((treatvar + controlvar)/2)

paired.exactmatch.std.diff <- cbind(stand.diff.before, stand.diff.after.paired.exact)
paired.exactmatch.std.diff




# full matching w/o restriction
matchvec.full.exact <- fullmatch(distmat4)

treated.subject.index.full.exact <- rep(0,sum(z==1))
# The subject indices in the order of matchvec
matchedset.index.full.exact <- substr(matchvec.full.exact, start=3, stop=10) # group name
matchedset.index.full.exact.numeric <- as.numeric(matchedset.index.full.exact) # group name (numeric)
subjects.match.order.full.exact <- as.numeric(names(matchvec.full.exact)) # subject indices

# Create a numeric variable for which stratum each unit belongs to
# 0 denotes that the unit was not matched
# there is not units not matched b/c we use full.exact matching
stratum.short.full.exact <- substr(matchvec.full.exact, start=3, stop=10)
stratum.full.exact.numeric <- as.numeric(stratum.short.full.exact) # i-th element : group name of unit i

# Reassign numbers to each stratum that go from 1,..., no. of straum
sort.unique.stratum.full.exact <- sort(unique(stratum.full.exact.numeric)) # 136 groups
stratum.myindex.matchvecorder.full.exact <- rep(0,length(stratum.full.exact.numeric))
for(i in 1:length(sort.unique.stratum.full.exact)){
  stratum.myindex.matchvecorder.full.exact[stratum.full.exact.numeric==sort.unique.stratum.full.exact[i]] <- i  # 136 matched group sense indices / i-th element : group name (1~136) of unit i
}

stratum.myindex.full.exact <- rep(0, length(stratum.myindex.matchvecorder.full.exact)) # length 406
stratum.myindex.full.exact[subjects.match.order.full.exact] <- stratum.myindex.matchvecorder.full.exact # i-th element : matched group index for unit i (1~136)

stratumStructure(matchvec.full.exact)



# Calculate the standardized differences
Xmat.matchorder.full.exact <- Xmat[subjects.match.order.full.exact,]

std.diff.before <- rep(0, ncol(Xmat.matchorder.full.exact))
std.diff.after.full.exact <- rep(0, ncol(Xmat.matchorder.full.exact))
names(std.diff.before) <- names(Xmat.matchorder.full.exact[1,])
names(std.diff.after.full.exact) <- names(Xmat.matchorder.full.exact[1,])

for(i in 1:ncol(Xmat.matchorder.full.exact)){
  temp.stand.diff <- standardized.diff.func(Xmat[,i], z, stratum.myindex.full.exact)
  std.diff.before[i] <- temp.stand.diff$std.diff.before.matching
  std.diff.after.full.exact[i] <- temp.stand.diff$std.diff.after.matching
}

full.exactmatch.std.diff <- cbind(std.diff.before, std.diff.after.full.exact)
full.exactmatch.std.diff
'''


################################################################
################################################################
###################### USE ALL DATA ############################
################################################################
################################################################
rankmahal <- function(z,X){
  X <- as.matrix(X)
  n <- dim(X)[1]
  rownames(X) <- 1:n
  k <- dim(X)[2]
  m <- sum(z)
  
  for(j in 1:k){
    X[,j] <- rank(X[,j])
  }
  cv <- cov(X)
  vuntied <- var(1:n)
  rat <- sqrt(vuntied/diag(cv))
  cv <- diag(rat)%*%cv%*%diag(rat)
  out <- matrix(NA,m,n-m)
  Xc <- X[z==0,]
  Xt <- X[z==1,]
  rownames(out) <- rownames(X)[z==1]
  colnames(out) <- rownames(X)[z==0]
  
  icov <- ginv(cv)
  for(i in 1:m){
    out[i,] <- mahalanobis(Xc, Xt[i,], icov, inverted=T)
  } 
  
  return(out)
}



addcaliper <- function(dmat, z, logitp, calipersd=.2, penalty=1000){
  sd.logitp <- sd(logitp)
  adif <- abs(outer(logitp[z==1],logitp[z==0],"-"))
  adif <- (adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
  dmat <- dmat+adif*penalty
  return(dmat)
}

Xmat.all <- (propscore.model$x[,-1])
z.all <- x_data$AD

# Rank based Mahalanobis distance
distmat.all <- rankmahal(z.all,Xmat.all)
# Add caliper
logit.propscore.all <- (predict(propscore.model))
distmat2.all <- addcaliper(distmat.all, z.all, logit.propscore.all) # 186 x 223


#######################
## paired matching
#######################
matchvec.paired.all <- pairmatch(distmat2.all)

summary(matchvec.paired.all)
stratumStructure(matchvec.paired.all)

# assess covariate balance
treated.subject.index.paired.all <- rep(0,sum(z.all==1))
matched.control.subject.index.paired.all <- rep(0,length(treated.subject.index.paired.all))
matchedset.index.paired.all <- substr(matchvec.paired.all, start=3, stop=10) # group name
matchedset.index.paired.numeric.all <- as.numeric(matchedset.index.paired.all) # group name (numeric)
subjects.match.order.paired.all <- as.numeric(names(matchvec.paired.all)) # The subject indices 

# the order of matchvec
for(i in 1:length(treated.subject.index.paired.all)){
  matched.set.temp <- which(matchedset.index.paired.numeric.all==i) # ith matched group
  matched.set.temp.indices <- subjects.match.order.paired.all[matched.set.temp] # matched data indices within the entire data
  if(z.all[matched.set.temp.indices[1]]==1){ # if the first element is treated
    treated.subject.index.paired.all[i] <- matched.set.temp.indices[1] # treated subject index for the ith matched group
    matched.control.subject.index.paired.all[i] <- matched.set.temp.indices[2] # controlled subject index for the ith matched group
  }
  if(z.all[matched.set.temp.indices[2]]==1){ # if the second element is treated
    treated.subject.index.paired.all[i] <- matched.set.temp.indices[2]
    matched.control.subject.index.paired.all[i] <- matched.set.temp.indices[1]
  }
}


# Standardized differences before paired matching
controlmean.before.all <- apply(Xmat.all[z.all==0,], 2, mean)
treatmean.all <- apply(Xmat.all[z.all==1,], 2, mean)
treatvar.all <- apply(Xmat.all[z.all==1,], 2, var)
controlvar.all <- apply(Xmat.all[z.all==0,], 2, var)

stand.diff.before.all <- (treatmean.all - controlmean.before.all) / sqrt((treatvar.all + controlvar.all)/2)


# Standardized differences after paired matching
controlmean.after.all <- apply(Xmat.all[matched.control.subject.index.paired.all,], 2 ,mean)

# Standardized differences after matching
stand.diff.after.paired.all <- (treatmean.all - controlmean.after.all) / sqrt((treatvar.all + controlvar.all)/2)

pairedmatch.std.diff.all <- cbind(stand.diff.before.all, stand.diff.after.paired.all)
pairedmatch.std.diff.all


####################################
# full matching w/o restriction
####################################
matchvec.full.all <- fullmatch(distmat2.all)

treated.subject.index.full.all <- rep(0,sum(z.all==1))
# The subject indices in the order of matchvec
matchedset.index.full.all <- substr(matchvec.full.all, start=3, stop=10) # group name
matchedset.index.full.numeric.all <- as.numeric(matchedset.index.full.all) # group name (numeric)
subjects.match.order.full.all <- as.numeric(names(matchvec.full.all)) # subject indices

# Create a numeric variable for which stratum each unit belongs to
# 0 denotes that the unit was not matched
# there is not units not matched b/c we use full matching
stratum.short.full.all <- substr(matchvec.full.all, start=3, stop=10)
stratum.full.numeric.all <- as.numeric(stratum.short.full.all) # i-th element : group name of unit i

# Reassign numbers to each stratum that go from 1,..., no. of straum
sort.unique.stratum.full.all <- sort(unique(stratum.full.numeric.all)) # 135 groups
stratum.myindex.matchvecorder.full.all <- rep(0,length(stratum.full.numeric.all))
for(i in 1:length(sort.unique.stratum.full.all)){
  stratum.myindex.matchvecorder.full.all[stratum.full.numeric.all==sort.unique.stratum.full.all[i]] <- i  # 135 matched group sense indices / i-th element : group name (1~135) of unit i
}

stratum.myindex.full.all <- rep(0, length(stratum.myindex.matchvecorder.full.all)) # length 409
stratum.myindex.full.all[subjects.match.order.full.all] <- stratum.myindex.matchvecorder.full.all # i-th element : matched group index for unit i (1~135)

stratumStructure(matchvec.full.all)
# 8:1 6:1 5:1 4:1 3:1 2:1 1:1 1:2 1:3 1:4 1:5 1:6 1:7 
# 1   1   1   3   7  12  74  11  11   7   3   2   2 


# Calculate standardized difference before and after a full match
# Drop observations with missing values from the calculations
# stratum.myindex should contain strata for each subject with numbers 1:I (=
# number of strata), 0 means a unit was not matched
standardized.diff.func <- function(x, treatment, stratum.myindex, missing=rep(0,length(x))){
  xtreated <- x[treatment==1 & missing==0]
  xcontrol <- x[treatment==0 & missing==0]
  var.xtreated <- var(xtreated)
  var.xcontrol <- var(xcontrol)
  combinedsd <- sqrt(.5*(var.xtreated + var.xcontrol))
  std.diff.before.matching <- (mean(xtreated) - mean(xcontrol))/combinedsd
  nostratum <- length(unique(stratum.myindex)) - 1*min(stratum.myindex==0)
  diff.in.stratum <- rep(0, nostratum)
  treated.in.stratum <- rep(0,nostratum)
  
  for(i in 1:nostratum){
    diff.in.stratum[i] <- mean(x[stratum.myindex==i & treatment==1 & missing==0]) - mean(x[stratum.myindex==i & treatment==0 & missing==0])
    treated.in.stratum[i] <- sum(stratum.myindex==i & treatment==1 & missing==0)
    if(sum(stratum.myindex==i & treatment==0 & missing==0)==0){
      treated.in.stratum[i] <- 0
      diff.in.stratum[i] <- 0
    }
  }
  std.diff.after.matching <- (sum(treated.in.stratum*diff.in.stratum)/sum(treated.in.stratum))/combinedsd
  return(list(std.diff.before.matching = std.diff.before.matching, std.diff.after.matching = std.diff.after.matching))
}


# Calculate the standardized differences
Xmat.matchorder.full.all <- Xmat.all[subjects.match.order.full.all,]

std.diff.before.all <- rep(0, ncol(Xmat.matchorder.full.all))
std.diff.after.full.all <- rep(0, ncol(Xmat.matchorder.full.all))
names(std.diff.before.all) <- names(Xmat.matchorder.full.all[1,])
names(std.diff.after.full.all) <- names(Xmat.matchorder.full.all[1,])

for(i in 1:ncol(Xmat.matchorder.full.all)){
  temp.stand.diff <- standardized.diff.func(Xmat.all[,i], z.all, stratum.myindex.full.all)
  std.diff.before.all[i] <- temp.stand.diff$std.diff.before.matching
  std.diff.after.full.all[i] <- temp.stand.diff$std.diff.after.matching
}

fullmatch.std.diff.all <- cbind(std.diff.before.all, std.diff.after.full.all)
fullmatch.std.diff.all


#############################################################
## full matching with max.control = 2, min.control = 1/2
#############################################################
matchvec.full2.all <- fullmatch(distmat2.all, max.controls = 2, min.controls = 1/2)

treated.subject.index.full2.all <- rep(0,sum(z.all==1))
# The subject indices in the order of matchvec
matchedset.index.full2.all <- substr(matchvec.full2.all, start=3, stop=10) # group name
matchedset.index.full2.numeric.all <- as.numeric(matchedset.index.full2.all) # group name (numeric)
subjects.match.order.full2.all <- as.numeric(names(matchvec.full2.all)) # subject indices

# Create a numeric variable for which stratum each unit belongs to
# 0 denotes that the unit was not matched
# there is not units not matched b/c we use full matching
stratum.short.full2.all <- substr(matchvec.full2.all, start=3, stop=10)
stratum.full2.numeric.all <- as.numeric(stratum.short.full2.all) # i-th element : group name of unit i

# Reassign numbers to each stratum that go from 1,..., no. of straum
sort.unique.stratum.full2.all <- sort(unique(stratum.full2.numeric.all)) # 148 groups
stratum.myindex.matchvecorder.full2.all <- rep(0,length(stratum.full2.numeric.all))
for(i in 1:length(sort.unique.stratum.full2.all)){
  stratum.myindex.matchvecorder.full2.all[stratum.full2.numeric.all==sort.unique.stratum.full2.all[i]] <- i  # 148 matched group sense indices / i-th element : group name (1~148) of unit i
}

stratum.myindex.full2.all <- rep(0, length(stratum.myindex.matchvecorder.full2.all)) # length 409
stratum.myindex.full2.all[subjects.match.order.full2.all] <- stratum.myindex.matchvecorder.full2.all # i-th element : matched group index for unit i (1~148)

stratumStructure(matchvec.full2.all)
# 2:1 1:1 1:2 
# 38  35  75 


# Calculate the standardized differences
Xmat.matchorder.full2.all <- Xmat.all[subjects.match.order.full2.all,]

std.diff.before.all <- rep(0, ncol(Xmat.matchorder.full2.all))
std.diff.after.full2.all <- rep(0, ncol(Xmat.matchorder.full2.all))
names(std.diff.before.all) <- names(Xmat.matchorder.full2.all[1,])
names(std.diff.after.full2.all) <- names(Xmat.matchorder.full2.all[1,])

for(i in 1:ncol(Xmat.matchorder.full2.all)){
  temp.stand.diff <- standardized.diff.func(Xmat.all[,i], z.all, stratum.myindex.full2.all)
  std.diff.before.all[i] <- temp.stand.diff$std.diff.before.matching
  std.diff.after.full2.all[i] <- temp.stand.diff$std.diff.after.matching
}

fullmatch.std.diff2.all <- cbind(std.diff.before.all, std.diff.after.full2.all)
fullmatch.std.diff2.all


#############################################################
## full matching with max.control = 3, min.control = 1/3
#############################################################
matchvec.full3.all <- fullmatch(distmat2.all, max.controls = 3, min.controls = 1/3)

treated.subject.index.full3.all <- rep(0,sum(z.all==1))
# The subject indices in the order of matchvec
matchedset.index.full3.all <- substr(matchvec.full3.all, start=3, stop=10) # group name
matchedset.index.full3.numeric.all <- as.numeric(matchedset.index.full3.all) # group name (numeric)
subjects.match.order.full3.all <- as.numeric(names(matchvec.full3.all)) # subject indices

# Create a numeric variable for which stratum each unit belongs to
# 0 denotes that the unit was not matched
# there is not units not matched b/c we use full matching
stratum.short.full3.all <- substr(matchvec.full3.all, start=3, stop=10)
stratum.full3.numeric.all <- as.numeric(stratum.short.full3.all) # i-th element : group name of unit i

# Reassign numbers to each stratum that go from 1,..., no. of straum
sort.unique.stratum.full3.all <- sort(unique(stratum.full3.numeric.all)) # 140 groups
stratum.myindex.matchvecorder.full3.all <- rep(0,length(stratum.full3.numeric.all))
for(i in 1:length(sort.unique.stratum.full3.all)){
  stratum.myindex.matchvecorder.full3.all[stratum.full3.numeric.all==sort.unique.stratum.full3.all[i]] <- i  # 140 matched group sense indices / i-th element : group name (1~140) of unit i
}

stratum.myindex.full3.all <- rep(0, length(stratum.myindex.matchvecorder.full3.all)) # length 409
stratum.myindex.full3.all[subjects.match.order.full3.all] <- stratum.myindex.matchvecorder.full3.all # i-th element : matched group index for unit i (1~140)

stratumStructure(matchvec.full3.all)
# 3:1 2:1 1:1 1:2 1:3 
# 16  14  63  11  36 


# Calculate the standardized differences
Xmat.matchorder.full3.all <- Xmat.all[subjects.match.order.full3.all,]

std.diff.before.all <- rep(0, ncol(Xmat.matchorder.full3.all))
std.diff.after.full3.all <- rep(0, ncol(Xmat.matchorder.full3.all))
names(std.diff.before.all) <- names(Xmat.matchorder.full3.all[1,])
names(std.diff.after.full3.all) <- names(Xmat.matchorder.full3.all[1,])

for(i in 1:ncol(Xmat.matchorder.full3.all)){
  temp.stand.diff <- standardized.diff.func(Xmat.all[,i], z.all, stratum.myindex.full3.all)
  std.diff.before.all[i] <- temp.stand.diff$std.diff.before.matching
  std.diff.after.full3.all[i] <- temp.stand.diff$std.diff.after.matching
}

fullmatch.std.diff3.all <- cbind(std.diff.before.all, std.diff.after.full3.all)
fullmatch.std.diff3.all


#############################################################
## full matching with max.control = 4, min.control = 1/4
#############################################################
matchvec.full4.all <- fullmatch(distmat2.all, max.controls = 4, min.controls = 1/4)

treated.subject.index.full4.all <- rep(0,sum(z.all==1))
# The subject indices in the order of matchvec
matchedset.index.full4.all <- substr(matchvec.full4.all, start=3, stop=10) # group name
matchedset.index.full4.numeric.all <- as.numeric(matchedset.index.full4.all) # group name (numeric)
subjects.match.order.full4.all <- as.numeric(names(matchvec.full4.all)) # subject indices

# Create a numeric variable for which stratum each unit belongs to
# 0 denotes that the unit was not matched
# there is not units not matched b/c we use full matching
stratum.short.full4.all <- substr(matchvec.full4.all, start=3, stop=10)
stratum.full4.numeric.all <- as.numeric(stratum.short.full4.all) # i-th element : group name of unit i

# Reassign numbers to each stratum that go from 1,..., no. of straum
sort.unique.stratum.full4.all <- sort(unique(stratum.full4.numeric.all)) # 137 groups
stratum.myindex.matchvecorder.full4.all <- rep(0,length(stratum.full4.numeric.all))
for(i in 1:length(sort.unique.stratum.full4.all)){
  stratum.myindex.matchvecorder.full4.all[stratum.full4.numeric.all==sort.unique.stratum.full4.all[i]] <- i  # 137 matched group sense indices / i-th element : group name (1~137) of unit i
}

stratum.myindex.full4.all <- rep(0, length(stratum.myindex.matchvecorder.full4.all)) # length 409
stratum.myindex.full4.all[subjects.match.order.full4.all] <- stratum.myindex.matchvecorder.full4.all # i-th element : matched group index for unit i (1~137)

stratumStructure(matchvec.full4.all)
# 4:1 3:1 2:1 1:1 1:2 1:3 1:4 
# 7   7  14  69  11  12  17 


# Calculate the standardized differences
Xmat.matchorder.full4.all <- Xmat.all[subjects.match.order.full4.all,]

std.diff.before.all <- rep(0, ncol(Xmat.matchorder.full4.all))
std.diff.after.full4.all <- rep(0, ncol(Xmat.matchorder.full4.all))
names(std.diff.before.all) <- names(Xmat.matchorder.full4.all[1,])
names(std.diff.after.full4.all) <- names(Xmat.matchorder.full4.all[1,])

for(i in 1:ncol(Xmat.matchorder.full4.all)){
  temp.stand.diff <- standardized.diff.func(Xmat.all[,i], z.all, stratum.myindex.full4.all)
  std.diff.before.all[i] <- temp.stand.diff$std.diff.before.matching
  std.diff.after.full4.all[i] <- temp.stand.diff$std.diff.after.matching
}

fullmatch.std.diff4.all <- cbind(std.diff.before.all, std.diff.after.full4.all)
fullmatch.std.diff4.all


#########################
## almost exact matching
#########################
max.distmat2.all <- max(distmat2.all)
diff.all <- abs(outer(Xmat.all[z.all==1,8], Xmat.all[z.all==0,8], "-")) # 186 x 223
distmat3.all <- distmat2.all + diff.all*(max.distmat2.all + 10)


# paired matching
matchvec.paired.almostexact.all <- pairmatch(distmat3.all)

summary(matchvec.paired.almostexact.all)


# assess covariate balance
treated.subject.index.paired.almostexact.all <- rep(0,sum(z.all==1))
matched.control.subject.index.paired.almostexact.all <- rep(0,length(treated.subject.index.paired.almostexact.all))
matchedset.index.paired.almostexact.all <- substr(matchvec.paired.almostexact.all, start=3, stop=10) # group name
matchedset.index.paired.almostexact.numeric.all <- as.numeric(matchedset.index.paired.almostexact.all) # group name (numeric)
subjects.match.order.paired.almostexact.all <- as.numeric(names(matchvec.paired.almostexact.all)) # The subject indices 

# the order of matchvec
for(i in 1:length(treated.subject.index.paired.almostexact.all)){
  matched.set.temp <- which(matchedset.index.paired.almostexact.numeric.all==i) # ith matched group
  matched.set.temp.indices <- subjects.match.order.paired.almostexact.all[matched.set.temp] # matched data indices within the entire data
  if(z.all[matched.set.temp.indices[1]]==1){ # if the first element is treated
    treated.subject.index.paired.almostexact.all[i] <- matched.set.temp.indices[1] # treated subject index for the ith matched group
    matched.control.subject.index.paired.almostexact.all[i] <- matched.set.temp.indices[2] # controlled subject index for the ith matched group
  }
  if(z.all[matched.set.temp.indices[2]]==1){ # if the second element is treated
    treated.subject.index.paired.almostexact.all[i] <- matched.set.temp.indices[2]
    matched.control.subject.index.paired.almostexact.all[i] <- matched.set.temp.indices[1]
  }
}


# Standardized differences before paired.almostexact matching
controlmean.before.all <- apply(Xmat.all[z.all==0,], 2, mean)
treatmean.all <- apply(Xmat.all[z.all==1,], 2, mean)
treatvar.all <- apply(Xmat.all[z.all==1,], 2, var)
controlvar.all <- apply(Xmat.all[z.all==0,], 2, var)

stand.diff.before.all <- (treatmean.all - controlmean.before.all) / sqrt((treatvar.all + controlvar.all)/2)


# Standardized differences after paired.almostexact matching
controlmean.after.all <- apply(Xmat.all[matched.control.subject.index.paired.almostexact.all,], 2 ,mean)

# Standardized differences after matching
stand.diff.after.paired.almostexact.all <- (treatmean.all - controlmean.after.all) / sqrt((treatvar.all + controlvar.all)/2)

paired.almostexactmatch.std.diff.all <- cbind(stand.diff.before.all, stand.diff.after.paired.almostexact.all)
paired.almostexactmatch.std.diff.all




# full matching w/o restriction
matchvec.full.almostexact.all <- fullmatch(distmat3.all)

treated.subject.index.full.almostexact.all <- rep(0,sum(z.all==1))
# The subject indices in the order of matchvec
matchedset.index.full.almostexact.all <- substr(matchvec.full.almostexact.all, start=3, stop=10) # group name
matchedset.index.full.almostexact.numeric.all <- as.numeric(matchedset.index.full.almostexact.all) # group name (numeric)
subjects.match.order.full.almostexact.all <- as.numeric(names(matchvec.full.almostexact.all)) # subject indices

# Create a numeric variable for which stratum each unit belongs to
# 0 denotes that the unit was not matched
# there is not units not matched b/c we use full.almostexact matching
stratum.short.full.almostexact.all <- substr(matchvec.full.almostexact.all, start=3, stop=10)
stratum.full.almostexact.numeric.all <- as.numeric(stratum.short.full.almostexact.all) # i-th element : group name of unit i

# Reassign numbers to each stratum that go from 1,..., no. of straum
sort.unique.stratum.full.almostexact.all <- sort(unique(stratum.full.almostexact.numeric.all)) # 117 groups
stratum.myindex.matchvecorder.full.almostexact.all <- rep(0,length(stratum.full.almostexact.numeric.all))
for(i in 1:length(sort.unique.stratum.full.almostexact.all)){
  stratum.myindex.matchvecorder.full.almostexact.all[stratum.full.almostexact.numeric.all==sort.unique.stratum.full.almostexact.all[i]] <- i  # 117 matched group sense indices / i-th element : group name (1~117) of unit i
}

stratum.myindex.full.almostexact.all <- rep(0, length(stratum.myindex.matchvecorder.full.almostexact.all)) # length 409
stratum.myindex.full.almostexact.all[subjects.match.order.full.almostexact.all] <- stratum.myindex.matchvecorder.full.almostexact.all # i-th element : matched group index for unit i (1~117)

stratumStructure(matchvec.full.almostexact.all)
# 8:1  6:1  5:1  4:1  3:1  2:1  1:1  1:2  1:3  1:4  1:5  1:6  1:7  1:9 1:12 
# 2    2    4    5    2   10   54   11   11    7    4    1    2    1    1


# Calculate the standardized differences
Xmat.matchorder.full.almostexact.all <- Xmat.all[subjects.match.order.full.almostexact.all,]

std.diff.before.all <- rep(0, ncol(Xmat.matchorder.full.almostexact.all))
std.diff.after.full.almostexact.all <- rep(0, ncol(Xmat.matchorder.full.almostexact.all))
names(std.diff.before.all) <- names(Xmat.matchorder.full.almostexact.all[1,])
names(std.diff.after.full.almostexact.all) <- names(Xmat.matchorder.full.almostexact.all[1,])

for(i in 1:ncol(Xmat.matchorder.full.almostexact.all)){
  temp.stand.diff <- standardized.diff.func(Xmat.all[,i], z.all, stratum.myindex.full.almostexact.all)
  std.diff.before.all[i] <- temp.stand.diff$std.diff.before.matching
  std.diff.after.full.almostexact.all[i] <- temp.stand.diff$std.diff.after.matching
}

full.almostexactmatch.std.diff.all <- cbind(std.diff.before.all, std.diff.after.full.almostexact.all)
full.almostexactmatch.std.diff.all



#########################
## exact matching
#########################
max.distmat2.all <- max(distmat2.all)
diff2.all <- abs(outer(Xmat.all[z.all==1,8], Xmat.all[z.all==0,8], "-")) # 186 x 223
diff2.all <- ifelse(diff2.all!=0, 1, 0)
distmat4.all <- distmat2.all + diff.all*(max.distmat2.all + 10000)

# paired matching
matchvec.paired.exact.all <- pairmatch(distmat4.all)

summary(matchvec.paired.exact.all)


# assess covariate balance
treated.subject.index.paired.exact.all <- rep(0,sum(z.all==1))
matched.control.subject.index.paired.exact.all <- rep(0,length(treated.subject.index.paired.exact.all))
matchedset.index.paired.exact.all <- substr(matchvec.paired.exact.all, start=3, stop=10) # group name
matchedset.index.paired.exact.numeric.all <- as.numeric(matchedset.index.paired.exact.all) # group name (numeric)
subjects.match.order.paired.exact.all <- as.numeric(names(matchvec.paired.exact.all)) # The subject indices 

# the order of matchvec
for(i in 1:length(treated.subject.index.paired.exact.all)){
  matched.set.temp <- which(matchedset.index.paired.exact.numeric.all==i) # ith matched group
  matched.set.temp.indices <- subjects.match.order.paired.exact.all[matched.set.temp] # matched data indices within the entire data
  if(z.all[matched.set.temp.indices[1]]==1){ # if the first element is treated
    treated.subject.index.paired.exact.all[i] <- matched.set.temp.indices[1] # treated subject index for the ith matched group
    matched.control.subject.index.paired.exact.all[i] <- matched.set.temp.indices[2] # controlled subject index for the ith matched group
  }
  if(z.all[matched.set.temp.indices[2]]==1){ # if the second element is treated
    treated.subject.index.paired.exact.all[i] <- matched.set.temp.indices[2]
    matched.control.subject.index.paired.exact.all[i] <- matched.set.temp.indices[1]
  }
}


# Standardized differences before paired.exact matching
controlmean.before.all <- apply(Xmat.all[z.all==0,], 2, mean)
treatmean.all <- apply(Xmat.all[z.all==1,], 2, mean)
treatvar.all <- apply(Xmat.all[z.all==1,], 2, var)
controlvar.all <- apply(Xmat.all[z.all==0,], 2, var)

stand.diff.before.all <- (treatmean.all - controlmean.before.all) / sqrt((treatvar.all + controlvar.all)/2)


# Standardized differences after paired.exact matching
controlmean.after.all <- apply(Xmat.all[matched.control.subject.index.paired.exact.all,], 2 ,mean)

# Standardized differences after matching
stand.diff.after.paired.exact.all <- (treatmean.all - controlmean.after.all) / sqrt((treatvar.all + controlvar.all)/2)

paired.exactmatch.std.diff.all <- cbind(stand.diff.before.all, stand.diff.after.paired.exact.all)
paired.exactmatch.std.diff.all




# full matching w/o restriction
matchvec.full.exact.all <- fullmatch(distmat4.all)

treated.subject.index.full.exact.all <- rep(0,sum(z.all==1))
# The subject indices in the order of matchvec
matchedset.index.full.exact.all <- substr(matchvec.full.exact.all, start=3, stop=10) # group name
matchedset.index.full.exact.numeric.all <- as.numeric(matchedset.index.full.exact.all) # group name (numeric)
subjects.match.order.full.exact.all <- as.numeric(names(matchvec.full.exact.all)) # subject indices

# Create a numeric variable for which stratum each unit belongs to
# 0 denotes that the unit was not matched
# there is not units not matched b/c we use full.exact matching
stratum.short.full.exact.all <- substr(matchvec.full.exact.all, start=3, stop=10)
stratum.full.exact.numeric.all <- as.numeric(stratum.short.full.exact.all) # i-th element : group name of unit i

# Reassign numbers to each stratum that go from 1,..., no. of straum
sort.unique.stratum.full.exact.all <- sort(unique(stratum.full.exact.numeric.all)) # 117 groups
stratum.myindex.matchvecorder.full.exact.all <- rep(0,length(stratum.full.exact.numeric.all))
for(i in 1:length(sort.unique.stratum.full.exact.all)){
  stratum.myindex.matchvecorder.full.exact.all[stratum.full.exact.numeric.all==sort.unique.stratum.full.exact.all[i]] <- i  # 117 matched group sense indices / i-th element : group name (1~117) of unit i
}

stratum.myindex.full.exact.all <- rep(0, length(stratum.myindex.matchvecorder.full.exact.all)) # length 406
stratum.myindex.full.exact.all[subjects.match.order.full.exact.all] <- stratum.myindex.matchvecorder.full.exact.all # i-th element : matched group index for unit i (1~117)

stratumStructure(matchvec.full.exact.all)
# 8:1  6:1  5:1  4:1  3:1  2:1  1:1  1:2  1:3  1:4  1:5  1:6  1:7  1:9 1:12 
# 2    2    4    5    2   10   53   13   11    5    5    1    2    1    1 


# Calculate the standardized differences
Xmat.matchorder.full.exact.all <- Xmat.all[subjects.match.order.full.exact.all,]

std.diff.before.all <- rep(0, ncol(Xmat.matchorder.full.exact.all))
std.diff.after.full.exact.all <- rep(0, ncol(Xmat.matchorder.full.exact.all))
names(std.diff.before.all) <- names(Xmat.matchorder.full.exact.all[1,])
names(std.diff.after.full.exact.all) <- names(Xmat.matchorder.full.exact.all[1,])

for(i in 1:ncol(Xmat.matchorder.full.exact.all)){
  temp.stand.diff <- standardized.diff.func(Xmat.all[,i], z.all, stratum.myindex.full.exact.all)
  std.diff.before.all[i] <- temp.stand.diff$std.diff.before.matching
  std.diff.after.full.exact.all[i] <- temp.stand.diff$std.diff.after.matching
}

full.exactmatch.std.diff.all <- cbind(std.diff.before.all, std.diff.after.full.exact.all)
full.exactmatch.std.diff.all



######################################
## Experiments: corpus callosum data
######################################

source("efficiency.R")
source("GeodRegr.R")
source("normal.R")
source("others.R")

## initializations

manifold <- 'kendall'

boundary_points <- 50
dim <- 2 * boundary_points - 4
embed <- boundary_points
# k <- 1
# estimator <- 'l2'
# estimator <- 'l1'
# estimator <- 'tukey'

x_data2 <- t(t(read.table(file="data/x_data")))
y_data2 <- read.table(file="data/y_data")


calculateWeight <- function(z, strata, beta.T, beta.C){
  N <- length(z)
  L <- length(unique(strata))
  m.T <- as.numeric(table(strata[z==1]))
  m.C <- as.numeric(table(strata[z==0]))
  lambda.hat <- as.numeric(table(strata)) / N
  # beta.T <- sum(z) / N
  # beta.C <- 1-beta.T
  wt <- beta.T*((lambda.hat/m.T)[strata])*z + beta.C*((lambda.hat/m.C)[strata])*(1-z)
  return(wt)
}

calculateWeight.T <- function(z, strata){
  N <- length(z)
  L <- length(unique(strata))
  m.T <- as.numeric(table(strata[z==1]))
  lambda.hat <- as.numeric(table(strata)) / N
  wt <- ((lambda.hat/m.T)[strata])*z 
  return(wt)
}

calculateWeight.C <- function(z, strata){
  N <-length(z)
  L <- length(unique(strata))
  m.C <- as.numeric(table(strata[z==0]))
  lambda.hat <- as.numeric(table(strata)) / N
  wt <- ((lambda.hat/m.C)[strata])*(1-z)
  return(wt)
}

getpval <- function(z, strata, T.obs, B=1000, method, estimator, p_tol=1e-10, v_tol=1e-10, max_iter=200){
  N <- length(z)
  L <- length(unique(strata))
  strata.list <- list()
  z.list <- list()
  for(i in 1:L){
    strata.list[[i]] <- which(strata==i)
    z.list[[i]] <- z[which(strata==i)]
  }
  
  T.shuffle <- vector(length=B)
  if(method=="regression"){
    for(j in 1:B){
      print(paste("iter:", j))
      z.shuffle <- vector(length=N)
      z.shuffle.list <- lapply(z.list, function(x) x[sample(1:length(x))])
      z.shuffle[unlist(strata.list)] <- unlist(z.shuffle.list)
      
      tmp <- geo_reg(manifold, x_data2[,9], y_data2, w = calculateWeight(z.shuffle, strata, sum(z.shuffle)/N, 1-sum(z.shuffle)/N), 
                     estimator, max_iter = max_iter)
      T.shuffle[j] <- sqrt(abs(sum(Conj(tmp$V)*tmp$V)))
    }
  } else if(method=="dist"){
    for(j in 1:B){
      print(paste("iter:", j))
      z.shuffle <- vector(length=N)
      z.shuffle.list <- lapply(z.list, function(x) x[sample(1:length(x))])
      z.shuffle[unlist(strata.list)] <- unlist(z.shuffle.list)
      
      T.shuffle[j] <- geo_dist(manifold,intrinsic_location(manifold, y_data2, w = calculateWeight.T(z.shuffle, strata), estimator, p_tol=p_tol, V_tol=v_tol),
                               intrinsic_location(manifold, y_data2, w = calculateWeight.C(z.shuffle, strata), estimator, p_tol=p_tol, V_tol=v_tol))
    }
  }
  
  p.val <- mean(T.shuffle >= T.obs)
  
  res <- list()
  res$T.shuffle <- T.shuffle
  res$pval <- p.val
  return(res)  
}

getpval_with_x <- function(z, strata, T.obs, B=1000, method, estimator, p_tol=1e-10, v_tol=1e-10, max_iter=200){
  N <- length(z)
  L <- length(unique(strata))
  strata.list <- list()
  z.list <- list()
  for(i in 1:L){
    strata.list[[i]] <- which(strata==i)
    z.list[[i]] <- z[which(strata==i)]
  }
  
  T.shuffle <- vector(length=B)
  if(method=="regression"){
    for(j in 1:B){
      # print(paste("iter:", j))
      z.shuffle <- vector(length=N)
      z.shuffle.list <- lapply(z.list, function(x) x[sample(1:length(x))])
      z.shuffle[unlist(strata.list)] <- unlist(z.shuffle.list)
      
      x_data2.tmp <- x_data2
      x_data2.tmp[,9] <- z.shuffle
      tmp <- geo_reg(manifold, x_data2.tmp, y_data2, w = calculateWeight(z.shuffle, strata, sum(z.shuffle)/N, 1-sum(z.shuffle)/N), 
                     estimator, max_iter = max_iter)
      T.shuffle[j] <- sqrt(abs(sum(Conj(tmp$V)*tmp$V)))
      print(paste("B:", j, "iteration:", tmp$iteration))
    }
  } else if(method=="dist"){
    for(j in 1:B){
      print(paste("iter:", j))
      z.shuffle <- vector(length=N)
      z.shuffle.list <- lapply(z.list, function(x) x[sample(1:length(x))])
      z.shuffle[unlist(strata.list)] <- unlist(z.shuffle.list)
      
      T.shuffle[j] <- geo_dist(manifold,intrinsic_location(manifold, y_data2, w = calculateWeight.T(z.shuffle, strata), estimator, p_tol=p_tol, V_tol=v_tol),
                               intrinsic_location(manifold, y_data2, w = calculateWeight.C(z.shuffle, strata), estimator, p_tol=p_tol, V_tol=v_tol))
    }
  }
  
  p.val <- mean(T.shuffle >= T.obs)
  
  res <- list()
  res$T.shuffle <- T.shuffle
  res$pval <- p.val
  return(res)  
}


group.use <- stratum.myindex.full.all

estimator <- 'l2'

# distance 
T1.l2 <- geo_dist(manifold,intrinsic_location(manifold, y_data2, w = calculateWeight.T(z.all, group.use), estimator, p_tol=1e-6, V_tol=1e-6),
                  intrinsic_location(manifold, y_data2, w = calculateWeight.C(z.all, group.use), estimator, p_tol=1e-6, V_tol=1e-6))

write.csv(intrinsic_location(manifold, y_data2, w = calculateWeight.T(z.all, group.use), estimator, p_tol=1e-6, V_tol=1e-6), 
          file="res/T1_L2_trt.csv", row.names=FALSE)
write.csv(intrinsic_location(manifold, y_data2, w = calculateWeight.C(z.all, group.use), estimator, p_tol=1e-6, V_tol=1e-6), 
          file="res/T1_L2_ctrl.csv", row.names=FALSE)

# ttmp <- read.csv("T1_L2_trt.csv", header=TRUE)

# regression with only z
T2.res.l2 <- geo_reg(manifold, x_data2[,9], y_data2, w = calculateWeight(z.all, group.use, sum(z.all)/length(z.all), 1-sum(z.all)/length(z.all)), 
                     estimator, max_iter = 200)

T2.l2 <- sqrt(abs(sum(Conj(T2.res.l2$V)*T2.res.l2$V)))

# regression with z & x
# T3.res.l2 <- geo_reg(manifold, x_data2, y_data2, w = calculateWeight(z.all, group.use, sum(z.all)/length(z.all), 1-sum(z.all)/length(z.all)), estimator)
# 
# T3.l2 <- sqrt(abs(sum(Conj(T3.res.l2$V)*T3.res.l2$V)))


res.dist.l2 <- getpval(z.all, group.use, T1.l2, B=1000, method="dist", estimator = 'l2', p_tol=1e-6, v_tol=1e-6)
# res.reg.l2 <- getpval_with_x(z.all, group.use, T3.l2, B=1000, method="regression", estimator = 'l2', max_iter=200)


estimator <- 'l1'

# distance
T1.l1 <- geo_dist(manifold,intrinsic_location(manifold, y_data2, w = calculateWeight.T(z.all, group.use), estimator, p_tol=1e-6, V_tol=1e-6),
                  intrinsic_location(manifold, y_data2, w = calculateWeight.C(z.all, group.use), estimator, p_tol=1e-6, V_tol=1e-6))

write.csv(intrinsic_location(manifold, y_data2, w = calculateWeight.T(z.all, group.use), estimator, p_tol=1e-6, V_tol=1e-6), 
          file="res/T1_L1_trt.csv", row.names=FALSE)
write.csv(intrinsic_location(manifold, y_data2, w = calculateWeight.C(z.all, group.use), estimator, p_tol=1e-6, V_tol=1e-6), 
          file="res/T1_L1_ctrl.csv", row.names=FALSE)


T2.res.l1 <- geo_reg(manifold, x_data2[,9], y_data2, w = calculateWeight(z.all, group.use, sum(z.all)/length(z.all), 1-sum(z.all)/length(z.all)), estimator, p_tol=1e-6, V_tol=1e-6)

T2.l1 <- sqrt(abs(sum(Conj(T2.res.l1$V)*T2.res.l1$V)))

# T3.res.l1 <- geo_reg(manifold, x_data2, y_data2, w = calculateWeight(z.all, group.use, sum(z.all)/length(z.all), 1-sum(z.all)/length(z.all)), estimator)
# 
# T3.l1 <- sqrt(abs(sum(Conj(T3.res.l1$V)*T3.res.l1$V)))



res.dist.l1 <- getpval(z.all, group.use, T1.l1, B=1000, method="dist", estimator = 'l1', p_tol=1e-6, v_tol=1e-6)
# res.reg.l1 <- getpval_with_x(z.all, group.use, T3.l1, B=1000, method="regression", estimator = 'l1', max_iter=200)


# AATE result
hist(res.dist.l2$T.shuffle, main=expression(bold(paste("Permutation Test Distribution of AATE with", ~ L[2], ~"loss (B=1000)"))), 
     xlim = c(min(res.dist.l2$T.shuffle), T1.l2), xaxt = "n", xlab="AATE", breaks = seq(from=0, to=0.02, by=0.001))
axis(1, at = seq(0, 0.02, by=0.002))
# segments(x0=T1.l2, y0=0, x1=T1.l2, y1=180,col="red", lwd=2)
abline(v=T1.l2,col="red",lwd=2)
text(0.016, 180, expression(paste(T[paste(L[2],",", 1)], ~ "(P-value=0)")), pos = 3, srt = 0, cex=1.2)

# AMTE result
hist(res.dist.l1$T.shuffle, main=expression(bold(paste("Permutation Test Distribution of AMTE with", ~ L[1], ~"loss (B=1000)"))), 
     xlim = c(min(res.dist.l1$T.shuffle), T1.l1), xaxt = "n", xlab="AMTE", breaks=seq(from=0, to=0.02, by=0.001))
axis(1, at = seq(0, 0.02, by=0.002))
abline(v=T1.l1,col="red",lwd=2)
text(0.0166, 163, expression(paste(T[paste(L[1],",", 1)], ~ "(P-value=0)")), pos = 3, srt = 0, cex=1.2)


# N <- length(z.all)
# L <- length(unique(strata))
# B <- 1
# method <- "dist"
# p_tol <- 1e-10
# v_tol <- 1e-10
# strata.list <- list()
# z.list <- list()
# for(i in 1:L){
#   strata.list[[i]] <- which(strata==i)
#   z.list[[i]] <- z.all[which(strata==i)]
# }
# 
# T.shuffle <- vector(length=B)
# if(method=="regression"){
#   for(j in 1:B){
#     z.shuffle <- vector(length=N)
#     z.shuffle.list <- lapply(z.list, function(x) x[sample(1:length(x))])
#     z.shuffle[unlist(strata.list)] <- unlist(z.shuffle.list)
#     
#     T.shuffle[j] <- geo_reg(manifold, x_data2[,9], y_data2, w = calculateWeight(z.shuffle, strata, sum(z.shuffle)/N, 1-sum(z.shuffle)/N), 
#                             estimator, p_tol=p_tol, V_tol=v_tol)
#   }
# } else if(method=="dist"){
#   for(j in 1:B){
#     z.shuffle <- vector(length=N)
#     z.shuffle.list <- lapply(z.list, function(x) x[sample(1:length(x))])
#     z.shuffle[unlist(strata.list)] <- unlist(z.shuffle.list)
#     
#     T.shuffle[j] <- geo_dist(manifold,intrinsic_location(manifold, y_data2, w = calculateWeight.T(z.shuffle, strata), estimator, p_tol=p_tol, V_tol=v_tol),
#                              intrinsic_location(manifold, y_data2, w = calculateWeight.C(z.shuffle, strata), estimator, p_tol=p_tol, V_tol=v_tol))
#   }
# }
# 
# z.shuffle
# z.shuffle.list
# 
# p.val <- mean(T.shuffle >= T.obs)
# 
# res <- list()
# res$T.shuffle <- T.shuffle
# res$pval <- p.val



## Visualization ## 


# projecting data onto 2d shape space

thing <- list(ADF_y_data, ADM_y_data, NCF_y_data, NCM_y_data)
kss <- vector('list',4)

for (k in 1:4) {
  r_data <- thing[[k]]
  r_data <- aperm(r_data, c(2, 1, 3))
  
  for (i in 1:dim(r_data)[3]) { ## remove translation
    for (j in 1:2) {
      r_data[j, , i] <- r_data[j, , i] - mean(r_data[j, , i])
    }
  }
  
  for (i in 1:dim(r_data)[3]) { ## remove scaling
    r_data[, , i] <- r_data[, , i] / ((sum(r_data[, , i] * r_data[, , i])) ^ 0.5)
  }
  
  r_data <- r_data[1, , ] + r_data[2, , ] * (0 + 1i)
  kss[[k]] <- r_data
}

r_data <- cbind(kss[[1]], kss[[2]], kss[[3]], kss[[4]])


## corpus callosum illustrations

l2c <- as.matrix(read.csv('res/T1_L2_ctrl.csv'))
l2t <- as.matrix(read.csv('res/T1_L2_trt.csv'))
l1c <- as.matrix(read.csv('res/T1_L1_ctrl.csv'))
l1t <- as.matrix(read.csv('res/T1_L1_trt.csv'))

cc <- list(l2c,l2t,l1c,l1t)

for (i in 1:4) {
  a <- colSums(l2c*Conj(cc[[i]]))
  cc[[i]] <- cc[[i]]*(a/abs(a))
}

par(mfrow=c(1,1))

pal <- colorRampPalette(c("blue", "red"))(2)

plot(Re(cc[[1]]), Im(cc[[1]]), type='n', xaxt='n', yaxt='n', ann=FALSE, asp=1)

for (i in 1:2) {
  lines(Re(cc[[i]]), Im(cc[[i]]), col=pal[i])
  lines(c(Re(cc[[i]])[50,], Re(cc[[i]])[1,]), c(Im(cc[[i]])[50,], Im(cc[[i]])[1,]), col=pal[i])
}


plot(Re(cc[[3]]), Im(cc[[3]]), type='n', xaxt='n', yaxt='n', ann=FALSE, asp=1)
for (i in 3:4) {
  lines(Re(cc[[i]]), Im(cc[[i]]), col=pal[i-2])
  lines(c(Re(cc[[i]])[50,], Re(cc[[i]])[1,]), c(Im(cc[[i]])[50,], Im(cc[[i]])[1,]), col=pal[i-2])
}


# gray scale
ltys <- c(2,1)

plot(Re(cc[[1]]), Im(cc[[1]]), type='n', xaxt='n', yaxt='n', ann=FALSE, asp=1)

for (i in 1:2) {
  lines(Re(cc[[i]]), Im(cc[[i]]), lty=ltys[i])
  lines(c(Re(cc[[i]])[50,], Re(cc[[i]])[1,]), c(Im(cc[[i]])[50,], Im(cc[[i]])[1,]), lty=ltys[i])
}


plot(Re(cc[[3]]), Im(cc[[3]]), type='n', xaxt='n', yaxt='n', ann=FALSE, asp=1)
for (i in 3:4) {
  lines(Re(cc[[i]]), Im(cc[[i]]), lty=ltys[i-2])
  lines(c(Re(cc[[i]])[50,], Re(cc[[i]])[1,]), c(Im(cc[[i]])[50,], Im(cc[[i]])[1,]), lty=ltys[i-2])
}



## confidence interval 
M <- 1000 # number of bootstrap samples
N <- nrow(x_data) # number of data
T2.l2.bstrp.list <- c()
T2.l1.bstrp.list <- c()
set.seed(1)
for(i in 1:M){
  print(i)
  
  bstrp.sample <- sort(sample(1:N, N, replace=TRUE))
  x_data.bstrp <- x_data[bstrp.sample,]
  
  propscore.model.bstrp <- glm(AD ~ gender + handedness + marital1 + marital2 + marital3 + 
                           educationlength + retirement + age,
                         family=binomial, x=TRUE, data=x_data.bstrp)
  
  Xmat.all.bstrp <- (propscore.model.bstrp$x[,-1])
  z.all.bstrp <- x_data.bstrp$AD
  
  # Rank based Mahalanobis distance
  distmat.all.bstrp <- rankmahal(z.all.bstrp, Xmat.all.bstrp)
  # Add caliper
  logit.propscore.all.bstrp <- (predict(propscore.model.bstrp))
  distmat2.all.bstrp <- addcaliper(distmat.all.bstrp, z.all.bstrp, logit.propscore.all.bstrp) # 186 x 223
  
  matchvec.full3.all.bstrp <- fullmatch(distmat2.all.bstrp, max.controls = 3, min.controls = 1/3)
  
  treated.subject.index.full3.all.bstrp <- rep(0,sum(z.all.bstrp==1))
  # The subject indices in the order of matchvec
  matchedset.index.full3.all.bstrp <- substr(matchvec.full3.all.bstrp, start=3, stop=10) # group name
  matchedset.index.full3.numeric.all.bstrp <- as.numeric(matchedset.index.full3.all.bstrp) # group name (numeric)
  subjects.match.order.full3.all.bstrp <- as.numeric(names(matchvec.full3.all.bstrp)) # subject indices
  
  # Create a numeric variable for which stratum each unit belongs to
  # 0 denotes that the unit was not matched
  # there is not units not matched b/c we use full matching
  stratum.short.full3.all.bstrp <- substr(matchvec.full3.all.bstrp, start=3, stop=10)
  stratum.full3.numeric.all.bstrp <- as.numeric(stratum.short.full3.all.bstrp) # i-th element : group name of unit i
  
  # Reassign numbers to each stratum that go from 1,..., no. of straum
  sort.unique.stratum.full3.all.bstrp <- sort(unique(stratum.full3.numeric.all.bstrp)) # 140 groups
  stratum.myindex.matchvecorder.full3.all.bstrp <- rep(0,length(stratum.full3.numeric.all.bstrp))
  for(i in 1:length(sort.unique.stratum.full3.all.bstrp)){
    stratum.myindex.matchvecorder.full3.all.bstrp[stratum.full3.numeric.all.bstrp==sort.unique.stratum.full3.all.bstrp[i]] <- i  # 140 matched group sense indices / i-th element : group name (1~140) of unit i
  }
  
  stratum.myindex.full3.all.bstrp <- rep(0, length(stratum.myindex.matchvecorder.full3.all.bstrp)) # length 409
  stratum.myindex.full3.all.bstrp[subjects.match.order.full3.all.bstrp] <- stratum.myindex.matchvecorder.full3.all.bstrp # i-th element : matched group index for unit i (1~140)
  
  stratumStructure(matchvec.full3.all.bstrp)
  # 3:1 2:1 1:1 1:2 1:3 
  # 16  14  63  11  36 
  
  
  # Calculate the standardized differences
  Xmat.matchorder.full3.all.bstrp <- Xmat.all[subjects.match.order.full3.all.bstrp,]
  
  std.diff.before.all.bstrp <- rep(0, ncol(Xmat.matchorder.full3.all.bstrp))
  std.diff.after.full3.all.bstrp <- rep(0, ncol(Xmat.matchorder.full3.all.bstrp))
  names(std.diff.before.all.bstrp) <- names(Xmat.matchorder.full3.all.bstrp[1,])
  names(std.diff.after.full3.all.bstrp) <- names(Xmat.matchorder.full3.all.bstrp[1,])
  
  for(i in 1:ncol(Xmat.matchorder.full3.all.bstrp)){
    temp.stand.diff.bstrp <- standardized.diff.func(Xmat.all.bstrp[,i], z.all.bstrp, stratum.myindex.full3.all.bstrp)
    std.diff.before.all.bstrp[i] <- temp.stand.diff.bstrp$std.diff.before.matching
    std.diff.after.full3.all.bstrp[i] <- temp.stand.diff.bstrp$std.diff.after.matching
  }
  
  fullmatch.std.diff3.all.bstrp <- cbind(std.diff.before.all.bstrp, std.diff.after.full3.all.bstrp)
  
  group.use.bstrp <- stratum.myindex.full3.all.bstrp
  
  estimator <- 'l2'
  # regression with only z
  T2.res.l2.bstrp <- geo_reg(manifold, x_data.bstrp$AD, y_data2[,bstrp.sample], 
                                 w = calculateWeight(z.all.bstrp, group.use.bstrp, sum(z.all.bstrp)/length(z.all.bstrp), 1-sum(z.all.bstrp)/length(z.all.bstrp)), 
                                 estimator, max_iter = 200)
  
  T2.l2.bstrp <- sqrt(abs(sum(Conj(T2.res.l2.bstrp$V)*T2.res.l2.bstrp$V)))
  T2.l2.bstrp.list <- c(T2.l2.bstrp.list, T2.l2.bstrp)
  
  estimator <- 'l1'
  T2.res.l1.bstrp <- geo_reg(manifold, x_data.bstrp$AD, y_data2[,bstrp.sample], 
                             w = calculateWeight(z.all.bstrp, group.use.bstrp, sum(z.all.bstrp)/length(z.all.bstrp), 1-sum(z.all.bstrp)/length(z.all.bstrp)), 
                             estimator, p_tol=1e-6, V_tol=1e-6)
  
  T2.l1.bstrp <- sqrt(abs(sum(Conj(T2.res.l1.bstrp$V)*T2.res.l1.bstrp$V)))
  T2.l1.bstrp.list <- c(T2.l1.bstrp.list, T2.l1.bstrp)
}


bstrp.pivot.ci.l2 <- as.vector(c(2*T2.l2 - quantile(T2.l2.bstrp.list, 0.975), 2*T2.l2 - quantile(T2.l2.bstrp.list, 0.025)))
bstrp.pivot.ci.l1 <- as.vector(c(2*T2.l1 - quantile(T2.l1.bstrp.list, 0.975), 2*T2.l1 - quantile(T2.l1.bstrp.list, 0.025)))

res.bstrp <- cbind(rbind(bstrp.pivot.ci.l2, bstrp.pivot.ci.l1), c(T2.l2, T2.l1))
colnames(res.bstrp) <- c("Lower", "Upper", "Estimate")
rownames(res.bstrp) <- c("L2", "L1")
res.bstrp
