# This script runs 100 bootstrapped regressions for each sample in the provided expression
# dataset, calculating the expected activity of the RBP stability programs based on the
# provided stability code. For each sample, the results are written in a separate file,
# with each row of the file corresponding to one covariate, and each column corresponding
# to one bootstrap iteration.


library(glmnet)

# change this line to indicate a different expression dataset
exp <- read.csv("CAGEKID.ccRCC.matching_pair.tumour_vs_normal_log10.txt",sep="\t")

exp[is.na(exp)] <- 0
numExp=ncol(exp)-1
exp$maxAbs_logFC <- apply( abs(exp[,2:ncol(exp)]), 1, max )
exp <- exp[exp$maxAbs_logFC>=1,]

# change this line to indicate a different stability code
code <- read.csv("StabilityCode.high_confidence.RBP_3UTR.txt",sep="\t")

code[is.na(code)] <- 0
numRBP=ncol(code)-1
code$sumBinding <- apply( abs(code[,2:ncol(code)]), 1, sum )
code <- code[code$sumBinding>=1,]

data <- merge( exp, code, by="Gene" )

covariates <- as.matrix(data[,(numExp+3):(numExp+numRBP+2)])

for(sampleIndex in 1:numExp)
{
  position <- sampleIndex+1
  sampleName <- colnames(data)[position]
  print( paste("Analyzing sample ",sampleName,sep="") )

  coefs <- data.frame()
  
  alpha=1
  
  for(iter in 1:100)
  {
    id <- sample(1:nrow(data),nrow(data),replace=TRUE)
    list <- 1:nrow(data)
    predicted <- data.frame()
    observed <- data.frame()
    
    trainingset <- data[id,]
    trainingcovariates <- covariates[id,]
    
    testset <- subset(data, !(list %in% id))
    testcovariates <- subset(covariates, !(list %in% id))
    
    model <- glmnet( trainingcovariates, trainingset[,position], family="gaussian",alpha=alpha )
    
    temp <- predict(model, testcovariates, s=0.0)
    
    this_coefs <- as.matrix(coef(model,s=0.0))
    this_coefs <- rbind( this_coefs, cor(temp,testset[,position]) )
    rownames(this_coefs)[dim(this_coefs)[1]] <- "Correlation"
    
    if( iter==1 )
      coefs <- this_coefs
    else
      coefs <- cbind( coefs, this_coefs )
    
    print(iter)
    
  }
	
  # change this line to indicate a different output path
  write.table(coefs,paste("Ridge.bootstrap.",sampleName,".coefs.txt",sep=""),sep="\t",col.names=F,quote=F)
}
