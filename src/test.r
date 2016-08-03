if(is.loaded("FSToolboxR", type = "External")){
  dyn.unload("libFSToolboxR.so")
}
dyn.load("libFSToolboxR.so")

FSToolboxR <- function(algorithm,maxSelectedFeatures,featureMatrix,initialPreselected=NULL,target,optionalParam1=0.0,optionalParam2=0.0) {
  allowedAlgorihtms<-c("MIFS","mRMR","CMIM","MIM","JMI","DISR","JMI","CIFE","ICAP","CondRed","BetaGamma","CMI","MIM")
  
  #Algorithm checks
  if(is.null(algorithm))
    stop("Algorithm cannot be null or empty")
  if(!is.character(algorithm)) 
    stop(paste("Algorithm class(",class(Algorithm),") is not a character"))
  if(nchar(algorithm)==0)
    stop("Algorithm cannot be empty")
  if(!is.element(algorithm,allowedAlgorihtms))
    stop(paste("Unknown algorithm, accepted values:",paste(allowedAlgorihtms,collapse=" ")))
  algorithmChar=as.character(algorithm)
  
  #Feature matrix
  if(is.null(featureMatrix))
    stop("Feature matrix cannot be null or empty")
  if(!is.matrix(featureMatrix)) 
    stop(paste("featureMatrix class(",class(featureMatrix),") is not a matrix"))
  fcols=as.integer(ncol(featureMatrix))
  frows=as.integer(nrow(featureMatrix))
  vectorizedfeatureMatrix=as.vector(featureMatrix,mode="double")
  
  #maxSelectedFeatures checks
  if(is.null(maxSelectedFeatures))
    stop("maxSelectedFeatures cannot be null")
  if(!is.numeric(maxSelectedFeatures))
    stop(paste("maxSelectedFeatures class(",class(maxSelectedFeatures),") has to be a number(numeric)"))
  if(maxSelectedFeatures<=0)
    stop("maxSelectedFeatures cannot be <=0")
  if(maxSelectedFeatures>fcols)
    stop(paste("maxSelectedFeatures(",maxSelectedFeatures,"), bigger than the number of available features -> featureMatrix cols(",fcols,")"))
  maxSelectedFeatures=as.integer(maxSelectedFeatures)
  
  #initialPreselected 
  if(!is.null(initialPreselected)&&!is.vector(initialPreselected)) 
    stop(paste("initialPreselected class(",class(initialPreselected),")is not a  vector"))
  if(length(unique(initialPreselected))!=length(initialPreselected))
    stop("initialPreselected contains repeated items")
  if(!is.null(initialPreselected)&&length(initialPreselected)>fcols)
    stop(paste("Cannot select (",length(initialPreselected),") features because the maximum is (",fcols,")"))
  if(is.null(initialPreselected)){
    initialPreselected=as.vector(mode="double",0)
    numberOfPreselected=0;
  }else{
    initialPreselected=as.vector(c(initialPreselected),mode="double")
    numberOfPreselected=as.integer(nrow(initialPreselected))
  }
  
  #target matrix-vector
  if(is.null(target)) stop("target matrix cannot be null")
  if((!is.matrix(target)&&!is.vector(target))) stop("target is not a matrix or vector")
  if(is.matrix(target)&&ncol(target)>1) stop("target has to be one column matrix or a vector")
  if(is.matrix(target)&&ncol(target)<=0) stop("target cannot be empty")
  if(is.vector(target)&&length(target)<=0) stop("target cannot be empty")
  if(is.matrix(target)&&!(nrow(target)==frows)) stop(paste("target rows(",nrow(target),") has to be same number of featureMatrix(",frows,")"))
  if(is.vector(target)&&!(length(target)==frows)) stop(paste("target rows(",length(target),") has to be same number of featureMatrix(",frows,")"))
  if(!is.vector(target))
    targetVect=as.double(as.vector(target))
  else
    targetVect=as.double(target)
  
  #optionalParam
  if(!is.double(optionalParam1))
    stop("optionalParam1 has to be a number(double)")
  if(!is.double(optionalParam2))
    stop("optionalParam2 has to be a number(double)")
  optionalParam1=as.double(optionalParam1)
  optionalParam2=as.double(optionalParam2)
  
  outputFeatures=(integer(maxSelectedFeatures)-1)
  
  #void FSToolboxR(char **algorithm,int *k, double *vectorizedfeatureMatrix,int *featureAsColumn,int *numberOfFeatures,int *numberOfSamples,double *initialFeatures,
  #  int *noInitialFeatures, double *vectorizedClassColumn, int *numberOfTargets, double *optionalParam1,double *optionalParam2,double *selectedFeatures,
  #  int *nselectedFeatures,int *error_code,char **error_msg)
  
  
  result <- .C("FSToolboxR",
               algorithm=algorithmChar,
               maxSelectedFeatures=as.integer(maxSelectedFeatures),
               vectorizedfeatureMatrix=vectorizedfeatureMatrix,
               featureAsColumn=as.integer(0),
               numberOfFeatures=fcols,
               numberOfSamples=frows,
               vectorizedPreselected=initialPreselected,
               numberOfPreselected=numberOfPreselected,
               vectorizedClassColumn=targetVect,
               numberOfTargets=frows,
               optionalParam1=optionalParam1,
               optionalParam2=optionalParam2,
               outputFeatures=as.integer(outputFeatures),#reserva memoria
               noutputFeatures=as.integer(0),#reserva memoria
               error_code=as.integer(0),
               error_msg=as.character(""))
  if(!is.null(result$error_code)&&result$error_code != 0) stop(result$error_msg)
  if(result$noutputFeatures<=0) stop("No selected features")
  return(result$outputFeatures[1:result$noutputFeatures])
}

load("test.RData")
FSToolboxR("CMIM",10,features,NULL,alarm,0.0,0.0)
