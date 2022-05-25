# Bivariate Twin Model Resources
#######################
####    WARNING    ####
#######################
# Functions in this script are merely intended to make the code body of AutomatedTwinModel.R more readable!
# Functions are NOT intended to be reused as they rely HEAVILY on global environment variables.
# Also, loading this will overwrite the default mxCompare(),
#   so there's yet another reason not to load this in any script other than AutomatedTwinModel.R

# From miFunctions2.R
# Note some of these have been changed slightly in order to make them return instead of print results
# Functions to assign labels
labLower  <- function(lab,nv) { paste(lab,rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="_") }
labSdiag  <- function(lab,nv) { paste(lab,rev(nv+1-sequence(1:(nv-1))),rep(1:(nv-1),(nv-1):1),sep="_") }
labFullSq <- function(lab,nv) { paste(lab,1:nv,rep(1:nv,each=nv),sep="_") }
labDiag   <- function(lab,nv) { paste(lab,1:nv,1:nv,sep="_") }
labSymm   <- function(lab,nv) { paste(lab,rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="_") }
labFull   <- function(lab,nr,nc) { paste(lab,1:nr,rep(1:nc,each=nr),sep="_") }
# Function "parameterSpecifations()" prints labels of a MxMatrix with
# square brackets surrounding free parameters; returns a matrix of strings
# -----------------------------------------------------------------------
parameterSpecifications <- function(model) {
  resultsList <- .collectParameterSpecifications(model)
  if(length(resultsList) > 0) {
    resultsNames <- names(resultsList)
    for(i in 1:length(resultsList)) {
      cat(resultsNames[[i]],'\n')
      print(resultsList[[i]], quote=FALSE)
      cat('\n')
    }
  }
}

.collectParameterSpecifications <- function(model) {
  listReturn <- list()
  if(length(model@matrices) > 0) {
    for(i in 1:length(model@matrices)) {
      current <- model@matrices[[i]]
      extract <- is(current, "FullMatrix") ||
        is(current, "LowerMatrix") ||
        is(current, "DiagMatrix") ||
        is(current, "SymmMatrix") ||
        is(current, "StandMatrix")
      if(extract) {
        retval <- mapply(.parameterSpecificationsHelper,
                         current@labels, current@free, current@values)
        retval <- matrix(retval, nrow(current), ncol(current))
        dimnames(retval) <- dimnames(current)
        storeName <- paste('model:', model@name,', matrix:', current@name, sep='')
        listReturn[[storeName]] <- retval
      }
    }
  }
}
.parameterSpecificationsHelper <- function(label, free, value) {
  if(free) return(paste('[', label, ']', sep = ''))
  else return(value)
}
# Function "formatOutputMatrices()" prints matrix with specified labels and
# number of decimals
# -----------------------------------------------------------------------
#parse(text=matricesList[k]) == matricesList[[k]]
formatOutputMatrices <- function(fittedModel,matricesList,labelsList,vars,digits) {
  result <- list()
  if(length(matricesList) > 0) {
    for(k in 1:length(matricesList)) {
      result[[matricesList[[k]]]] <- formatOutputMatrix(
        evalQuote(matricesList[[k]], fittedModel),
        labelsList[[k]],vars,digits)
    }
  }
  return(result)
}
formatOutputMatrix <- function(matrix,label,vars,digits) {
  #table <- round(eval(substitute(mxEval(Matrix,Model))),ND)
  # matrix <- apply(matrix, c(1,2), round, digits = digits)
  retval <- apply(matrix, c(1,2), format, scientific=FALSE, nsmall = 10)

  cols <- character(ncol(retval))
  for(i in 1:ncol(retval)) {paste(label,i,sep="")} -> cols[i]
  colnames(retval) <- cols
  if (nrow(retval) == length(vars)) {
    rownames(retval) <- vars
  } else {
    rows <- character(nrow(retval))
    for(j in 1:nrow(retval)) {paste("LP",j,sep="")} -> rows[j]
    rownames(retval) <- rows
  }
  return(retval)
}
# Function "formatMatrix()" returns a matrix with specified dimnames and # of decimal places
# -----------------------------------------------------------------------
formatMatrix <- function(matrix, dimnames, digits) {
  retval <- apply(matrix, c(1,2), round, digits)
  dimnames(retval) <- dimnames
  return(retval)
}

evalQuote <- function(expstring, model, compute = FALSE, show = FALSE) {
  return(eval(substitute(mxEval(x, model, compute, show),
                         list(x = parse(text=expstring)[[1]]))))
}
# Functions to generate output

fitGofs   <- function(fit) {
  summ <- summary(fit)
  cat(paste0("Mx:", fit$name,"  os=", summ$ob,"  ns=", summ$nu,"   ep=", summ$es,
            "   co=", sum(summ$cons),"  df=", summ$de, "  ll=", round(summ$Mi,4),
            "  cpu=", round(summ$cpu,4),"  opt=", summ$op,"  ver=", summ$mx,
            "  stc=", fit$output$status$code, "\n"))
}

fitGofS   <- function(fit) {
  summ <- summary(fit)
  cat(paste0("Mx:", fit$name,"  #statistics=", summ$ob,"  #records=", summ$nu,"   #parameters=", summ$es,
            "   #constraints=", sum(summ$cons),"  df=", summ$de, "  -2LL=", round(summ$Mi,4),
            "  cpu=", round(summ$cpu,4),"  optim=", summ$op,"  version=", summ$mx,
            "  code=", fit$output$status$code, "\n"))
}

# Custom functions
mxCompare <- function(base, comparison, ...) {
  x <- OpenMx::mxCompare(base, comparison, ...)
  for (y in which(is.na(x[2:nrow(x), "p"]))+1) {
    warning("Model comparison of ", x[y, "base"], " and ", x[y, "comparison"], " failed, returning 0s instead.")
    x[y, c(7,8,9)] <- 0
  }
  return(x)
}

Saturated_Model <- function() {
  # Saturated Model
  # Create Algebra for expected Mean Matrices
  meanMZf    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMeMZf, name="meanMZf" )
  meanDZf    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMeDZf, name="meanDZf" )
  meanMZm    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMeMZm, name="meanMZm" )
  meanDZm    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMeDZm, name="meanDZm" )
  if (!(samesex_only)) {meanDZo <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMeDZo, name="meanDZo")}


  # Create Algebra for expected Variance/Covariance Matrices
  #cholMZ    <- mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=svVas, lbound=lbVas, labels=labCvMZ, name="cholMZ" )
  #cholDZ    <- mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=svVas, lbound=lbVas, labels=labCvDZ, name="cholDZ" )
  #covMZ     <- mxAlgebra( expression=cholMZ %*% t(cholMZ), name="covMZ" )
  #covDZ     <- mxAlgebra( expression=cholDZ %*% t(cholDZ), name="covDZ" )

  cholMZf    <- mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=svVas, lbound=lbVas, labels=labCvMZf, name="cholMZf" )
  cholDZf    <- mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=svVas, lbound=lbVas, labels=labCvDZf, name="cholDZf" )
  cholMZm    <- mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=svVas, lbound=lbVas, labels=labCvMZm, name="cholMZm" )
  cholDZm    <- mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=svVas, lbound=lbVas, labels=labCvDZm, name="cholDZm" )
  if (!(samesex_only)) {cholDZo    <- mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=svVas, lbound=lbVas, labels=labCvDZo, name="cholDZo" )}


  covMZf     <- mxAlgebra( expression=cholMZf %*% t(cholMZf), name="covMZf" )
  covDZf     <- mxAlgebra( expression=cholDZf %*% t(cholDZf), name="covDZf" )
  covMZm     <- mxAlgebra( expression=cholMZm %*% t(cholMZm), name="covMZm" )
  covDZm     <- mxAlgebra( expression=cholDZm %*% t(cholDZm), name="covDZm" )
  if (!(samesex_only)) {covDZo     <- mxAlgebra( expression=cholDZo %*% t(cholDZo), name="covDZo" )}



  # Create Data Objects for Multiple Groups
  #dataMZ    <- mxData( observed=mzData, type="raw" )
  #dataDZ    <- mxData( observed=dzData, type="raw" )

  dataMZf    <- mxData( observed=mzfData, type="raw" )
  dataDZf    <- mxData( observed=dzfData, type="raw" )
  dataMZm    <- mxData( observed=mzmData, type="raw" )
  dataDZm    <- mxData( observed=dzmData, type="raw" )
  if (!(samesex_only)) {dataDZo    <- mxData( observed=dzoData, type="raw" )}



  # Create Expectation Objects for Multiple Groups
  #expMZ     <- mxExpectationNormal( covariance="covMZ", means="meanMZ", dimnames=selVars )
  #expDZ     <- mxExpectationNormal( covariance="covDZ", means="meanDZ", dimnames=selVars )

  expMZf     <- mxExpectationNormal( covariance="covMZf", means="meanMZf", dimnames=selVars )
  expDZf     <- mxExpectationNormal( covariance="covDZf", means="meanDZf", dimnames=selVars )
  expMZm     <- mxExpectationNormal( covariance="covMZm", means="meanMZm", dimnames=selVars )
  expDZm     <- mxExpectationNormal( covariance="covDZm", means="meanDZm", dimnames=selVars )
  if (!(samesex_only)) {expDZo     <- mxExpectationNormal( covariance="covDZo", means="meanDZo", dimnames=selVars )}

  funML     <- mxFitFunctionML()

  corMZf     <- mxAlgebra( cov2cor(covMZf), name="corMZf" )
  corDZf     <- mxAlgebra( cov2cor(covDZf), name="corDZf" )
  corMZm     <- mxAlgebra( cov2cor(covMZm), name="corMZm" )
  corDZm     <- mxAlgebra( cov2cor(covDZm), name="corDZm" )
  if (!(samesex_only)) {corDZo     <- mxAlgebra( cov2cor(covDZo), name="corDZo" )}


  # Create Model Objects for Multiple Groups
  #modelMZ   <- mxModel( "MZ", meanMZ, cholMZ, covMZ, dataMZ, expMZ, funML )
  #modelDZ   <- mxModel( "DZ", meanDZ, cholDZ, covDZ, dataDZ, expDZ, funML )

  modelMZf   <- mxModel( "MZf", meanMZf, cholMZf, covMZf, dataMZf, expMZf, funML, corMZf )
  modelDZf   <- mxModel( "DZf", meanDZf, cholDZf, covDZf, dataDZf, expDZf, funML, corDZf )
  modelMZm   <- mxModel( "MZm", meanMZm, cholMZm, covMZm, dataMZm, expMZm, funML, corMZm )
  modelDZm   <- mxModel( "DZm", meanDZm, cholDZm, covDZm, dataDZm, expDZm, funML, corDZm )
  if (!(samesex_only)) {
    modelDZo   <- mxModel( "DZo", meanDZo, cholDZo, covDZo, dataDZo, expDZo, funML, corDZo )
    multi     <- mxFitFunctionMultigroup( c("MZf","DZf","MZm","DZm", "DZo") )
  } else {
    multi     <- mxFitFunctionMultigroup( c("MZf","DZf","MZm","DZm") )
  }

  # Create Confidence Interval Objects
  #ciCov     <- mxCI( c('MZ.covMZ','DZ.covDZ') )
  #ciMean    <- mxCI( c('MZ.meanMZ','DZ.meanDZ') )

  # Build Saturated Model with Confidence Intervals
  if (!(samesex_only)) {
    return(mxModel( "mulSATc", modelMZf, modelDZf, modelMZm, modelDZm, modelDZo, multi))#, ciCov, ciMean )
  } else{
    return(mxModel( "mulSATc", modelMZf, modelDZf, modelMZm, modelDZm, multi))#, ciCov, ciMean )
  }

}

ADE_Model <- function() {
  # Create Algebra for expected Mean Matrices
  meanG     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe[1:length(labMe)], labels=labMe, name="meanG" )
  meanGf     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMef, name="meanGf" )
  meanGm     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMem, name="meanGm" )
  meanGo     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMeo, name="meanGo" )

  # Create Matrices for Path Coefficients
  pathAf     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPaD, label=labLower("af",nv), lbound=lbPaD, name="af" )
  pathDf     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPaD, label=labLower("df",nv), lbound=lbPaD, name="df" )
  pathEf     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPeD, label=labLower("ef",nv), lbound=lbPaD, name="ef" )

  pathAm     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPaD, label=labLower("am",nv), lbound=lbPaD, name="am" )
  pathDm     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPaD, label=labLower("dm",nv), lbound=lbPaD, name="dm" )
  pathEm     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPeD, label=labLower("em",nv), lbound=lbPaD, name="em" )

  #pathRg    <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=1, label=c("rg11"), lbound=0, ubound=1, name="rg" )

  # Create Algebra for Variance Components
  covAf      <- mxAlgebra( expression=af %*% t(af), name="Af" )
  covDf      <- mxAlgebra( expression=df %*% t(df), name="Df" )
  covEf      <- mxAlgebra( expression=ef %*% t(ef), name="Ef" )

  covAm      <- mxAlgebra( expression=am %*% t(am), name="Am" )
  covDm      <- mxAlgebra( expression=dm %*% t(dm), name="Dm" )
  covEm      <- mxAlgebra( expression=em %*% t(em), name="Em" )

  # Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
  covPf      <- mxAlgebra( expression= Af+Df+Ef, name="Vf" )
  covPm      <- mxAlgebra( expression= Am+Dm+Em, name="Vm" )
  if (include_sibs) {
    expCovMZf  <- mxAlgebra( expression= rbind(cbind(Af+Df+Ef,                               Af+Df,                                    0.5%x%(af%*%t(am))+0.25%x%(df%*%t(dm)),  0.5%x%Af+ 0.25%x%Df),
                                               cbind(Af+Df,                                  Af+Df+Ef,                                 0.5%x%(af%*%t(am))+0.25%x%(df%*%t(dm)),  0.5%x%Af+ 0.25%x%Df),
                                               cbind(0.5%x%(am%*%t(af))+0.25%x%(dm%*%t(df)), 0.5%x%(am %*% t(af))+0.25%x%(dm%*%t(df)), Am+Dm+Em,                                0.5%x%(am %*% t(af))+0.25%x%(dm%*%t(df))),
                                               cbind(0.5%x%Af+ 0.25%x%Df,                    0.5%x%Af+ 0.25%x%Df,                      0.5%x%(af%*%t(am))+0.25%x%(df%*%t(dm)),  Af+Df+Ef)), name="expCovMZf")

    expCovDZf  <- mxAlgebra( expression= rbind(cbind(Af+Df+Ef,                               0.5%x%Af+ 0.25%x%Df,                      0.5%x%(af%*%t(am))+0.25%x%(df%*%t(dm)),  0.5%x%Af+ 0.25%x%Df),
                                               cbind(0.5%x%Af+ 0.25%x%Df,                    Af+Df+Ef,                                 0.5%x%(af%*%t(am))+0.25%x%(df%*%t(dm)),  0.5%x%Af+ 0.25%x%Df),
                                               cbind(0.5%x%(am%*%t(af))+0.25%x%(dm%*%t(df)), 0.5%x%(am %*% t(af))+0.25%x%(dm%*%t(df)), Am+Dm+Em,                                0.5%x%(am%*%t(af))+0.25%x%(dm%*%t(df))),
                                               cbind(0.5%x%Af+ 0.25%x%Df,                    0.5%x%Af+ 0.25%x%Df,                      0.5%x%(af%*%t(am))+0.25%x%(df%*%t(dm)),  Af+Df+Ef)), name="expCovDZf")

    expCovMZm  <- mxAlgebra( expression= rbind(cbind(Am+Dm+Em,                               Am+Dm,                                  0.5%x%Am+ 0.25%x%Dm,                     0.5%x%(am%*%t(af))+0.25%x%(dm%*%t(df))),
                                               cbind(Am+Dm,                                  Am+Dm+Em,                               0.5%x%Am+ 0.25%x%Dm,                     0.5%x%(am%*%t(af))+0.25%x%(dm%*%t(df))),
                                               cbind(0.5%x%Am+ 0.25%x%Dm,                    0.5%x%Am+ 0.25%x%Dm,                    Am+Dm+Em,                                0.5%x%(am%*%t(af))+0.25%x%(dm%*%t(df))),
                                               cbind(0.5%x%(af%*%t(am))+0.25%x%(df%*%t(dm)), 0.5%x%(af%*%t(am))+0.25%x%(df%*%t(dm)), 0.5%x%(af%*%t(am))+0.25%x%(df%*%t(dm)),  Af+Df+Ef)), name="expCovMZm")

    expCovDZm  <- mxAlgebra( expression= rbind(cbind(Am+Dm+Em,                               0.5%x%Am+ 0.25%x%Dm,                    0.5%x%Am+ 0.25%x%Dm,                     0.5%x%(am%*%t(af))+0.25%x%(dm%*%t(df))),
                                               cbind(0.5%x%Am+ 0.25%x%Dm,                    Am+Dm+Em,                               0.5%x%Am+ 0.25%x%Dm,                     0.5%x%(am%*%t(af))+0.25%x%(dm%*%t(df))),
                                               cbind(0.5%x%Am+ 0.25%x%Dm,                    0.5%x%Am+ 0.25%x%Dm,                    Am+Dm+Em,                                0.5%x%(am%*%t(af))+0.25%x%(dm%*%t(df))),
                                               cbind(0.5%x%(af%*%t(am))+0.25%x%(df%*%t(dm)), 0.5%x%(af%*%t(am))+0.25%x%(df%*%t(dm)), 0.5%x%(af%*%t(am))+0.25%x%(df%*%t(dm)),  Af+Df+Ef)), name="expCovDZm")
    if (!(samesex_only)) {
        expCovDZo  <- mxAlgebra( expression= rbind(cbind(Am+Dm+Em,                                       0.5%x%(am%*%t(af))+ 0.25%x%(dm %*% t(df)), 0.5%x%Am+ 0.25%x%Dm,                     0.5%x%(am%*%t(af))+0.25%x%(dm%*%t(df))),
                                               cbind(0.5%x%(af%*%t(am))+ 0.25%x%(df %*% t(dm)),      Af+Df+Ef,                                  0.5%x%(af%*%t(am))+0.25%x%(df%*%t(dm)),  0.5%x%Af+ 0.25%x%Df),
                                               cbind(0.5%x%Am+ 0.25%x%Dm,                            0.5%x%Am+ 0.25%x%Dm,                       Am+Dm+Em,                                0.5%x%(am%*%t(af))+0.25%x%(dm%*%t(df))),
                                               cbind(0.5%x%(af%*%t(am))+0.25%x%(df%*%t(dm)),         0.5%x%Af+ 0.25%x%Df,                       0.5%x%(af%*%t(am))+0.25%x%(df%*%t(dm)),  Af+Df+Ef)), name="expCovDZo")
    }
  } else {
    expCovMZf  <- mxAlgebra( expression= rbind(cbind(Af+Df+Ef,                               Af+Df),
                                               cbind(Af+Df,                                  Af+Df+Ef)), name="expCovMZf")

    expCovDZf  <- mxAlgebra( expression= rbind(cbind(Af+Df+Ef,                               0.5%x%Af+ 0.25%x%Df),
                                               cbind(0.5%x%Af+ 0.25%x%Df,                    Af+Df+Ef)), name="expCovDZf")

    expCovMZm  <- mxAlgebra( expression= rbind(cbind(Am+Dm+Em,                               Am+Dm),
                                               cbind(Am+Dm,                                  Am+Dm+Em)), name="expCovMZm")

    expCovDZm  <- mxAlgebra( expression= rbind(cbind(Am+Dm+Em,                               0.5%x%Am+ 0.25%x%Dm),
                                               cbind(0.5%x%Am+ 0.25%x%Dm,                    Am+Dm+Em)), name="expCovDZm")
    if (!(samesex_only)) {
        expCovDZo  <- mxAlgebra( expression= rbind(cbind(Am+Dm+Em,                                       0.5%x%(am%*%t(af))+ 0.25%x%(dm %*% t(df))),
                                               cbind(0.5%x%(af%*%t(am))+ 0.25%x%(df %*% t(dm)),      Af+Df+Ef)), name="expCovDZo")
    }
  }



  # Create Algebra for Standardization
  matI      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I")
  invSDf     <- mxAlgebra( expression=solve(sqrt(I*Vf)), name="iSDf")
  invSDm     <- mxAlgebra( expression=solve(sqrt(I*Vm)), name="iSDm")

  # Calculate genetic and environmental correlations
  corAf      <- mxAlgebra( expression=solve(sqrt(I*Af))%&%Af, name ="rAf" ) #cov2cor()
  corDf      <- mxAlgebra( expression=solve(sqrt(I*Df))%&%Df, name ="rDf" )
  corEf      <- mxAlgebra( expression=solve(sqrt(I*Ef))%&%Ef, name ="rEf" )

  corAm      <- mxAlgebra( expression=solve(sqrt(I*Am))%&%Am, name ="rAm" ) #cov2cor()
  corDm      <- mxAlgebra( expression=solve(sqrt(I*Dm))%&%Dm, name ="rDm" )
  corEm      <- mxAlgebra( expression=solve(sqrt(I*Em))%&%Em, name ="rEm" )
  # Create Data Objects for Multiple Groups
  dataMZf    <- mxData( observed=mzfData, type="raw" )
  dataDZf    <- mxData( observed=dzfData, type="raw" )
  dataMZm    <- mxData( observed=mzmData, type="raw" )
  dataDZm    <- mxData( observed=dzmData, type="raw" )
  if (!(samesex_only)) {dataDZo    <- mxData( observed=dzoData, type="raw" )}

  # Create Expectation Objects for Multiple Groups
  expMZf     <- mxExpectationNormal( covariance="expCovMZf", means="meanGf", dimnames=selVars )
  expDZf     <- mxExpectationNormal( covariance="expCovDZf", means="meanGf", dimnames=selVars )
  expMZm     <- mxExpectationNormal( covariance="expCovMZm", means="meanGm", dimnames=selVars )
  expDZm     <- mxExpectationNormal( covariance="expCovDZm", means="meanGm", dimnames=selVars )
  if (!(samesex_only)) {expDZo     <- mxExpectationNormal( covariance="expCovDZo", means="meanGo", dimnames=selVars )}
  funML     <- mxFitFunctionML()

  # Create Model Objects for Multiple Groups
  parsf      <- list(meanGf, matI, invSDf,
                     pathAf, pathDf, pathEf, covAf, covDf, covEf, covPf, corAf, corDf, corEf)
  parsm      <- list(meanGm, matI, invSDm,
                     pathAm, pathDm, pathEm, covAm, covDm, covEm, covPm, corAm, corDm, corEm)
  #parso      <- list(pathRg)
  modelMZf   <- mxModel( name="MZf", parsf, parsm, meanGf,expCovMZf, dataMZf, expMZf, funML )
  modelDZf   <- mxModel( name="DZf", parsf, parsm, meanGf,expCovDZf, dataDZf, expDZf, funML )
  modelMZm   <- mxModel( name="MZm", parsm,  parsf,meanGm,expCovMZm, dataMZm, expMZm, funML )
  modelDZm   <- mxModel( name="DZm", parsm,  parsf, meanGm,expCovDZm, dataDZm, expDZm, funML )
  if (!(samesex_only)) {
    modelDZo   <- mxModel( name="DZo", parsf, parsm,  meanGo, expCovDZo, dataDZo, expDZo, funML )
    multi     <- mxFitFunctionMultigroup( c("MZf","DZf","MZm","DZm","DZo") )
  } else {
    multi     <- mxFitFunctionMultigroup( c("MZf","DZf","MZm","DZm") )
  }

  # Create Algebra for Variance Components
  colVC     <- rep('VC',nv)#vars
  rowVCf     <- rep(c('Af','Df','Ef','SAf','SDf','SEf'),each=nv)
  estVCf     <- mxAlgebra( expression=rbind(Af,Df,Ef,Af/Vf,Df/Vf,Ef/Vf), name="VCf", dimnames=list(rowVCf,colVC))

  rowVCm     <- rep(c('Am','Dm','Em','SAm','SDm','SEm'),each=nv)
  estVCm     <- mxAlgebra( expression=rbind(Am,Dm,Em,Am/Vm,Dm/Vm,Em/Vm), name="VCm", dimnames=list(rowVCm,colVC))

  # Build Model with Confidence Intervals
  if (!(samesex_only)) {
    return(mxModel( "mulADEc", parsf, parsm, modelMZf, modelDZf, modelMZm, modelDZm, modelDZo, multi, estVCf, estVCm ))
  } else {
    return(mxModel( "mulADEc", parsf, parsm, modelMZf, modelDZf, modelMZm, modelDZm, multi, estVCf, estVCm ))
  }

}

ACE_Model <- function() {
  # ACE Model
  # Create Algebra for expected Mean Matrices
  meanG     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMe, name="meanG" )
  meanGf     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMef, name="meanGf" )
  meanGm     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMem, name="meanGm" )
  meanGo     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMeo, name="meanGo" )

  # Create Matrices for Path Coefficients
  pathAf     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPaD, label=labLower("af",nv), lbound=lbPaD, name="af" )
  pathDf     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPaD, label=labLower("df",nv), lbound=lbPaD, name="df" )
  pathEf     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPeD, label=labLower("ef",nv), lbound=lbPaD, name="ef" )

  pathAm     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPaD, label=labLower("am",nv), lbound=lbPaD, name="am" )
  pathDm     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPaD, label=labLower("dm",nv), lbound=lbPaD, name="dm" )
  pathEm     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=svPeD, label=labLower("em",nv), lbound=lbPaD, name="em" )

  #pathRg    <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=1, label=c("rg11"), lbound=0, ubound=1, name="rg" )

  # Create Algebra for Variance Components
  covAf      <- mxAlgebra( expression=af %*% t(af), name="Af" )
  covDf      <- mxAlgebra( expression=df %*% t(df), name="Df" )
  covEf      <- mxAlgebra( expression=ef %*% t(ef), name="Ef" )

  covAm      <- mxAlgebra( expression=am %*% t(am), name="Am" )
  covDm      <- mxAlgebra( expression=dm %*% t(dm), name="Dm" )
  covEm      <- mxAlgebra( expression=em %*% t(em), name="Em" )

  # Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
  covPf      <- mxAlgebra( expression= Af+Df+Ef, name="Vf" )
  covPm      <- mxAlgebra( expression= Am+Dm+Em, name="Vm" )
  if (include_sibs) {
    expCovMZf  <- mxAlgebra( expression= rbind(cbind(Af+Df+Ef,                               Af+Df,                                    0.5%x%(af%*%t(am))+(df%*%t(dm)),  0.5%x%Af+ Df),
                                               cbind(Af+Df,                                  Af+Df+Ef,                                 0.5%x%(af%*%t(am))+(df%*%t(dm)),  0.5%x%Af+ Df),
                                               cbind(0.5%x%(am%*%t(af))+(dm%*%t(df)), 0.5%x%(am %*% t(af))+(dm%*%t(df)), Am+Dm+Em,                                0.5%x%(am %*% t(af))+(dm%*%t(df))),
                                               cbind(0.5%x%Af+ Df,                    0.5%x%Af+ Df,                      0.5%x%(af%*%t(am))+(df%*%t(dm)),  Af+Df+Ef)), name="expCovMZf")

    expCovDZf  <- mxAlgebra( expression= rbind(cbind(Af+Df+Ef,                               0.5%x%Af+ Df,                      0.5%x%(af%*%t(am))+(df%*%t(dm)),  0.5%x%Af+ Df),
                                               cbind(0.5%x%Af+ Df,                    Af+Df+Ef,                                 0.5%x%(af%*%t(am))+(df%*%t(dm)),  0.5%x%Af+ Df),
                                               cbind(0.5%x%(am%*%t(af))+(dm%*%t(df)), 0.5%x%(am %*% t(af))+(dm%*%t(df)), Am+Dm+Em,                                0.5%x%(am%*%t(af))+(dm%*%t(df))),
                                               cbind(0.5%x%Af+ Df,                    0.5%x%Af+ Df,                      0.5%x%(af%*%t(am))+(df%*%t(dm)),  Af+Df+Ef)), name="expCovDZf")

    expCovMZm  <- mxAlgebra( expression= rbind(cbind(Am+Dm+Em,                               Am+Dm,                                  0.5%x%Am+ Dm,                     0.5%x%(am%*%t(af))+(dm%*%t(df))),
                                               cbind(Am+Dm,                                  Am+Dm+Em,                               0.5%x%Am+ Dm,                     0.5%x%(am%*%t(af))+(dm%*%t(df))),
                                               cbind(0.5%x%Am+ Dm,                    0.5%x%Am+ Dm,                    Am+Dm+Em,                                0.5%x%(am%*%t(af))+(dm%*%t(df))),
                                               cbind(0.5%x%(af%*%t(am))+(df%*%t(dm)), 0.5%x%(af%*%t(am))+(df%*%t(dm)), 0.5%x%(af%*%t(am))+(df%*%t(dm)),  Af+Df+Ef)), name="expCovMZm")

    expCovDZm  <- mxAlgebra( expression= rbind(cbind(Am+Dm+Em,                               0.5%x%Am+ Dm,                    0.5%x%Am+ Dm,                     0.5%x%(am%*%t(af))+(dm%*%t(df))),
                                               cbind(0.5%x%Am+ Dm,                    Am+Dm+Em,                               0.5%x%Am+ Dm,                     0.5%x%(am%*%t(af))+(dm%*%t(df))),
                                               cbind(0.5%x%Am+ Dm,                    0.5%x%Am+ Dm,                    Am+Dm+Em,                                0.5%x%(am%*%t(af))+(dm%*%t(df))),
                                               cbind(0.5%x%(af%*%t(am))+(df%*%t(dm)), 0.5%x%(af%*%t(am))+(df%*%t(dm)), 0.5%x%(af%*%t(am))+(df%*%t(dm)),  Af+Df+Ef)), name="expCovDZm")
    if (!(samesex_only)) {
      expCovDZo  <- mxAlgebra( expression= rbind(cbind(Am+Dm+Em,                              0.5%x%(am%*%t(af))+ (dm %*% t(df)),  0.5%x%Am+ Dm,                     0.5%x%(am%*%t(af))+(dm%*%t(df))),
                                               cbind(0.5%x%(af%*%t(am))+ (df %*% t(dm)),      Af+Df+Ef,                            0.5%x%(af%*%t(am))+(df%*%t(dm)),  0.5%x%Af+ Df),
                                               cbind(0.5%x%Am+ Dm,                            0.5%x%Am+ Dm,                        Am+Dm+Em,                         0.5%x%(am%*%t(af))+(dm%*%t(df))),
                                               cbind(0.5%x%(af%*%t(am))+(df%*%t(dm)),         0.5%x%Af+ Df,                        0.5%x%(af%*%t(am))+(df%*%t(dm)),  Af+Df+Ef)), name="expCovDZo")
    }
  } else {
    expCovMZf  <- mxAlgebra( expression= rbind(cbind(Af+Df+Ef,                               Af+Df),
                                               cbind(Af+Df,                                  Af+Df+Ef)), name="expCovMZf")

    expCovDZf  <- mxAlgebra( expression= rbind(cbind(Af+Df+Ef,                               0.5%x%Af+ Df),
                                               cbind(0.5%x%Af+ Df,                    Af+Df+Ef)), name="expCovDZf")

    expCovMZm  <- mxAlgebra( expression= rbind(cbind(Am+Dm+Em,                               Am+Dm),
                                               cbind(Am+Dm,                                  Am+Dm+Em)), name="expCovMZm")

    expCovDZm  <- mxAlgebra( expression= rbind(cbind(Am+Dm+Em,                               0.5%x%Am+ Dm),
                                               cbind(0.5%x%Am+ Dm,                    Am+Dm+Em)), name="expCovDZm")
    if (!(samesex_only)) {
      expCovDZo  <- mxAlgebra( expression= rbind(cbind(Am+Dm+Em,                              0.5%x%(am%*%t(af))+ (dm %*% t(df))),
                                               cbind(0.5%x%(af%*%t(am))+ (df %*% t(dm)),      Af+Df+Ef)), name="expCovDZo")
    }
  }


  # Create Algebra for Standardization
  matI      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I")
  invSDf     <- mxAlgebra( expression=solve(sqrt(I*Vf)), name="iSDf")
  invSDm     <- mxAlgebra( expression=solve(sqrt(I*Vm)), name="iSDm")

  # Calculate genetic and environmental correlations
  corAf      <- mxAlgebra( expression=solve(sqrt(I*Af))%&%Af, name ="rAf" ) #cov2cor()
  corDf      <- mxAlgebra( expression=solve(sqrt(I*Df))%&%Df, name ="rDf" )
  corEf      <- mxAlgebra( expression=solve(sqrt(I*Ef))%&%Ef, name ="rEf" )

  corAm      <- mxAlgebra( expression=solve(sqrt(I*Am))%&%Am, name ="rAm" ) #cov2cor()
  corDm      <- mxAlgebra( expression=solve(sqrt(I*Dm))%&%Dm, name ="rDm" )
  corEm      <- mxAlgebra( expression=solve(sqrt(I*Em))%&%Em, name ="rEm" )
  # Create Data Objects for Multiple Groups
  dataMZf    <- mxData( observed=mzfData, type="raw" )
  dataDZf    <- mxData( observed=dzfData, type="raw" )
  dataMZm    <- mxData( observed=mzmData, type="raw" )
  dataDZm    <- mxData( observed=dzmData, type="raw" )
  if (!(samesex_only)) {dataDZo    <- mxData( observed=dzoData, type="raw" )}


  # Create Expectation Objects for Multiple Groups
  expMZf     <- mxExpectationNormal( covariance="expCovMZf", means="meanGf", dimnames=selVars )
  expDZf     <- mxExpectationNormal( covariance="expCovDZf", means="meanGf", dimnames=selVars )
  expMZm     <- mxExpectationNormal( covariance="expCovMZm", means="meanGm", dimnames=selVars )
  expDZm     <- mxExpectationNormal( covariance="expCovDZm", means="meanGm", dimnames=selVars )
  if (!(samesex_only)) {expDZo     <- mxExpectationNormal( covariance="expCovDZo", means="meanGo", dimnames=selVars )}
  funML     <- mxFitFunctionML()

  # Create Model Objects for Multiple Groups
  parsf      <- list(meanGf, matI, invSDf,
                     pathAf, pathDf, pathEf, covAf, covDf, covEf, covPf, corAf, corDf, corEf)
  parsm      <- list(meanGm, matI, invSDm,
                     pathAm, pathDm, pathEm, covAm, covDm, covEm, covPm, corAm, corDm, corEm)
  #parso      <- list(pathRg)
  modelMZf   <- mxModel( name="MZf", parsf, parsm, meanGf,expCovMZf, dataMZf, expMZf, funML )
  modelDZf   <- mxModel( name="DZf", parsf, parsm, meanGf,expCovDZf, dataDZf, expDZf, funML )
  modelMZm   <- mxModel( name="MZm", parsm,  parsf,meanGm,expCovMZm, dataMZm, expMZm, funML )
  modelDZm   <- mxModel( name="DZm", parsm,  parsf, meanGm,expCovDZm, dataDZm, expDZm, funML )
  if (!(samesex_only)) {
    modelDZo   <- mxModel( name="DZo", parsf, parsm,  meanGo, expCovDZo, dataDZo, expDZo, funML )
    multi     <- mxFitFunctionMultigroup( c("MZf","DZf","MZm","DZm","DZo") )
  } else {
    multi     <- mxFitFunctionMultigroup( c("MZf","DZf","MZm","DZm") )
  }

  # Create Algebra for Variance Components
  colVC     <- rep('VC',nv)#vars
  rowVCf     <- rep(c('Af','Df','Ef','SAf','SDf','SEf'),each=nv)
  estVCf     <- mxAlgebra( expression=rbind(Af,Df,Ef,Af/Vf,Df/Vf,Ef/Vf), name="VCf", dimnames=list(rowVCf,colVC))

  rowVCm     <- rep(c('Am','Dm','Em','SAm','SDm','SEm'),each=nv)
  estVCm     <- mxAlgebra( expression=rbind(Am,Dm,Em,Am/Vm,Dm/Vm,Em/Vm), name="VCm", dimnames=list(rowVCm,colVC))

  # Build Model with Confidence Intervals
  if (!(samesex_only)) {
   return(mxModel( "mulACEc", parsf, parsm, modelMZf, modelDZf, modelMZm, modelDZm, modelDZo, multi, estVCf, estVCm ))
  } else {
    return(mxModel( "mulACEc", parsf, parsm, modelMZf, modelDZf, modelMZm, modelDZm, multi, estVCf, estVCm ))
  }

}

Test_Correlation <- function(base_model, droplabel) {
  modelAE_rGE   <- mxModel( base_model, name=paste0("no", droplabel, collapse = ''))
  modelAE_rGE   <- omxSetParameters( modelAE_rGE, labels=droplabel, free=FALSE, values=0 )
  if (openmx_quiet) {
    quiet(fitAE_rGE     <- mxTryHard( modelAE_rGE, intervals=T, silent=TRUE, verbose=0 ), all=TRUE)
  } else {
    fitAE_rGE     <- mxTryHard( modelAE_rGE, intervals=T, silent=TRUE, verbose=0 )
  }
  return(mxCompare( base_model, fitAE_rGE))
}

Test_Multiple_Correlations <- function(base_model, base_model_name, labels_to_test, label_names, output_df, p_cutoff, rawresultsfilepath) {
  for (n in 1:length(labels_to_test)) {
    compare <- Test_Correlation(base_model, labels_to_test[n])
    if (save_extensive_raw_results) {
      write.table(data.frame(c('', paste0('########## ', base_model_name, ' vs dropped ', labels_to_test[n],'############'))),
                  rawresultsfilepath, sep='\t', append=TRUE, row.names = FALSE, col.names=FALSE, quote=F)
      suppressWarnings(write.table(compare, rawresultsfilepath, sep='\t', append=TRUE, quote=F))
    }
    compare_p <- compare$p[2]
    if (compare_p > p_cutoff) {
      output_df[, paste0(label_names[n], "_significant")] <- FALSE
    } else {
      output_df[, paste0(label_names[n], "_significant")] <- TRUE
    }
    output_df[, paste0(label_names[n], "_significant_p")] <- compare_p
  }
  return(output_df)
}

Test_Cor_Paths <- function(base_model, base_model_name, rA_paths, rE_paths, output_df, rawresultsfilepath, p_cutoff, sexlab) {
  comparison <- Test_Correlation(base_model, rA_paths)
  if (save_extensive_raw_results) {
    write.table(data.frame(c('', paste0('########## ', base_model_name,' vs dropped genetic pleiotropy paths', sexlab, ' ############'))), rawresultsfilepath, sep='\t', append=TRUE, row.names = FALSE, col.names=FALSE, quote=F)
    suppressWarnings(write.table(comparison, rawresultsfilepath, sep='\t', append=TRUE, quote=F))
  }
  compare_p <- comparison$p[2]
  if (compare_p > p_cutoff) {
    output_df[, paste0(sexlab, "_genetic_pleoitropy_paths_significant")] <- FALSE
  } else {
    output_df[, paste0(sexlab, "_genetic_pleoitropy_paths_significant")] <- TRUE
  }
  output_df[, paste0(sexlab, "_genetic_pleoitropy_paths_significant_p")] <- compare_p

  comparison <- Test_Correlation(base_model, rE_paths)
  if (save_extensive_raw_results) {
      write.table(data.frame(c('', paste0('########## ', base_model_name,' vs dropped rE paths ############'))), rawresultsfilepath, sep='\t', append=TRUE, row.names = FALSE, col.names=FALSE, quote=F)
  suppressWarnings(write.table(comparison, rawresultsfilepath, sep='\t', append=TRUE, quote=F))
  }
  compare_p <- comparison$p[2]
  if (compare_p > p_cutoff) {
    output_df[, paste0(sexlab, "_rE_paths_significant")] <- FALSE
  } else {
    output_df[, paste0(sexlab, "_rE_paths_significant")] <- TRUE
  }
  output_df[, paste0(sexlab, "_rE_paths_significant_p")] <- compare_p
  if ((output_df[, paste0(sexlab, "_rE_paths_significant")] ) & (output_df[, paste0(sexlab, "_genetic_pleoitropy_paths_significant")])) {
    output_df[, paste0(sexlab, "_all_correlations_significant")] <- TRUE
  } else {
    output_df[, paste0(sexlab, "_all_correlations_significant")] <- FALSE
  }
  return(output_df)
}