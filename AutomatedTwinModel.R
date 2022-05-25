rm(list=ls(all=TRUE)) # Clear working memory
##############################
###        README          ###
##############################
# Moved to https://github.com/matthijsz/AutomatedBivariateTwinModels

##############################
### GLOBAL SCRIPT SETTINGS ###
##############################
# For a description of these settings please see
# https://github.com/matthijsz/AutomatedBivariateTwinModels
xvars <- read.csv("BivarXnames.csv")[, "x"]
yvars <- read.csv("BivarYnames.csv")[, "x"]
wd <- "."
datapath <- "DataWide.csv"
output_file_name <- "Bivariate_model_results.csv"

include_sibs <- FALSE
preferred_model <- "ACE"
p_cutoff <- 0.01
save_extensive_raw_results <- FALSE
samesex_only <- FALSE
univar <- FALSE
force_AE <- FALSE
print_status <- TRUE

openmx_n_threads <- 2
openmx_optimizer <- "NPSOL"
openmx_quiet <- TRUE

##############################
###     DEPENDENCIES       ###
##############################

setwd(wd)
library(OpenMx)
library(psych)
library(ddpcr)
source("BivariateTwinModelsResources.R")

##############################
###   MAIN ANALYSIS CODE   ###
##############################

main <- function(twinData, xname, yname, modelno) {
  mxOption(NULL, key="Number of Threads", value=openmx_n_threads)
  mxOption( NULL, "Default optimizer", openmx_optimizer)
  # Load Data
  if (univar) {nv        <<- 1} else {nv        <<- 2}       # number of variables
  if (include_sibs) {nind <<- 4} else {nind <<- 2}
  ntv       <<- nv*nind    # number of total variables

  if (!(univar)) {
    rawresultsfilepath <<- paste0("RAW_OUTPUT_Bivariate_",yname, "_", xname,".txt")
    selVars <<- c(paste0(yname,"twin1"), paste0(xname,"twin1"),
                  paste0(yname,"twin2"), paste0(xname,"twin2"))
    labMe     <<- paste("mean",c(paste0(yname), paste0(xname)),sep="_")
    labMef     <<- c(paste0("meanf_", yname), paste0("meanf_", xname), paste0("meanf_", yname), paste0("meanf_", xname))
    labMem     <<- c(paste0("meanm_", yname), paste0("meanm_", xname), paste0("meanm_", yname), paste0("meanm_", xname))
    labMeo     <<- c(paste0("meanm_", yname), paste0("meanm_", xname), paste0("meanf_", yname), paste0("meanf_", xname))
    if (include_sibs) {
      selVars <<- c(selVars, paste0(yname,"sibm1"), paste0(xname,"sibm1"), paste0(yname,"sibf2"), paste0(xname,"sibf2"))
      labMef <<- c(labMef, paste0("meanm_", yname), paste0("meanm_", xname), paste0("meanf_", yname), paste0("meanf_", xname))
      labMem <<- c(labMem, paste0("meanm_", yname), paste0("meanm_", xname), paste0("meanf_", yname), paste0("meanf_", xname))
      labMeo <<- c(labMeo, paste0("meanm_", yname), paste0("meanm_", xname), paste0("meanf_", yname), paste0("meanf_", xname))
    }
  } else {
    rawresultsfilepath <<- paste0("RAW_OUTPUT_Univariate_", xname,".txt")
    selVars <<- c(paste0(xname,"twin1"), paste0(xname,"twin2"))
    labMe     <<- paste("mean",paste0(xname),sep="_")
    labMef     <<- c(paste0("meanf_", xname), paste0("meanf_", xname))
    labMem     <<- c(paste0("meanm_", xname), paste0("meanm_", xname))
    labMeo     <<- c(paste0("meanm_", xname), paste0("meanf_", xname))
    if (include_sibs) {
      selVars <<- c(selVars, paste0(xname,"sibm1"), paste0(xname,"sibf2"))
      labMef <<- c(labMef, paste0("meanm_", xname), paste0("meanf_", xname))
      labMem <<- c(labMem, paste0("meanm_", xname), paste0("meanf_", xname))
      labMeo <<- c(labMeo, paste0("meanm_", xname), paste0("meanf_", xname))
    }
  }


  if (save_extensive_raw_results) {
    if (!(univar)) {
      write.table(paste0("# ", yname," ~ ", xname), rawresultsfilepath, sep='\t', row.names = FALSE, col.names=FALSE, quote=F)
    } else {
      write.table(paste0("# ", xname), rawresultsfilepath, sep='\t', row.names = FALSE, col.names=FALSE, quote=F)
    }
  }

  # Select Data for Analysis
  assign("mzfData", subset(twinData, zyg==3, selVars), envir = .GlobalEnv)
  assign("dzfData", subset(twinData, zyg==4, selVars), envir = .GlobalEnv)
  assign("mzmData", subset(twinData, zyg==1, selVars), envir = .GlobalEnv)
  assign("dzmData", subset(twinData, zyg==2, selVars), envir = .GlobalEnv)
  if (!(samesex_only)) {
    assign("dzoData", subset(twinData, zyg==5, selVars), envir = .GlobalEnv)
  }

  # Guesstimate starting values
  if (!(samesex_only)) {
    i <<- rbind(round(colMeans(mzfData,na.rm=TRUE),4), round(colMeans(dzfData,na.rm=TRUE),4),
                round(colMeans(mzmData,na.rm=TRUE),4), round(colMeans(dzmData,na.rm=TRUE),4),
                round(colMeans(dzoData,na.rm=TRUE),4))
    j <<- c(diag(round(cov(mzfData,use="pairwise"),4)), diag(round(cov(dzfData,use="pairwise"),4)),
            diag(round(cov(mzmData,use="pairwise"),4)), diag(round(cov(dzmData,use="pairwise"),4)),
            diag(round(cov(dzoData,use="pairwise"),4)))
  } else {
    i <<- rbind(round(colMeans(mzfData,na.rm=TRUE),4), round(colMeans(dzfData,na.rm=TRUE),4),
                round(colMeans(mzmData,na.rm=TRUE),4), round(colMeans(dzmData,na.rm=TRUE),4))
    j <<- c(diag(round(cov(mzfData,use="pairwise"),4)), diag(round(cov(dzfData,use="pairwise"),4)),
            diag(round(cov(mzmData,use="pairwise"),4)), diag(round(cov(dzmData,use="pairwise"),4)))
  }
  svMe <<- unname(round(colMeans(i, na.rm=TRUE), 4))
  svVa <<- unname(round(mean(j, na.rm=TRUE), 4))

  # start value for variances
  svVas     <<- diag(svVa,ntv,ntv)    # assign start values to diagonal of matrix
  lbVa      <<- .0001                 # start value for lower bounds
  lbVas     <<- diag(lbVa,ntv,ntv)    # assign lower bounds values to diagonal of matrix
  lbVas[lower.tri(lbVas)] <<- -10     # lower bounds for below diagonal elements
  lbVas[upper.tri(lbVas)] <<- NA      # lower bounds for above diagonal elements

  # Create Labels
  labMeMZ   <<- paste("meanMZ",selVars,sep="_")
  labMeDZ   <<- paste("meanDZ",selVars,sep="_")
  labMeMZf   <<- paste("meanMZf",selVars,sep="_")
  labMeDZf   <<- paste("meanDZf",selVars,sep="_")
  labMeMZm   <<- paste("meanMZm",selVars,sep="_")
  labMeDZm   <<- paste("meanDZm",selVars,sep="_")
  labCvMZ   <<- labLower("covMZ",ntv)
  labCvMZf   <<- labLower("covMZf",ntv)
  labCvDZf   <<- labLower("covDZf",ntv)
  labCvMZm   <<- labLower("covMZm",ntv)
  labCvDZm   <<- labLower("covDZm",ntv)
  labCvDZ   <<- labLower("covDZ",ntv)
  labCvZ    <<- labLower("covZ",ntv)
  labVaMZf   <<- labDiag("covMZf",ntv)
  labVaDZf   <<- labDiag("covDZf",ntv)
  labVaMZm   <<- labDiag("covMZm",ntv)
  labVaDZm   <<- labDiag("covDZm",ntv)
  labVaDZo   <<- labDiag("covDZo",ntv)
  labVaZ    <<- labDiag("covZ",ntv)
  labVaMZ    <<- labDiag("covMZ",ntv)
  labVaDZ    <<- labDiag("covDZ",ntv)
  if (!(samesex_only)) {
    labCvDZo   <<- labLower("covDZo",ntv)
    labMeDZo   <<- paste("meanDZo",selVars,sep="_")
  }
  # ------------------------------------------------------------------------------
  # PREPARE MODEL
  if (print_status) {
    cat(paste0("Bivar model ", modelno," ",xname, " ~ ", yname, ":Fitting saturated model (1/8).\n"))
  }
  model <<- Saturated_Model()

  # ------------------------------------------------------------------------------
  # Run Saturated Model
  if (openmx_quiet) {
    quiet(fit       <<- mxTryHard( model, intervals=F, silent=TRUE, verbose=0 ), all=TRUE)
  } else {
    fit       <<- mxTryHard( model, intervals=F, silent=TRUE, verbose=0 )
  }

  assign("sum", summary( fit ), , envir = .GlobalEnv)

  ##----------------------------------------------------------
  #ADE model
  # Set Starting Values
  svPa      <<- .4                        # start value for path coefficient
  svPaD     <<- vech(diag(svPa,nv,nv))    # start values for diagonal of covariance matrix
  svPe      <<- .8                        # start value for path coefficient for e
  svPeD     <<- vech(diag(svPe,nv,nv))    # start values for diagonal of covariance matrix
  lbPa      <<- .0001                     # start value for lower bounds
  lbPaD     <<- diag(lbPa,nv,nv)          # lower bounds for diagonal of covariance matrix
  lbPaD[lower.tri(lbPaD)] <<- NA         # lower bounds for below diagonal elements
  lbPaD[upper.tri(lbPaD)] <<- NA          # lower bounds for above diagonal elements
  # ADE Model
  if (print_status) {
    if (univar) {cat(paste0("Univar model ",xname, ":Fitting ADE model (2/8).\n"))} else {cat(paste0("Bivar model ", modelno," ",xname, " ~ ", yname, ":Fitting ADE model (2/8).\n"))}
  }
  modelADE  <<- omxAssignFirstParameters(ADE_Model())
  # Run ADE Model
  if (openmx_quiet) {
    quiet(fitADE    <<- mxTryHard( modelADE, intervals=F, silent=TRUE, verbose=0), all=TRUE)
  } else {
    fitADE    <<- mxTryHard( modelADE, intervals=F, silent=TRUE, verbose=0)
  }
  sumADE    <<- summary( fitADE )
  parameterSpecifications(fitADE)
  # Compare with Saturated Model
  if (save_extensive_raw_results) {
    write.table(data.frame(c('', '########## SATURATED VS ADE ############')), rawresultsfilepath, sep='\t', append=TRUE, row.names = FALSE, col.names=FALSE, quote=F)
    suppressWarnings(write.table(mxCompare( fit, fitADE ), rawresultsfilepath, sep='\t', append=TRUE, quote=F))
  }

  # ------------------------------------------------------------------------------
  # ACE Model
  if (print_status) {
    cat(paste0("Bivar model ", modelno," ",xname, " ~ ", yname, ":Fitting ACE model (3/8).\n"))
  }
  modelACE  <<- omxAssignFirstParameters(ACE_Model())
  # Run ACE Model
  if (openmx_quiet) {
    quiet(fitACE    <<- mxTryHard( modelACE, intervals=F, silent=TRUE, verbose=0), all=TRUE)
  } else {
    fitADE    <<- mxTryHard( modelADE, intervals=F, silent=TRUE, verbose=0)
  }
  sumACE    <<- summary( fitACE )
  parameterSpecifications(fitACE)
  # Compare with Saturated Model
  if (save_extensive_raw_results) {
    write.table(data.frame(c('', '########## SATURATED VS ACE ############')), rawresultsfilepath, sep='\t', append=TRUE, row.names = FALSE, col.names=FALSE, quote=F)
    suppressWarnings(write.table(mxCompare( fit, fitACE), rawresultsfilepath, sep='\t', append=TRUE, quote=F))
  }

  # Select ADE or ACE for further analysis first based on p value, or else based on preference
  if (mxCompare( fit, fitADE )$p[2] > p_cutoff) {
    base_model_sig <<- TRUE
    if ((mxCompare( fit, fitACE )$p[2] > p_cutoff) & (preferred_model == "ACE")) {
      base_model_sig_p <<- mxCompare(fit, fitACE )$p[2]
      base_model_used <<- "ACE"
      fitADE <<- fitACE
    } else {
      base_model_sig_p <<- mxCompare(fit, fitADE )$p[2]
      base_model_used <<- "ADE"
      fitADE <<- fitADE
    }
  } else if (mxCompare( fit, fitACE )$p[2] > p_cutoff) {
    base_model_sig <<- TRUE
    base_model_sig_p <<- mxCompare(fit, fitACE )$p[2]
    base_model_used <<- "ACE"
    fitADE <<- fitACE
  } else if (preferred_model == "ADE") {
    base_model_sig <<- FALSE
    base_model_sig_p <<- mxCompare(fit, fitADE )$p[2]
    base_model_used <<- "ADE"
    fitADE <<- fitADE
  } else if (preferred_model == "ACE") {
    base_model_sig <<- FALSE
    base_model_sig_p <<- mxCompare(fit, fitACE )$p[2]
    base_model_used <<- "ACE"
    fitADE <<- fitACE
  }

  # ------------------------------------------------------------------------------
  # RUN SUBMODELS
  modelADEs   <<- mxModel( fitADE, name="mulADEc" )
  modelADEs   <<- omxSetParameters( modelADEs, labels=c("af_1_1", "am_1_1"), free=TRUE, values=.5, newlabels="a_1_1" )
  modelADEs   <<- omxSetParameters( modelADEs, labels=c("df_1_1", "dm_1_1"), free=TRUE, values=.5, newlabels="d_1_1" )
  modelADEs   <<- omxSetParameters( modelADEs, labels=c("ef_1_1", "em_1_1"), free=TRUE, values=.5, newlabels="e_1_1" )
  if (!(univar)) {
    modelADEs   <<- omxSetParameters( modelADEs, labels=c("af_2_1", "am_2_1"), free=TRUE, values=.5, newlabels="a_2_1" )
    modelADEs   <<- omxSetParameters( modelADEs, labels=c("af_2_2", "am_2_2"), free=TRUE, values=.5, newlabels="a_2_2" )

    modelADEs   <<- omxSetParameters( modelADEs, labels=c("df_2_1", "dm_2_1"), free=TRUE, values=.5, newlabels="d_2_1" )
    modelADEs   <<- omxSetParameters( modelADEs, labels=c("df_2_2", "dm_2_2"), free=TRUE, values=.5, newlabels="d_2_2" )

    modelADEs   <<- omxSetParameters( modelADEs, labels=c("ef_2_1", "em_2_1"), free=TRUE, values=.5, newlabels="e_2_1" )
    modelADEs   <<- omxSetParameters( modelADEs, labels=c("ef_2_2", "em_2_2"), free=TRUE, values=.5, newlabels="e_2_2" )
  }

  if (print_status) {
    if (!(univar)) {cat(paste0("Bivar model ", modelno," ",xname, " ~ ", yname, ":Fitting equal sexes model (4/8).\n"))} else {cat(paste0("Bivar model ",xname, " ~ ", yname, ":Fitting equal sexes model (4/8).\n"))}
  }
  # Run eqs model
  if (openmx_quiet) {
    quiet(fitADEs     <<- mxTryHard( modelADEs, intervals=F, silent=TRUE, verbose=0 ), all=TRUE)
  } else {
    fitADEs     <<- mxTryHard( modelADEs, intervals=F, silent=TRUE, verbose=0 )
  }
  if (save_extensive_raw_results) {
    write.table(data.frame(c('', paste0('########## ',base_model_used,' vs equal path-values for sexes ############'))), rawresultsfilepath, sep='\t', append=TRUE, row.names = FALSE, col.names=FALSE, quote=F)
    suppressWarnings(write.table(mxCompare( fitADE, fitADEs), rawresultsfilepath, sep='\t', append=TRUE, quote=F))
  }
  equal_sexes_ade_p <<- mxCompare( fitADE, fitADEs)$p[2]
  if (equal_sexes_ade_p > p_cutoff) {
    equal_sexes_ade <<- TRUE
    fitADE <<- fitADEs
  } else {
    equal_sexes_ade <<- FALSE
  }
  # Run AE model
  modelAE   <<- mxModel( fitADE, name="mulAEc" )
  if (equal_sexes_ade) {
    modelAE   <<- omxSetParameters( modelAE, labels=labLower("d",nv), free=FALSE, values=0 )
  } else {
    modelAE   <<- omxSetParameters( modelAE, labels=labLower("df",nv), free=FALSE, values=0 )
    modelAE   <<- omxSetParameters( modelAE, labels=labLower("dm",nv), free=FALSE, values=0 )
  }
  if (print_status) {
    cat(paste0("Bivar model ", modelno," ",xname, " ~ ", yname, ":Fitting AE model(5/8).\n"))
  }
  if (openmx_quiet) {
    quiet(fitAE     <<- mxTryHard( modelAE, intervals=F, silent=TRUE, verbose=0 ), all=TRUE)
  } else {
    fitAE     <<- mxTryHard( modelAE, intervals=F, silent=TRUE, verbose=0 )
  }
  if (save_extensive_raw_results) {
    write.table(data.frame(c('', paste0('########## ',base_model_used,' vs AE ############'))), rawresultsfilepath, sep='\t', append=TRUE, row.names = FALSE, col.names=FALSE, quote=F)
    suppressWarnings(write.table(mxCompare( fitADE, fitAE ), rawresultsfilepath, sep='\t', append=TRUE, quote=F))
  }
  # Run DE model
  modelDE   <<- mxModel( fitADE, name="mulDEc" )
  if (equal_sexes_ade) {
    modelDE   <<- omxSetParameters( modelDE, labels=labLower("a",nv), free=FALSE, values=0 )
  } else {
    modelDE   <<- omxSetParameters( modelDE, labels=labLower("af",nv), free=FALSE, values=0 )
    modelDE   <<- omxSetParameters( modelDE, labels=labLower("am",nv), free=FALSE, values=0 )
  }
  if (openmx_quiet) {
    quiet(fitDE     <<- mxTryHard( modelDE, intervals=F, silent=TRUE, verbose=0 ), all=TRUE)
  } else {
    fitDE     <<- mxTryHard( modelDE, intervals=F, silent=TRUE, verbose=0 )
  }
  if (save_extensive_raw_results) {
    if (base_model_used == "ADE") {
      write.table(data.frame(c('', paste0('########## ADE vs DE ############'))), rawresultsfilepath, sep='\t', append=TRUE, row.names = FALSE, col.names=FALSE, quote=F)
    } else if (base_model_used == "ACE") {
      write.table(data.frame(c('', paste0('########## ACE vs CE ############'))), rawresultsfilepath, sep='\t', append=TRUE, row.names = FALSE, col.names=FALSE, quote=F)
    }
    suppressWarnings(write.table(mxCompare( fitADE, fitDE ), rawresultsfilepath, sep='\t', append=TRUE, quote=F))
  }

  if (print_status) {
    cat(paste0("Bivar model ", modelno," ",xname, " ~ ", yname, ":Fiting DE model(6/8).\n"))
  }
  could_drop_a_p <<- mxCompare(fitADE, fitDE)$p[2]
  if (could_drop_a_p > p_cutoff) {
    could_drop_a <<- TRUE
  } else {
    could_drop_a <<- FALSE
  }
  drop_d_p <<- mxCompare( fitADE, fitAE )$p[2]
  if (univar) {
    if ((drop_d_p > p_cutoff) | (force_AE)) {
      # If D (or C) can be dropped, use AE as base model
      drop_d <<- TRUE
      # Get important estimates
      if (equal_sexes_ade) {
        female_ests <<- formatOutputMatrices(fitAE, c("Af","Ef","Vf","Af/Vf","Ef/Vf"), c("covAf","covEf","Var","stCovA","stCovE"), selVars,4)
        important_values <- data.frame(matrix(as.numeric(c(female_ests$`Af/Vf`[1,1], female_ests$`Ef/Vf`[1,1])), nrow=1))
        colnames(important_values) <<- c("Af_v1", "Ef_v1")
        important_values[, c("Am_v1", "Em_v1","Df_v1", "Dm_v1", "Cf_v1", "Cm_v1")] <<- NA
      } else {
        female_ests <<- formatOutputMatrices(fitAE, c("Af","Ef","Vf","Af/Vf","Ef/Vf"), c("covAf","covEf","Var","stCovA","stCovE"), selVars,4)
        male_ests <<- formatOutputMatrices(fitAE, c("Am","Em","Vm","Am/Vm","Em/Vm"), c("covA","covE","Var","stCovA","stCovE"), selVars,4)
        important_values <- data.frame(matrix(as.numeric(c(female_ests$`Af/Vf`[1,1], female_ests$`Ef/Vf`[1,1],
                                                           male_ests$`Am/Vm`[1,1], male_ests$`Em/Vm`[1,1])), nrow=1))
        colnames(important_values) <<- c("Af_v1", "Ef_v1",
                                         "Am_v1", "Em_v1")
        important_values[, c("Df_v1", "Dm_v1", "Cf_v1", "Cm_v1")] <<- NA
      }
      # Test significance of correlations in the model
    } else {
      drop_d <<- FALSE
      if (equal_sexes_ade) {
        female_ests <<- formatOutputMatrices(fitADE, c("Af","Df","Ef","Vf","Af/Vf","Df/Vf","Ef/Vf"), c("covAf","covDf","covEf","Var","stCovA","stCovD","stCovE"), selVars,4)
        important_values <- data.frame(matrix(as.numeric(c(female_ests$`Af/Vf`[1,1], female_ests$`Df/Vf`[1,1], female_ests$`Ef/Vf`[1,1])), nrow=1))
        if (base_model_used == "ADE") {
          colnames(important_values) <<- c("Af_v1","Df_v1", "Ef_v1")
          important_values[, c("Cf_v1", "Cm_v1")] <<- NA
        } else {
          colnames(important_values) <<- c("Af_v1","Cf_v1", "Ef_v1")
          important_values[, c("Df_v1", "Dm_v1")] <<- NA
        }
      } else {
        female_ests <<- formatOutputMatrices(fitADE, c("Af","Df","Ef","Vf","Af/Vf","Df/Vf","Ef/Vf"), c("covAf","covDf","covEf","Var","stCovA","stCovD","stCovE"), selVars,4)
        male_ests <<- formatOutputMatrices(fitADE, c("Am","Dm","Em","Vm","Am/Vm","Dm/Vm","Em/Vm"), c("covA","covD","covE","Var","stCovA","stCovD","stCovE"), selVars,4)
        important_values <- data.frame(matrix(as.numeric(c(female_ests$`Af/Vf`[1,1], female_ests$`Df/Vf`[1,1], female_ests$`Ef/Vf`[1,1],
                                                           male_ests$`Am/Vm`[1,1], male_ests$`Dm/Vm`[1,1], male_ests$`Em/Vm`[1,1])), nrow=1))
        if (base_model_used == "ADE") {
          colnames(important_values) <<- c("Af_v1","Df_v1", "Ef_v1", "Am_v1", "Dm_v1", "Em_v1")
          important_values[, c("Cf_v1", "Cm_v1")] <<- NA
        } else {
          colnames(important_values) <<- c("Af_v1","Cf_v1", "Ef_v1", "Am_v1", "Cm_v1", "Em_v1")
          important_values[, c("Df_v1", "Dm_v1")] <<- NA
        }
      }
    }
  } else {
    if ((drop_d_p > p_cutoff) | (force_AE)) {
      # If D (or C) can be dropped, use AE as base model
      drop_d <<- TRUE
      # Get important estimates
      if (equal_sexes_ade) {
        female_ests <<- formatOutputMatrices(fitAE, c("Af","Ef","Vf","Af/Vf","Ef/Vf"), c("covAf","covEf","Var","stCovA","stCovE"), selVars,4)
        female_cors <<- formatOutputMatrices(fitAE, c("solve(sqrt(I*Af)) %&% Af","solve(sqrt(I*Ef)) %&% Ef"), c("corA","corE"), selVars, 4)
        important_values <- data.frame(matrix(as.numeric(c(
          female_ests$`Af/Vf`[1,1], female_ests$`Af/Vf`[2,2], female_ests$`Af/Vf`[1,2],
          female_ests$`Ef/Vf`[1,1], female_ests$`Ef/Vf`[2,2], female_ests$`Ef/Vf`[1,2],
          female_cors$`solve(sqrt(I*Af)) %&% Af`[1, 2], female_cors$`solve(sqrt(I*Ef)) %&% Ef`[1, 2]
        )), nrow=1))
        cor_labels_in_model <<- c("a_2_1", "e_2_1")
        cor_labels_in_result <<- c("Af_bivar", "Ef_bivar")
        colnames(important_values) <<- c("Af_v1", "Af_v2", "Af_bivar", "Ef_v1", "Ef_v2", "Ef_bivar",
                                         "rAf", "rEf")
        important_values[, c("Am_v1", "Am_v2", "Am_bivar", "Em_v1", "Em_v2", "Em_bivar","rAm", "rEm","Df_v1", "Df_v2", "Df_bivar", "Dm_v1", "Dm_v2", "Dm_bivar", "rDm", "rDf", "Cf_v1", "Cf_v2", "Cf_bivar", "Cm_v1", "Cm_v2", "Cm_bivar", "rCm", "rCf")] <<- NA
        if (print_status) {
          cat(paste0("Bivar model ", modelno," ",xname, " ~ ", yname, ":Testing genetic and environmental cors (7/8).\n"))
        }
        important_values <- Test_Cor_Paths(base_model=fitAE, base_model_name='AE',
                                           rA_paths = 'a_2_1', rE_paths= 'e_2_1',
                                           output_df=important_values,
                                           rawresultsfilepath=rawresultsfilepath, p_cutoff=p_cutoff, sexlab='females')
      } else {
        female_ests <<- formatOutputMatrices(fitAE, c("Af","Ef","Vf","Af/Vf","Ef/Vf"), c("covAf","covEf","Var","stCovA","stCovE"), selVars,4)
        female_cors <<- formatOutputMatrices(fitAE, c("solve(sqrt(I*Af)) %&% Af","solve(sqrt(I*Ef)) %&% Ef"), c("corA","corE"), selVars, 4)
        male_ests <<- formatOutputMatrices(fitAE, c("Am","Em","Vm","Am/Vm","Em/Vm"), c("covA","covE","Var","stCovA","stCovE"), selVars,4)
        male_cors <<- formatOutputMatrices(fitAE, c("solve(sqrt(I*Am)) %&% Am","solve(sqrt(I*Em)) %&% Em"),  c("corA","corE"), selVars, 4)

        important_values <- data.frame(matrix(as.numeric(c(
          female_ests$`Af/Vf`[1,1], female_ests$`Af/Vf`[2,2], female_ests$`Af/Vf`[1,2],
          female_ests$`Ef/Vf`[1,1], female_ests$`Ef/Vf`[2,2], female_ests$`Ef/Vf`[1,2],
          male_ests$`Am/Vm`[1,1], male_ests$`Am/Vm`[2,2], male_ests$`Am/Vm`[1,2],
          male_ests$`Em/Vm`[1,1], male_ests$`Em/Vm`[2,2], male_ests$`Em/Vm`[1,2],
          female_cors$`solve(sqrt(I*Af)) %&% Af`[1, 2], female_cors$`solve(sqrt(I*Ef)) %&% Ef`[1, 2],
          male_cors$`solve(sqrt(I*Am)) %&% Am`[1, 2], male_cors$`solve(sqrt(I*Em)) %&% Em`[1, 2]
        )), nrow=1))
        cor_labels_in_model <<- c("af_2_1", "am_2_1", "ef_2_1", "em_2_1")
        cor_labels_in_result <<- c("Af_bivar", "Am_bivar", "Ef_bivar", "Em_bivar")
        colnames(important_values) <<- c("Af_v1", "Af_v2", "Af_bivar", "Ef_v1", "Ef_v2", "Ef_bivar",
                                         "Am_v1", "Am_v2", "Am_bivar", "Em_v1", "Em_v2", "Em_bivar",
                                         "rAf", "rEf", "rAm", "rEm")
        important_values[, c("Df_v1", "Df_v2", "Df_bivar", "Dm_v1", "Dm_v2", "Dm_bivar", "rDm", "rDf", "Cf_v1", "Cf_v2", "Cf_bivar", "Cm_v1", "Cm_v2", "Cm_bivar", "rCm", "rCf")] <<- NA
        if (print_status) {
          cat(paste0("Bivar model ", modelno," ",xname, " ~ ", yname, ":Testing genetic and environmental cors (7/8).\n"))
        }
        important_values <- Test_Cor_Paths(base_model=fitAE, base_model_name='AE',
                                           rA_paths = 'af_2_1', rE_paths= 'ef_2_1',
                                           output_df=important_values,
                                           rawresultsfilepath=rawresultsfilepath, p_cutoff=p_cutoff, sexlab='females')
        important_values <- Test_Cor_Paths(base_model=fitAE, base_model_name='AE',
                                           rA_paths = 'am_2_1', rE_paths= 'em_2_1',
                                           output_df=important_values,
                                           rawresultsfilepath=rawresultsfilepath, p_cutoff=p_cutoff, sexlab='males')
      }
      # Test significance of correlations in the model
      if (print_status) {
        cat(paste0("Bivar model ", modelno," ",xname, " ~ ", yname, ":Testing individual correlations (8/8).\n"))
      }
      important_values <- Test_Multiple_Correlations(base_model=fitAE, base_model_name='AE',
                                                     labels_to_test=cor_labels_in_model, label_names=cor_labels_in_result,
                                                     output_df=important_values, p_cutoff=p_cutoff, rawresultsfilepath=rawresultsfilepath)
    } else {
      drop_d <<- FALSE
      if (equal_sexes_ade) {
        female_ests <<- formatOutputMatrices(fitADE, c("Af","Df","Ef","Vf","Af/Vf","Df/Vf","Ef/Vf"), c("covAf","covDf","covEf","Var","stCovA","stCovD","stCovE"), selVars,4)
        female_cors <<- formatOutputMatrices(fitADE, c("solve(sqrt(I*Af)) %&% Af","solve(sqrt(I*Df)) %&% Df","solve(sqrt(I*Ef)) %&% Ef"), c("corA","corD","corE"), selVars, 4)
        important_values <<- data.frame(matrix(as.numeric(c(
          female_ests$`Af/Vf`[1,1], female_ests$`Af/Vf`[2,2], female_ests$`Af/Vf`[1,2], female_ests$`Df/Vf`[1,1], female_ests$`Df/Vf`[2,2], female_ests$`Df/Vf`[1,2], female_ests$`Ef/Vf`[1,1], female_ests$`Ef/Vf`[2,2], female_ests$`Ef/Vf`[1,2],
          female_cors$`solve(sqrt(I*Af)) %&% Af`[1, 2], female_cors$`solve(sqrt(I*Df)) %&% Df`[1, 2], female_cors$`solve(sqrt(I*Ef)) %&% Ef`[1, 2]
        )), nrow=1))
        cor_labels_in_model <<- c("a_2_1", "d_2_1", "e_2_1")
        if (base_model_used == "ADE") {
          cor_labels_in_result <<- c("Af_bivar", "Df_bivar", "Ef_bivar")
          colnames(important_values) <<- c("Af_v1", "Af_v2", "Af_bivar","Df_v1", "Df_v2", "Df_bivar", "Ef_v1", "Ef_v2", "Ef_bivar",
                                           "rAf", "rDf", "rEf"
          )
          important_values[, c("Cf_v1", "Cf_v2", "Cf_bivar", "Cm_v1", "Cm_v2", "Cm_bivar", "rCm", "rCf")] <<- NA
        } else {
          cor_labels_in_result <<- c("Af_bivar", "Cf_bivar", "Ef_bivar")
          colnames(important_values) <<- c("Af_v1", "Af_v2", "Af_bivar","Cf_v1", "Cf_v2", "Cf_bivar", "Ef_v1", "Ef_v2", "Ef_bivar",
                                           "rAf", "rCf", "rEf"
          )
          important_values[, c("Df_v1", "Df_v2", "Df_bivar", "Dm_v1", "Dm_v2", "Dm_bivar", "rDm", "rDf")] <<- NA
        }
      } else {
        female_ests <<- formatOutputMatrices(fitADE, c("Af","Df","Ef","Vf","Af/Vf","Df/Vf","Ef/Vf"), c("covAf","covDf","covEf","Var","stCovA","stCovD","stCovE"), selVars,4)
        female_cors <<- formatOutputMatrices(fitADE, c("solve(sqrt(I*Af)) %&% Af","solve(sqrt(I*Df)) %&% Df","solve(sqrt(I*Ef)) %&% Ef"), c("corA","corD","corE"), selVars, 4)
        male_ests <<- formatOutputMatrices(fitADE, c("Am","Dm","Em","Vm","Am/Vm","Dm/Vm","Em/Vm"), c("covA","covD","covE","Var","stCovA","stCovD","stCovE"), selVars,4)
        male_cors <<- formatOutputMatrices(fitADE, c("solve(sqrt(I*Am)) %&% Am","solve(sqrt(I*Dm)) %&% Dm","solve(sqrt(I*Em)) %&% Em"),  c("corA","corD","corE"), selVars, 4)
        important_values <<- data.frame(matrix(as.numeric(c(
          female_ests$`Af/Vf`[1,1], female_ests$`Af/Vf`[2,2], female_ests$`Af/Vf`[1,2], female_ests$`Df/Vf`[1,1], female_ests$`Df/Vf`[2,2], female_ests$`Df/Vf`[1,2], female_ests$`Ef/Vf`[1,1], female_ests$`Ef/Vf`[2,2], female_ests$`Ef/Vf`[1,2],
          male_ests$`Am/Vm`[1,1], male_ests$`Am/Vm`[2,2], male_ests$`Am/Vm`[1,2], male_ests$`Dm/Vm`[1,1], male_ests$`Dm/Vm`[2,2], male_ests$`Dm/Vm`[1,2], male_ests$`Em/Vm`[1,1], male_ests$`Em/Vm`[2,2], male_ests$`Em/Vm`[1,2],
          female_cors$`solve(sqrt(I*Af)) %&% Af`[1, 2], female_cors$`solve(sqrt(I*Df)) %&% Df`[1, 2], female_cors$`solve(sqrt(I*Ef)) %&% Ef`[1, 2],
          male_cors$`solve(sqrt(I*Am)) %&% Am`[1, 2], male_cors$`solve(sqrt(I*Dm)) %&% Dm`[1, 2], male_cors$`solve(sqrt(I*Em)) %&% Em`[1, 2]
        )), nrow=1))
        cor_labels_in_model <<- c("af_2_1", "am_2_1", "df_2_1", "dm_2_1", "ef_2_1", "em_2_1")
        if (base_model_used == "ADE") {
          cor_labels_in_result <<- c("Af_bivar", "Am_bivar", "Df_bivar", "Dm_bivar", "Ef_bivar", "Em_bivar")
          colnames(important_values) <- c("Af_v1", "Af_v2", "Af_bivar","Df_v1", "Df_v2", "Df_bivar", "Ef_v1", "Ef_v2", "Ef_bivar",
                                           "Am_v1", "Am_v2", "Am_bivar","Dm_v1", "Dm_v2", "Dm_bivar", "Em_v1", "Em_v2", "Em_bivar",
                                           "rAf", "rDf", "rEf", "rAm","rDm", "rEm"
          )
          important_values[, c("Cf_v1", "Cf_v2", "Cf_bivar", "Cm_v1", "Cm_v2", "Cm_bivar", "rCm", "rCf")] <<- NA
        } else {
          cor_labels_in_result <<- c("Af_bivar", "Am_bivar", "Cf_bivar", "Cm_bivar", "Ef_bivar", "Em_bivar")
          colnames(important_values) <- c("Af_v1", "Af_v2", "Af_bivar","Cf_v1", "Cf_v2", "Cf_bivar", "Ef_v1", "Ef_v2", "Ef_bivar",
                                           "Am_v1", "Am_v2", "Am_bivar","Cm_v1", "Cm_v2", "Cm_bivar", "Em_v1", "Em_v2", "Em_bivar",
                                           "rAf", "rCf", "rEf", "rAm","rCm", "rEm"
          )
          important_values[, c("Df_v1", "Df_v2", "Df_bivar", "Dm_v1", "Dm_v2", "Dm_bivar", "rDm", "rDf")] <<- NA
        }
      }
      if (print_status) {
        cat(paste0("Bivar model ", modelno," ",xname, " ~ ", yname, ":Testing genetic and environmental cors (7/8).\n"))
      }
      if (equal_sexes_ade) {
        if (base_model_used == 'ADE') {
          important_values <- Test_Cor_Paths(base_model=fitADE, base_model_name=base_model_used,
                                             rA_paths = c('a_2_1', 'd_2_1'), rE_paths= 'e_2_1',
                                             output_df=important_values,
                                             rawresultsfilepath=rawresultsfilepath, p_cutoff=p_cutoff, sexlab='females')
        } else {
          important_values <- Test_Cor_Paths(base_model=fitADE, base_model_name=base_model_used,
                                             rA_paths = 'a_2_1', rE_paths=c('d_2_1', 'e_2_1'),
                                             output_df=important_values,
                                             rawresultsfilepath=rawresultsfilepath, p_cutoff=p_cutoff, sexlab='females')
        }
      } else {
        if (base_model_used == 'ADE') {
          important_values <- Test_Cor_Paths(base_model=fitADE, base_model_name=base_model_used,
                                             rA_paths = c('af_2_1', 'df_2_1'), rE_paths= 'ef_2_1',
                                             output_df=important_values,
                                             rawresultsfilepath=rawresultsfilepath, p_cutoff=p_cutoff, sexlab='females')
          important_values <- Test_Cor_Paths(base_model=fitADE, base_model_name=base_model_used,
                                             rA_paths = c('am_2_1', 'dm_2_1'), rE_paths= 'em_2_1',
                                             output_df=important_values,
                                             rawresultsfilepath=rawresultsfilepath, p_cutoff=p_cutoff, sexlab='males')
        } else {
          important_values <- Test_Cor_Paths(base_model=fitADE, base_model_name=base_model_used,
                                             rA_paths = 'af_2_1', rE_paths=c('df_2_1', 'ef_2_1'),
                                             output_df=important_values,
                                             rawresultsfilepath=rawresultsfilepath, p_cutoff=p_cutoff, sexlab='females')
          important_values <- Test_Cor_Paths(base_model=fitADE, base_model_name=base_model_used,
                                             rA_paths = 'am_2_1', rE_paths=c('dm_2_1', 'em_2_1'),
                                             output_df=important_values,
                                             rawresultsfilepath=rawresultsfilepath, p_cutoff=p_cutoff, sexlab='males')
        }
      }
      if (print_status) {
        cat(paste0("Bivar model ", modelno," ",xname, " ~ ", yname, ":Testing individual correlations (8/8).\n"))
      }
      important_values <- Test_Multiple_Correlations(base_model=fitADE, base_model_name=base_model_used,
                                                     labels_to_test=cor_labels_in_model, label_names=cor_labels_in_result,
                                                     output_df=important_values, p_cutoff=p_cutoff, rawresultsfilepath=rawresultsfilepath)
    }
  }


  output <- cbind(data.frame(v1=xname, v2=yname,
                             base_model_used=base_model_used, base_model_sig=base_model_sig, base_model_sig_p=base_model_sig_p,
                             equal_sexes_ade=equal_sexes_ade, equal_sexes_ade_p=equal_sexes_ade_p,
                             drop_d=drop_d, drop_d_p=drop_d_p,
                             could_drop_a=could_drop_a, could_drop_a_p=could_drop_a_p), important_values)
  if (print_status) {
    cat(paste0("Bivar model ", modelno," ",xname, " ~ ", yname, ":FINISHED.\n"))
  }
  return(output)
}

##############################
###     SCRIPT START       ###
##############################
# Load Data
data <- read.csv(datapath, stringsAsFactors=FALSE)

# Verify data
if (!("zyg" %in% colnames(data))) {
  stop("Column 'zyg' not found in the data file")
}
missing_cols <- c()
for (x in c(xvars, yvars)) {
  if (include_sibs) {subjects <- c("twin1", "twin2", "sibm1", "sibf2")} else {subjects <- c("twin1", "twin2")}
  for (s in subjects) {
    if (!(paste0(x, s) %in% colnames(data))) {
      missing_cols <- c(missing_cols, paste0(x, s))
    }
  }
}
if (length(missing_cols) > 0) {
  stop(paste0("The following phenotype columns were expected but not found: ", paste0(missing_cols, collapse = ', ')))
}

if (print_status) {
  cat("Data seems okay, starting analyses...\n")
}
### Start all analyses
final_result_initialized <- FALSE
for (ii in 1:length(xvars)) {
  for (jj in 1:length(yvars)) {
    x <- xvars[ii]
    y <- yvars[jj]
    modelnumber <- paste0('(', (ii-1)*length(yvars) + jj, '/', length(xvars)*length(yvars),')')
    if (!(final_result_initialized)) {
      final_result <- main(data, x, y, modelnumber)
      final_result_initialized <- TRUE
    } else {
      result <- main(data, x, y, modelnumber)
      # Ensure all column names are shared
      if (any(!(colnames(final_result) %in% colnames(result)))) {
        result[, colnames(final_result)[!(colnames(final_result) %in% colnames(result))]] <- NA
      }
      if (any(!(colnames(result) %in% colnames(final_result)))) {
        final_result[, colnames(result)[!(colnames(result) %in% colnames(final_result))]] <- NA
      }
      # Ensure column order is equal
      result <- result[, colnames(final_result)]
      # Append
      final_result <- rbind(final_result, result)
    }
  }
}

write.csv(final_result, output_file_name, row.names=FALSE)