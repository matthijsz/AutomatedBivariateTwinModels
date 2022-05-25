# AutomatedBivariateTwinModels
#### Automate fitting of an arbitrary number of bivariate twin models using OpenMx  
This script is ment to make it easy to fit a bunch of bivariate models and get the important estimates without any interactive hassle. It is not ideal for finely detailed info/output/adaptations to single models, if that is the goal, extract the model specifications from the BivariateTwinModelsResources.R and the labels and stuff from AutomatedTwinModel.R and make another version... Then again if you want really fine control you may just want to make your own version.

The main script is AutomatedTwinModel.R which relies heavily on the functions in the accompanying BivariateTwinModelsResources.R. This script was made by me, based on script I got from a colleague, who in turn based it on scripts from other sources. So use at your own discretion. If you do use it please make sure you cite the original OpenMx publication(s).

## Usage
First change global script settings in AutomatedTwinModel.R to your liking, save your changes Then run AutomatedTwinModel.R.

## Script settings
### Basic settings
```xvars <- c("x1", "x2", "x3")```  
xvars should be a vector of the base names of x variables, if there are 100s I recommend storing them in a file and reading them from there to keep things clean, e.g. `xvars <- read.csv("BivarXnames.csv")[, "x"]`
```yvars <- c("y1", "y2", "y3")```  
yvars should be a vector of the base names of y variables, if there are 100s I recommend storing them in a file and reading them from there to keep things clean, e.g. `yvars <- read.csv("BivarYnames.csv")[, "x"]`

```wd <- "."```  
Working directory, this directory should also contain "BivariateTwinModelsResources.R". Leave this at "." if everything is in the same root directory

```datapath <- "DataWide.csv"```   
Path to the dataset file, see below.  
This file should have a columnd called 'zyg' to indicate the zygosity of the twin pair, according to the following:  
  1: MZM  
  2: DZM  
  3: MZF  
  4: DZF  
  5: DOS (Note the male twin in a DOS pair should always be twin1, the female twin should be twin2)  
The file should also adhere to the following naming scheme:  
```[variable]twin[twin-number]```  
For example:  x1twin1, x2twin1, x1twin2, x2twin2, y1twin1, y2twin1, y1twin2, y2twin2, etc.    
Order of these columns is not important, column name is crucial  
If 'include_sibs' is TRUE the file should have columns  
```[variable]sibm1 and [variable]sibf2```  
for a brother and sister of every variable respectivelly.

```output_file_name <- "Bivariate_model_results.csv"```  
Name of the final output file

### Model settings

```include_sibs <- FALSE```  
Expand the standard twinmodel to a twin-sibling model (note: requires different data, see below)

```preferred_model <- "ACE"```   
The model (ADE or ACE) that is preferred, see below.  
This script will start with a Saturated model and fit both an ADE and ACE model.
If both of these fit significantly worse, or if neither of these fit significantly worse, the script will revert to 'preferred_model'
Note that to spare myself a migraine the C component in an ACE model is internally still labelled 'D' (so I didn't have to write everything twice)
The labels in the final output file are changed appropriately, but nothing else. Just be aware of that when running this script manually, or checking the extensive raw results.

```p_cutoff <- 0.01```   
P-value cut-off for significant, for a standard multiple-testing correction use something like:
```p_cutoff <- 0.05/(length(xvars)*length(yvars))```  
though this is likely an overcorrection as then you assume all x and y combinations are totally independant (which is unlikely)

```save_extensive_raw_results <- FALSE```   
Should more extensive results be saved to a text file for each combination of phenotypes?
Note this setting will save 1 text file per combination of x and y, so it can become a bit messy.

```samesex_only <- FALSE```    
Should the script only use same-sex twin pairs?

```univar <- FALSE```    
Run univariate (in this case y-variables don't matter)

```force_AE <- FALSE```  
Should the script enforce the use of an AE model, regardless of how it fits
You can use this to, for example, ensure rA estimates across all phenotypes are comparable.
Though it of course comes with the note that it can result in sub-optimal models being used.

### Other settings

```print_status <- TRUE```    
Print updates on what the script is doing, this will print a line of which models are currently being fitted on which variables. 

```openmx_n_threads <- 2```   
CPU threads passed to OpenMx (it is still unclear to me if this even makes a meaningfull difference, but hey, here you go...)

```openmx_optimizer <- "NPSOL"```  
OpenMx optimizer, see the openmx documentation
```openmx_quiet <- TRUE```    
Forcably silence OpenMx and optimizer updates


## Pipeline
This script will do the following, in this order:
  1. Fit a saturated model.
  2. Fit an ADE model.
  3. Fit an ACE model.
  4. Select ADE or ACE model based on comparison with saturated model (or divert to preferred_model, see settings).      The selected model will then be used from here on out.
  5. Fit a model with equal A, D/C, E for both sexes. If the fit is not significantly (see settings) worse, this model will then be used from here on out.
  6. Fit an AE model.
  7. Fit an CE/DE model.
  8. Compare AE and CE/DE model to the current best model. If the fit of the AE model is not significantly worse, this model will then be used from here on out.
  9. Test genetic and environmental correlations globally, what happens here depends on the model and if sexes are equal:
      - With ACE model: test rA as genetic correlation, test rC and rE combined as environmental correlation
      -  With ADE model: test rA and rD combined as genetic correlation, test rE as environmental correlation
      - With AE model: test rA as genetic correlation, test rE as environmental correlation
      - If sexes are not equal, each test is performed twice (one for females, one for males)
 10. Test each correlation individually, this test is a bit redundant with AE models, but to maintain consistent column names and stuff I've left it as is.
 Intermediate "result" files can be saved that contain more detail about each model comparison, and more detailed information of the final model, see settings below

## Result file
The final result file is one dataframe with 1 row per combination of x and y with the following columns:
   - v1: Name of variable 1
   - v2: Name of variable 2
   - base_model_used: Which base model was used (ADE or ACE)
   - base_model_significant: Did the base model fit? (i.e. comparison_p of base model with saturated was below cutoff)
   - base_model_significant_p: p_value of the model_comparison used for the above
   - equal_sexes: Did the equal sexes model fit (i.e. comparison_p of equal sexes with ADE/ACE was below cutoff, I'm gonna stop adding this now)
   - equal_sexes_p: p_value of the model_comparison used for the above
   - drop_d: Can D be dropped?
   - drop_d_p: p_value of the model_comparison used for the above
   - could_drop_a: Could A have been dropped (instead of D)
   - could_drop_a_p: p_value of the model_comparison used for the above
  And the following columns for males and females (I'm only adding females to this list as I think the naming scheme is pretty obvious)
   - Af_v1: % variance explained by A in v1 for females
   - Af_v2: % variance explained by A in v2 for females
   - Af_bivar: bivariate A path value for females
   - Cf_v1: % variance explained by C in v1 for females
   - Cf_v2: % variance explained by C in v2 for females
   - Cf_bivar: bivariate C path value for females
   - Df_v1: % variance explained by D in v1 for females
   - Df_v2: % variance explained by D in v2 for females
   - Df_bivar: bivariate D path value for females
   - Ef_v1: % variance explained by E in v1 for females
   - Ef_v2: % variance explained by E in v2 for females
   - Ef_bivar: bivariate E path value for females
   - rAf: Additive genetic correlation for females
   - rCf: Common environmental correlation for females
   - rDf: Non-additive genetic correlation for females
   - rEf: Unique environmental correlation for females
   - females_genetic_pleoitropy_paths_significant: Are the genetic correlation paths (see step 9 above) significant for females?
   - females_genetic_pleoitropy_paths_significant_p: p-value, you get the idea
   - females_rE_paths_significant: Are the environmental correlation paths (see step 9 above) significant for females?
   - females_rE_paths_significant_p: p-value, you get the idea
   - females_all_correlations_significant: Are both genetic and environmental correlations independantly significant for females?
   - Af_bivar_significant_p: is path Af_bivar significant?
   - Cf_bivar_significant_p: is path Cf_bivar significant?
   - Df_bivar_significant_p: is path Df_bivar significant?
   - Ef_bivar_significant_p: is path Ef_bivar significant?  
 If sexes are equal, the female estimates contain the combined-sexes estiamtes, so the proper labelling for females would be "females (or all sexes)"  
 I did this on purpose to prevent making the output datafile even bigger than it already is.  
 Also, the order of these columns is not fixed, so if you have many varying results (AE, ADE and ACE models) you may want to order them after the fact.
