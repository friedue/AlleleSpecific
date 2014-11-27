ReadingInFeatureCountData.edgeR.genome <- function(FileName,  ##<< path to the featureCount output (can also be a vector with more than one file)
                                            Sample, ##<< Sample name, e.g. "Input"
                                            Cell ##<< CellType, e.g. "femES"
                                            )
{
  ### based on featureCount output, this function reads in the entire table and prepares it for edgeR
  ### difference to original function, e.g. found in f_edgeRWithFeatureCounts_2014-11.R: selection of columns instead of merging count results from different files
  ### does not return a composite DF, but only one for the selected column from the composite featureCount table
  
  library(data.table)
  
  #In = lapply(FileName, function(x){ z =  unique(fread(x, header=T))                                     
  #                                   setkey(z, Geneid) # requires the presence of a column named "Geneid" in the input (default of featureCounts)
  #                                   return(z)})
  
  In.out = list()
  
  In.select = grep(paste(Sample,".*", Cell,"_[129S1|CASTEiJ]", sep = ""),
                   names(fread(paste("head -3", FileName), header=T, skip = 1, sep = "\t"))) # extract column number with information for the indicated Sample
  In = as.data.frame(fread(FileName, header=T, skip = 1, sep = "\t", select = c(1:6, In.select)))
  In.out$Bins = In[c(1:4)]
  
  ### Prepare for EdgeR
  ### probably not very generic due to the string replacement for the namings
  ### (requires a rather strict naming order)
  ### and the format expectations based on featureCount output
  
  names(In) = gsub(".*BAMs/([a-zA-Z0-9]*)(A|B)_(femES|femNPC)_(.*)(129S1|CASTEiJ).*\\.bam",  "\\1_\\2_\\3_\\5", names(In))
  row.names(In) = In$Geneid
  In.out$forEdgeR = subset(In, select = grepl(gsub("(.*)[A|B]_.*", "\\1", Sample), names(In)))
  
  return(In.out)
}


PreparingFactorLists <- function(Names, ##<< Names should be the names of the factors, e.g. c("MOF", "MSL1", "MSL2")
                                 ReferenceLevel ##<< should indicate the factor that will serve as the "baseline" for all comparisons, e.g. "MOF"
){
  if(is.character(Names) & length(ReferenceLevel)==1){
    Out = factor(Names, levels= unique(Names))
    Out = relevel(Out, ref = ReferenceLevel)
  }else{
    stop("Check the input supplied to PreparingFactorLists()")
  }
  
  return(Out)
}

RunningEdgeR.dispEst <- function(EdgeR.dge, designMat){
  ## should be run _after_ dgeGeneration that generated the EdgeR.dge list
  ## designMat should contain the coefficients for the model
  EdgeR.dge <- estimateGLMCommonDisp(EdgeR.dge, designMat)
  EdgeR.dge <- estimateGLMTrendedDisp(EdgeR.dge, designMat)
  EdgeR.dge <- estimateGLMTagwiseDisp(EdgeR.dge, designMat)
  
  return(EdgeR.dge)
}

RunningEdgeR <- function(Counts.DF, ##<< DF of counts, should not contain any additional columns except the counts; row.names should be set
                         gtReference = "CASTEiJ",
                         Group , ##<< vector or factor, e.g. gsub("(.*)_([0-9AB]*)_(.*)_(.*)", "\\4", names(Counts.DF)
                                       ##<< indicating experimental group/condition for each sample/library:
                                       ##<< e.g.rep("1", length(Counts[c(2:length(Counts))])
                         calculateSizeFactor = TRUE, ##<< whether a library size factor should be calculated by edgeR
                                       ## can be BOOLEAN or a vector of numbers or
                                       ## a data.frame with row.names matching those of the samples and a column called norm.factors (= result of edgeR's DGEList generation)
                         ExpDesign = "Genotype",
                         PairwiseTest = FALSE
){
  require(edgeR)
  
  dge <-DGEList(Counts.DF, group = Group)
  
  ### filtering
  # countsPerMillion = cpm(dge)
  ### tag must have 2 CPM in at least 2 samples
  #keep = rowSums(countsPerMillion > 2) >= 2 
  keep = rowSums(Counts.DF > 2) >= 2 ## changed this to the raw counts because I was losing too many reads with filtering based on cpm
  dge = dge[keep,]
  
  if(is.logical(calculateSizeFactor)){
    if(calculateSizeFactor){
      dge<-calcNormFactors(dge)}
  }
  else if(is.numeric(calculateSizeFactor) & (length(calculateSizeFactor) == nrow(dge$samples))){
    dge$samples$norm.factors = calculateSizeFactor # This does not pay attention to the correct assignment; it's simply order-based
  }
  else if(is.data.frame(calculateSizeFactor) & ("norm.factors" %in% names(calculateSizeFactor))){
    dge$samples$norm.factors = sapply(gsub("(.*)_([0-9AB]*)_(.*)_(.*)", "\\1_\\2_\\3", row.names(dge$samples)), #need to extract the sample name without the genotyp info
                                      function(x) dge$samples[x,]$norm.factors = calculateSizeFactor[x,]$norm.factors) # fortunately, the row indexing works with just parts of the nam
  }
  else{
    print("The option indicated for size factor calculation is neither logical nor numeric nor a data.frame with a column named 'norm.factors'.
             No size factor was calculated.")
  }
  
  #return(dge)
  
  #Replicates = gsub("(.*)_([0-9AB]*)_(.*)", "\\2", names(Counts.DF))
  
  ### Factor lists
  #Replicate = PreparingFactorLists(gsub("(.*)_([0-9AB]*)_(.*)", "\\2", names(Counts.DF)), ReferenceLevel = 1)
  Genotype = PreparingFactorLists(gsub("(.*)_([0-9AB]*)_(.*)_(.*)", "\\4", names(Counts.DF)), ReferenceLevel = gtReference)
  
  dge$samples$group = relevel(dge$samples$group, ref = gtReference)
  ### making the design matrix
  design <- model.matrix(reformulate(ExpDesign))
  
  ### estimate dispersion
  dge = RunningEdgeR.dispEst(EdgeR.dge = dge, designMat=design)
  
  ### GLM fitting
  ModelFit <- glmFit(dge, design)
  
  Out = list(kept = keep, DGE = dge, Model = ModelFit, Design = design)
  ### exact test
  if(PairwiseTest){
    Out[["ExactTest"]] <- exactTest(dge)
    
    ### GLM 
    Out[["GLMtest"]] <- glmLRT(ModelFit, coef=2)
  }
  
  return(Out)
}

AIFromlog2FC.df <- function(DF ##<< must contain a column named "logFC"
){
  # Calculate the AI score based on the foldChange calculated by edgeR.
  if("logFC" %in% names(DF)){
    FC = 2^DF$logFC
    DF$AI = (FC-1)/(FC+1)
  }else{
    stop("The provided data.frame does not contain the accessor 'logFC'")
  }
  return(DF)
}
