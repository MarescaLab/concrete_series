"DESeq2_rlog_norm_qiime"<- function(input_path, out_path, DESeq_negatives_to_zero=NULL) {
    library(metagenomeSeq)
    library(DESeq2)
    library(biomformat)
    foo = read_biom(input_path)
    x = as(biom_data(foo), "matrix")
    # avoid zeros
    x = x + 1
    sampleTable <- data.frame(sampleName=colnames(x))
    #Add mock design: these should not influence normalization - just required for DESeqDataSetFromMatrix
    dds <- DESeqDataSetFromMatrix(x, sampleTable, design=~1)
    vsmat = assay(rlogTransformation(dds))
    if (!is.null(DESeq_negatives_to_zero)) {
      vsmat[vsmat<0]=0
    }
    vsmat = newMRexperiment(vsmat)
    colnames(vsmat) <- c(colnames(x))
    write_biom(MRexperiment2biom(vsmat), out_path)
}

"CSS_norm_qiime_log" <- function(input_path, out_path, output_CSS_statistics=NULL) {
    library(metagenomeSeq)
    library(DESeq2)
    library(biomformat)
    obj = load_biom(input_path)
    p = cumNormStat(obj)
    obj = cumNorm(obj, p = p)
    if (!is.null(output_CSS_statistics)) {
      exportStats(obj, p=p, file = file.path(output_CSS_statistics))
    }
    write_biom(MRexperiment2biom(obj, norm=TRUE, log=TRUE), out_path)
}

"CSS_norm_qiime_counts" <- function(input_path, out_path, output_CSS_statistics=NULL) {
  library(metagenomeSeq)
  library(DESeq2)
  library(biomformat)
  obj = load_biom(input_path)
  p = cumNormStat(obj)
  obj = cumNorm(obj, p = p)
  if (!is.null(output_CSS_statistics)) {
    exportStats(obj, p=p, file = file.path(output_CSS_statistics))
  }
  write_biom(MRexperiment2biom(obj, norm=TRUE, log=FALSE), out_path)
}