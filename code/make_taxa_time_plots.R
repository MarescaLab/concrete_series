process_for_plot <- function(otu_table,taxonomy_tsv,metadata, transformation){
  
  library(qiime2R)
  library(tidyverse)
  library(readr)
  library(data.table)
  library(reshape2)
  library(scales)
  library(lubridate)
  library(microbiome)
  
  #rel_abund <- function(x){x/sum(x)}
  
  otuTable <- read_qza(otu_table) %>% .$data %>% as.matrix() %>%  microbiome::transform(transform = transformation)
  
  tax_ranks <- c("kingdom","phylum","class","order","family","genus","species")
  
  taxonomy <- read_delim(taxonomy_tsv, "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
    select(-confidence) %>% 
    mutate(for_split = taxonomy) %>% 
    separate(for_split,tax_ranks,sep = "; ") %>% 
    mutate_at(tax_ranks,str_remove,pattern = "[a-z]__")
    
    #For default SILVA string parsing
    #separate(for_split,tax_ranks,sep = ";") %>% 
    #mutate_at(tax_ranks,str_remove,pattern = "D_[0-9]__")
    
  met <- read_delim(metadata, "\t", escape_double = FALSE, trim_ws = TRUE)
  
  summary_file <- t(otuTable) %>% 
    melt(value.name = "rel_abund") %>% 
    dplyr::rename("#SampleID" = Var1, "#OTU ID" = Var2) %>% 
    left_join(taxonomy) %>% 
    left_join(met) %>%
    group_by(`#SampleID`) %>%
    ungroup()
  
  return(summary_file)  
}

series_plot <- function(summary_table,taxa,agglom = FALSE,smooth_method = "loess",se = TRUE,point_alpha = 4/10,ncol = 2,facet_scales,temp_shade,filter_zeros){
  
  #Set taxa as ordered factors, so that they'e plotted in the right order
  taxa <- taxa %>% factor(.,levels = .,ordered = TRUE)
  
  #make table of interesting taxa function
  
  tax_ranks <- c("kingdom","phylum","class","order","family","genus","species","#OTU ID")  
  
  plot_table <- summary_table %>% mutate(otu_name = .$`#OTU ID`) %>%
    filter_at(tax_ranks, any_vars(. %in% taxa)) %>% 
    {if (filter_zeros==TRUE) filter(., rel_abund > 0) else .} %>% 
    select(`#SampleID`,ASR,Months,Collected,rel_abund,taxonomy,temp_avg,otu_name,one_of(tax_ranks)) %>%
    gather(key="Taxa_type",value = "Taxa",-c(`#SampleID`,ASR,Months,temp_avg,Collected,rel_abund,taxonomy,otu_name)) %>%
    filter(Taxa %in% taxa) %>% 
    mutate(for_split = taxonomy) %>% 
    separate(for_split,tax_ranks,sep = "; ") %>% 
    mutate_at(tax_ranks,str_remove,pattern = "[a-z]__") %>%
    mutate(og_taxa = Taxa, 
           #Taxa = ifelse(Taxa_type == "species", paste(genus, species),Taxa), 
           Taxa = ifelse(Taxa_type == "#OTU ID", paste0(taxonomy,"\n",otu_name),Taxa), 
           "#OTU ID" = otu_name, 
           tax_and_level= paste(Taxa," (",substr(Taxa_type, start = 1, stop = 1),")",sep = ""), 
           s_date = mdy(Collected)) %>%
    select(-otu_name) #%>%
    #mutate(Taxa = factor(og_taxa,levels = taxa,ordered = TRUE))
  
  temp_data <- plot_table %>% select(samp_day = s_date,temp_avg) %>% ungroup() %>% unique()

  months_between_label = 3
  
  #Make plots
  if (agglom == FALSE){
    plot <- plot_table %>% 
      ggplot(aes(s_date,rel_abund,color = ASR, shape = ASR, fill = NULL, linetype = ASR, label = `#SampleID`)) +
      #geom_rect(aes(xmin=date - days(30),xmax=date,ymin=min(rel_abund),ymax=max(rel_abund), fill=temp_avg, color = NULL),alpha = 0.1,inherit.aes = FALSE) +
      {if(temp_shade == TRUE)geom_rect(data = temp_data, aes(xmin=samp_day - days(30),xmax=samp_day,ymin=0,ymax=1, fill=temp_avg, color = NULL),alpha = 0.2,inherit.aes = FALSE)} +
      scale_fill_viridis_c() +
      geom_point(alpha = point_alpha, position = "jitter", size = 2) + 
      geom_smooth(method = smooth_method, se = se) +
      scale_y_log10() +
      #facet_wrap(facets = ~ `#OTU ID` + taxonomy, ncol = 3,scales = "free_y", labeller = label_context) + 
      theme(legend.key.size = unit(.5,"inches"), legend.title=element_text(size=16) , legend.text=element_text(size=16)) +
      theme(axis.text=element_text(size=16), axis.title=element_text(size=24)) +
      theme(strip.text.x = element_text(size = 24, colour = "black", angle = 0)) +
      labs(x= NULL, y= "Relative abundance (log)", shape = "ASR", color = "ASR", fill = "Temp. (F)", alpha = "Temp. (F)", linetype = "ASR") +
      theme(legend.position = c(0.9, 0.6),legend.background = element_rect(fill = "white", colour = NA)) +
      theme_bw() +
      scale_shape_manual(values=c(17, 21), name = "ASR") +
      scale_color_manual(values = c("darkorange2","grey40"), name = "ASR") +
      scale_x_date(labels=date_format("%b '%y"),date_breaks = paste(months_between_label,"months"),limits = c(min(mdy(summary_table$Collected)),max(mdy(summary_table$Collected)))) +
      theme(axis.text.x = element_text(angle=45,vjust = 1, hjust = 1)) +
      facet_wrap(~ tax_and_level,ncol = ncol,scales = facet_scales)
    print(plot)
  }
  else{
    plot <- plot_table %>% group_by(ASR,s_date,`#SampleID`,tax_and_level,temp_avg) %>% summarise(rel_abund = sum(rel_abund)) %>% 
      ggplot(aes(s_date,rel_abund,color = ASR, shape = ASR, fill = NULL, linetype = ASR, label = `#SampleID`)) +
      #geom_rect(aes(xmin=date - days(30),xmax=date,ymin=min(rel_abund),ymax=max(rel_abund), fill=temp_avg, color = NULL),alpha = 0.1,inherit.aes = FALSE) +
      {if(temp_shade == TRUE)geom_rect(data = temp_data, aes(xmin=samp_day - days(30),xmax=samp_day,ymin=0,ymax=1, fill=temp_avg, color = NULL),alpha = 0.2,inherit.aes = FALSE)} +
      scale_fill_viridis_c() +
      geom_point(alpha = point_alpha, position = "jitter", size = 2) + 
      geom_smooth(method = smooth_method, se = se) +
      scale_y_log10() +
      #facet_wrap(facets = ~ `#OTU ID` + taxonomy, ncol = 3,scales = "free_y", labeller = label_context) + 
      theme(legend.key.size = unit(.5,"inches"), legend.title=element_text(size=16) , legend.text=element_text(size=16)) +
      theme(axis.text=element_text(size=16), axis.title=element_text(size=24)) +
      theme(strip.text.x = element_text(size = 24, colour = "black", angle = 0)) +
      labs(x= NULL, y= "Relative abundance (log)", shape = "ASR", color = "ASR", fill = "Temp. (F)", alpha = "Temp. (F)", linetype = "ASR") +
      theme(legend.position = c(0.9, 0.6),legend.background = element_rect(fill = "white", colour = NA)) +
      theme_bw() +
      scale_shape_manual(values=c(17, 21), name = "ASR") +
      scale_color_manual(values = c("darkorange2", "grey40"), name = "ASR") +
      scale_x_date(labels=date_format("%b '%y"),date_breaks = paste(months_between_label,"months"),limits = c(min(mdy(summary_table$Collected)),max(mdy(summary_table$Collected)))) +
      theme(axis.text.x = element_text(angle=45,vjust = 1, hjust = 1)) +
      facet_wrap(~ tax_and_level,ncol = ncol,scales = facet_scales)
    print(plot)
  }
}



#' Plot taxa observed in the concrete cylinder series.
#' 
#' @param taxa The taxa to plot.
#' @param taxonomy The taxonomy file to use.
#' @param metadata The metadata file to use.
#' @param aggolm Should subtaxa be agglomerated prior to plotting.
#' @param smooth_method The smoothing method to apply.
#' @param se Should standard error shading be shown.
#' @param point_alpha The alpha value of points.
#' @param ncol Number of faceting columns.
#' @param facet_scales Scaling applied to facets, eg. same for all or free x/y.
#' @param temp_shade Visualize temerature throughout the series by shading the background.
#' @param transformation Transformation applied to OTU table prior to plotting.
#' @return Plots \code{taxa} relative abundance over time.
#' @examples
#' plot.taxa("Bacillus")
plot.taxa <- function(taxa,otu_table="d:/Google Drive/Documents/UDel/Research/Concrete_analysis/data/processed/biom_tables/concrete_only/all_w_taxonomy_decontam_table.biom_concrete.biom_table.qza",
                      taxonomy = "d:/Google Drive/Documents/UDel/Research/Concrete_analysis/data/processed/silva_taxonomy.tsv",
                      metadata = "d:/Google Drive/Documents/UDel/Research/Concrete_analysis/data/sample-metadata.tsv",
                      agglom=FALSE,
                      smooth_method = "loess",
                      se = TRUE,
                      point_alpha = 4/10,
                      ncol = 2,
                      facet_scales = "free_y", 
                      temp_shade = FALSE, 
                      transformation = "compositional",
                      filter_zeros = TRUE){
  
  summary <- process_for_plot(otu_table,taxonomy,metadata,transformation)
  series_plot(summary,taxa,agglom,smooth_method,se,point_alpha,ncol,facet_scales,temp_shade,filter_zeros)
}

###Usage: plot_taxa(taxa = c("Psychrobacter"),agglom = TRUE) + geom_smooth(method = "lm")







#Parameter definitions for testing
# met_col = "#SampleID"
# met_conditions = c('Type == "Concrete" & sample_id == "021715M_2"')
# otu_table = "d:/Google Drive/Documents/UDel/Research/Concrete_analysis/data/processed/biom_tables/decontaminated/all_w_taxonomy_decontam_table.biom_table.qza"
# tax_file = "d:/Google Drive/Documents/UDel/Research/Concrete_analysis/data/qiime2/silva_taxonomy.qza"
# metadata = "d:/Google Drive/Documents/UDel/Research/Concrete_analysis/data/sample-metadata.tsv"
# out_file = "c:/Users/eande/Desktop/Concrete_bacteria_overview.png"
# minimum_count = 50

series.metacoder <- function(met_conditions = 'Type == "Concrete"', 
                             otu_table = "d:/Google Drive/Documents/UDel/Research/Concrete_analysis/data/processed/biom_tables/decontaminated/all_w_taxonomy_decontam_table.biom_table.qza",
                             tax_file = "d:/Google Drive/Documents/UDel/Research/Concrete_analysis/data/qiime2/silva_taxonomy.qza",
                             metadata = "d:/Google Drive/Documents/UDel/Research/Concrete_analysis/data/sample-metadata.tsv",
                             out_file = "c:/Users/eande/Desktop/metacoder.png",
                             w = 10,
                             h = 10,
                             min_rel = 0.001,
                             otus_to_exclude = NA,
                             otus_to_include = NA){
  
  library(tidyverse)
  library(qiime2R)
  library(metacoder)
  library(viridis)
  
  
  if (is.na(otus_to_exclude)) {
    otu_mat <- qiime2R::read_qza(otu_table,tmp = "c:/TEMP/qiime2r") %>% .$data %>% 
      as.data.frame() %>% 
      mutate(otu_id = rownames(.))
  }
  if (!is.na(otus_to_exclude)) {
    otu_mat <- qiime2R::read_qza(otu_table,tmp = "c:/TEMP/qiime2r") %>% .$data %>% 
      as.data.frame() %>% 
      mutate(otu_id = rownames(.)) %>%
      filter(!otu_id %in% otus_to_exclude)
  }
  if (!is.na(otus_to_include)) {
    otu_mat <- qiime2R::read_qza(otu_table,tmp = "c:/TEMP/qiime2r") %>% .$data %>% 
      as.data.frame() %>% 
      mutate(otu_id = rownames(.)) %>%
      filter(otu_id %in% otus_to_include)
  }
  
  taxonomy <- qiime2R::read_qza(tax_file,tmp = "c:/TEMP/qiime2r") %>% .$data %>% 
    as.data.frame() %>% 
    mutate(otu_id = .$Feature.ID,lineage= .$Taxon) %>% 
    select(otu_id,lineage)
  
  met_samples <- read_delim(metadata, "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
    select(sample_id = `#SampleID`,everything()) %>% 
    filter(sample_id %in% colnames(otu_mat)) %>% 
    mutate(all_samples = "all_samples") %>%
    filter_(met_conditions)
  
  met_otus <- otu_mat %>% 
    select(one_of(met_samples$sample_id),otu_id) %>% 
    left_join(taxonomy) %>% 
    select(otu_id,lineage,everything()) %>% 
    filter(!is.na(lineage),lineage!="Unassigned", !str_detect(lineage,"Chloroplast"), !str_detect(lineage,"Mitochondria"))
  
  
  #Make the MetaCodeR obj
  ##Greengenes class_regex = "^(.+)__(.*)$"
  ##SILVA class_regex = "^D_(.+)__(.*)$"  
  
  obj <- parse_tax_data(met_otus, class_cols = "lineage", class_sep = ";",
                        class_key = c(tax_rank = "info", tax_name = "taxon_name"),
                        class_regex = "^D_(.+)__(.*)$")
  
  #filtering otus with no reads
  no_reads <- rowSums(obj$data$tax_data[, met_samples$sample_id]) == 0
  sum(no_reads)
  obj <- filter_obs(obj, "tax_data", ! no_reads, drop_taxa = TRUE)
  
 
  #summing per-taxon counts
  obj$data$rel_abd <- calc_obs_props(obj, "tax_data")
  obj$data$tax_abund <- calc_taxon_abund(obj, "rel_abd",groups = met_samples$all_samples)
  
  #To get in final relative abundance form
  obj$data$tax_abund$all_rel <- obj$data$tax_abund$all_samples / obj$data$tax_abund$all_samples[1]
  
  #To create human readable table of relative abundance of each taxonomic group
  # prop_w_tax <- obj$data$tax_abund %>% 
  #   left_join(obj$data$class_data %>% select(-input_index), by = "taxon_id") %>% 
  #   distinct()
  # 
  
  #Make the heat tree
  obj %>%
  taxa::filter_taxa(all_rel > min_rel) %>%
  heat_tree(node_label = taxon_names,
            node_size = all_rel,
            node_size_range = c(0.0005,0.05),
            node_label_size_range = c(.015,.025),
            node_size_axis_label = "OTU count",
            initial_layout = "reingold-tilford", layout = "davidson-harel",
            overlap_avoidance = 10,
            node_label_max = 60 ,
            node_color = all_rel,
            node_color_range = c("gray80","gray80","gray80"),
            node_color_axis_label = "Relative abundance") +
    ggsave(out_file, units = "in", dpi = 300, type = "cairo", width = w, height = h)
  
}




#Make plot of relative abundance in each precursor material


plot_precursor_abund <- function(taxa, 
                                 include_concrete = FALSE, 
                                 split_asr = TRUE,
                                 agglom = TRUE,
                                 transformation = "compositional",
                                 otu_table = "data/processed/biom_tables/decontaminated/all_w_taxonomy_decontam_table.biom_table.qza", 
                                 taxonomy_tsv = "data/processed/silva_taxonomy.tsv",
                                 metadata = "data/sample-metadata.tsv",
                                 point_alpha = 0.4,
                                 facet_scales = "free_y",
                                 ncol = 2,
                                 y_axis_label = "Relative abundance (log)",
                                 workdir = "d:/Google Drive/Documents/UDel/Research/Concrete_analysis/"
                                 ){
  
  library(reshape2)
  library(lubridate)
  
  
  otu_table <- file.path(workdir,otu_table) 
  taxonomy_tsv <- file.path(workdir,taxonomy_tsv)
  metadata <- file.path(workdir,metadata)
  
  otuTable <- read_qza(otu_table)$data %>% as.matrix() %>% 
    microbiome::transform(transform = transformation)
  
  tax_ranks <- c("kingdom","phylum","class","order","family","genus","species")
  
  taxonomy <- read_delim(taxonomy_tsv, "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
    select(-confidence) %>% 
    mutate(for_split = taxonomy) %>% 
    separate(for_split,tax_ranks,sep = "; ") %>% 
    mutate_at(tax_ranks,str_remove,pattern = "[a-z]__")
  
  #For default SILVA string parsing
  #separate(for_split,tax_ranks,sep = ";") %>% 
  #mutate_at(tax_ranks,str_remove,pattern = "D_[0-9]__")
  
  met <- read_delim(metadata, "\t", escape_double = FALSE, trim_ws = TRUE)
  
  summary_table <- t(otuTable) %>% 
    reshape2::melt(value.name = "rel_abund") %>% 
    filter(rel_abund > 0) %>%
    dplyr::rename("#SampleID" = Var1, "#OTU ID" = Var2) %>% 
    left_join(taxonomy) %>% 
    left_join(met) %>%
    group_by(`#SampleID`) %>%
    ungroup()
  
  
  #Set taxa as ordered factors, so that they'e plotted in the right order
  taxa <- taxa %>% factor(.,levels = .,ordered = TRUE)
  
  #make table of interesting taxa function
  
  tax_ranks <- c("kingdom","phylum","class","order","family","genus","species","#OTU ID")  
  
  plot_table <- summary_table %>% mutate(otu_name = .$`#OTU ID`) %>%
    filter_at(tax_ranks, any_vars(. %in% taxa & rel_abund >0)) %>%
    select(`#SampleID`,ASR,Months,Collected,Type,Type2,rel_abund,taxonomy,temp_avg,otu_name,one_of(tax_ranks)) %>%
    gather(key="Taxa_type",value = "Taxa",-c(`#SampleID`,ASR,Months,temp_avg,Collected,rel_abund,taxonomy,otu_name,Type,Type2)) %>%
    filter(Taxa %in% taxa) %>% 
    mutate(for_split = taxonomy) %>% 
    separate(for_split,tax_ranks,sep = "; ") %>% 
    mutate_at(tax_ranks,str_remove,pattern = "[a-z]__") %>%
    mutate(og_taxa = Taxa, Taxa = ifelse(Taxa_type == "species", paste(genus, species),Taxa), Taxa = ifelse(Taxa_type == "#OTU ID", paste0(taxonomy,"\n",otu_name),Taxa), "#OTU ID" = otu_name, tax_and_level= paste(Taxa," (",substr(Taxa_type, start = 1, stop = 1),")",sep = ""), s_date = mdy(Collected)) %>%
    select(-otu_name)

  if(agglom == TRUE) plot_table <- plot_table %>% 
                        group_by(Type, Type2, `#SampleID`,tax_and_level,temp_avg) %>% 
                        summarise(rel_abund = sum(rel_abund))
  
  months_between_label = 3
  
  if (split_asr == TRUE){x_group <- quo(Type2)
    } else {x_group <- quo(Type)}
  
  ##Make plots
    #({if(split_asr == TRUE) ggplot(plot_table, aes(Type2, rel_abund, label = `#SampleID`, group = Type2))} +
       #{if(split_asr == FALSE) ggplot(plot_table, aes(Type, rel_abund, label = `#SampleID`, group = Type))} +
      (ggplot(plot_table, aes(!!x_group, rel_abund, label = `#SampleID`, group = !!x_group)) + 
       scale_fill_viridis_c() +
       geom_boxplot() +
       geom_jitter(alpha = point_alpha, size = 2, width = 0.1) + 
       scale_y_log10() +
       theme(legend.key.size = unit(.5,"inches"), legend.title=element_text(size=16) , legend.text=element_text(size=16)) +
       theme(axis.text=element_text(size=16), axis.title=element_text(size=24)) +
       theme(strip.text.x = element_text(size = 24, colour = "black", angle = 0)) +
       labs(x= NULL, y= "Relative abundance (log)", shape = "ASR", color = "ASR", fill = "Temp. (F)", alpha = "Temp. (F)", linetype = "ASR") +
       theme(legend.position = c(0.9, 0.6),legend.background = element_rect(fill = "white", colour = NA)) +
       theme_bw() +
       theme(axis.text.x = element_text(angle=45,vjust = 1, hjust = 1)) +
       facet_wrap(~ tax_and_level,ncol = ncol,scales = facet_scales))
}




#plot prevalence over time

plot_prevalence <- function(taxa,
                      otu_table="d:/Google Drive/Documents/UDel/Research/Concrete_analysis/data/processed/biom_tables/concrete_only/all_w_taxonomy_decontam_table.biom_concrete.biom_table.qza",
                      taxonomy = "d:/Google Drive/Documents/UDel/Research/Concrete_analysis/data/processed/silva_taxonomy.tsv",
                      metadata = "d:/Google Drive/Documents/UDel/Research/Concrete_analysis/data/sample-metadata.tsv",
                      agglom=TRUE,
                      y_0_to_1 = TRUE,
                      smooth_method = "loess",
                      se = TRUE,
                      point_alpha = 4/10,
                      ncol = 2,
                      facet_scales = "free_y", 
                      temp_shade = FALSE, 
                      transformation = "compositional",
                      filter_zeros = FALSE){
  
  summary_table <- process_for_plot(otu_table,taxonomy,metadata,transformation)
  
  #Set taxa as ordered factors, so that they'e plotted in the right order
  taxa <- taxa %>% factor(.,levels = .,ordered = TRUE)
  
  #make table of interesting taxa function
  
  tax_ranks <- c("kingdom","phylum","class","order","family","genus","species","#OTU ID")  
  
  plot_table <- summary_table %>% mutate(otu_name = .$`#OTU ID`) %>%
    filter_at(tax_ranks, any_vars(. %in% taxa)) %>% 
    {if (filter_zeros==TRUE) filter(., rel_abund > 0) else .} %>% 
    select(`#SampleID`,ASR,Months,Collected,rel_abund,taxonomy,temp_avg,otu_name,one_of(tax_ranks)) %>%
    gather(key="Taxa_type",value = "Taxa",-c(`#SampleID`,ASR,Months,temp_avg,Collected,rel_abund,taxonomy,otu_name)) %>%
    filter(Taxa %in% taxa) %>% 
    mutate(for_split = taxonomy) %>% 
    separate(for_split,tax_ranks,sep = "; ") %>% 
    mutate_at(tax_ranks,str_remove,pattern = "[a-z]__") %>%
    mutate(og_taxa = Taxa, Taxa = ifelse(Taxa_type == "species", paste(genus, species),Taxa), Taxa = ifelse(Taxa_type == "#OTU ID", paste0(taxonomy,"\n",otu_name),Taxa), "#OTU ID" = otu_name, tax_and_level= paste(Taxa," (",substr(Taxa_type, start = 1, stop = 1),")",sep = ""), s_date = mdy(Collected)) %>%
    select(-otu_name) %>% 
    group_by(Months)
  
  temp_data <- plot_table %>% select(samp_day = s_date,temp_avg) %>% ungroup() %>% unique()
  
  if(agglom == TRUE) plot_table <- plot_table %>% 
    group_by(`#SampleID`,tax_and_level,s_date,temp_avg,Months) %>% 
    summarise(rel_abund = sum(rel_abund)) %>% 
    ungroup() %>% group_by(Months,tax_and_level)
  
  months_between_label = 3
  
  #Make plots
      plot <- plot_table %>% 
        mutate(prevalence = sum(rel_abund > 0) / length(unique(`#SampleID`))) %>% 
        group_by(s_date,tax_and_level,temp_avg) %>% 
        select(s_date,tax_and_level,temp_avg, prevalence) %>%  distinct() %>% 
      ggplot(aes(s_date,prevalence)) +
      scale_fill_viridis_c() +
      geom_point(alpha = point_alpha, size = 2) + 
      geom_smooth(method = smooth_method, se = se) +
      #facet_wrap(facets = ~ `#OTU ID` + taxonomy, ncol = 3,scales = "free_y", labeller = label_context) + 
      theme(legend.key.size = unit(.5,"inches"), legend.title=element_text(size=16) , legend.text=element_text(size=16)) +
      theme(axis.text=element_text(size=16), axis.title=element_text(size=24)) +
      theme(strip.text.x = element_text(size = 24, colour = "black", angle = 0)) +
      labs(x= NULL, y= "Prevalence") +
      theme(legend.position = c(0.9, 0.6),legend.background = element_rect(fill = "white", colour = NA)) +
      theme_bw() +
      scale_x_date(labels=date_format("%b '%y"),date_breaks = paste(months_between_label,"months"),limits = c(min(mdy(summary_table$Collected)),max(mdy(summary_table$Collected)))) +
      theme(axis.text.x = element_text(angle=45,vjust = 1, hjust = 1)) +
      {if(y_0_to_1 == TRUE) ylim(0,1)} +
      facet_wrap(~ tax_and_level,ncol = ncol,scales = facet_scales)
    print(plot)
  }




plot_emp_distribution <- function(taxa = NA,
                                  otu_of_interest = NA,
                                  scale_x_log10 = FALSE,
                                  prevalence_or_rel_freq = "rel_freq",
                                  points_for_no_ridge = FALSE,
                                  only_concrete_asvs = TRUE,
                                  workdir = "d:/Google Drive/Documents/UDel/Research/Concrete_analysis/"){
  setwd(workdir)
  library(ggridges)
  library(qiime2R)
  library(tidyverse)
  
  emp_tapwater_concrete_map <- read_delim("data/processed/emp_data/emp_tapwater_concrete_map.tsv", "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
    dplyr::rename(sample = "#SampleID") %>% filter(include_in_concrete_comparison == TRUE)
  
  #Colors for plotting
  empo3 <- levels(as.factor(emp_tapwater_concrete_map$empo_3)) %>% 
    data.frame(empo_3 = .) %>% 
    cbind( colors = c("grey80", "#ff0000", "#f27304", "#f27304", "#ff0000", "#ff0000", "black", "#fcc688", "maroon2", "maroon2", "#52ff3c", "#500261", "#80c99b", "#045b04", "#80c99b", "maroon2", "maroon2",  "#ff0000", "#a36f23", "#765f3d", "#6b440b", "#ffff00", "#ababab", "#5e5e5e", "#7cecf4", "#0000ff","#000098")) %>%
    mutate(empo_3 = as.character(empo_3), colors = as.character(colors)) #%>%
  #filter(!empo_3 %in% series_samples)
  
  taxonomy <- read_qza("data/processed/emp_data/emp_tapwater_and_concrete_silva_taxonomy.qza")$data %>% 
    select(otu = "Feature.ID", taxonomy = "Taxon")
  
  #tax_levels <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  
  rel_freq <- read_rds("data/processed/emp_data/emp_rel_freq_w_tax.rds") 
  
  concrete_otus <- read_rds("data/processed/emp_data/emp_otus_in_concrete.rds")
  
  tax_of_interest <- taxa
  
  otus_of_interest <- taxonomy %>% 
    do(if(!is.na(taxa)) filter(.,str_detect(taxonomy,tax_of_interest)) else .) %>%
    do(if(!is.na(otu_of_interest)) filter(.,str_detect(otu,otu_of_interest)) else .) %>% 
    do(if(only_concrete_asvs == FALSE) . else filter(., otu %in% concrete_otus)) %>% #toggle whether to only look at concrete OTUs or all, concrete_otus generated above
    pull(otu) %>% as.character()
  
  filt_freq <- rel_freq %>%
    filter(otu %in% otus_of_interest) %>%
    group_by(empo_3) %>% 
    mutate(mean_freq = mean(rel_freq)) %>%
    left_join(empo3) %>%
    ungroup() %>% mutate(empo_3 = as.character(empo_3)) 
  
  order <- filt_freq %>% select(empo_3,mean_freq) %>% distinct()
  
  filt_freq %>%
    ggplot(aes(y = reorder(empo_3,mean_freq),!!sym(prevalence_or_rel_freq), fill = colors, color = colors)) +
    {if(points_for_no_ridge == FALSE) geom_density_ridges(alpha = 0.5,jittered_points = TRUE)} +
    {if(points_for_no_ridge == TRUE) geom_density_ridges(alpha = 0.5,jittered_points = FALSE)} +
    {if(points_for_no_ridge == TRUE) geom_point()} +
    #geom_violin(alpha = 0.4) +
    lims(x = c(0,1)) +
    {if(scale_x_log10 == TRUE) scale_x_log10()} +
    scale_fill_identity() +
    scale_color_identity() +
    labs(y = NULL, x = "Relative frequency of observation",title = taxa) +
    {if(prevalence_or_rel_freq == "prevalence") labs(x = "Prevalence")} +
    theme_minimal()
  
}