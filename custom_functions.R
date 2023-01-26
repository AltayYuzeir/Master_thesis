top_N_featuers_per_cluster = function(FindAllMarkers_object, N = 20){
  require(tidyverse)
  clusters = as.character(unique(FindAllMarkers_object$cluster))
  all_features = data.frame(matrix(nrow = 1, ncol = N))
  
  for(i in (1:length(clusters)-1)){
    data = subset(FindAllMarkers_object, cluster == i)
    features = arrange(data, desc(avg_log2FC))
    features = features[1:N,]
    all_features[i+1,] = rownames(features)
    
  }
  
  all_features = as.data.frame(t(all_features))
  
  all_features = as.data.frame(lapply(all_features, function(y) gsub("NA.*", "---", y)))
  all_features = as.data.frame(lapply(all_features, function(y) gsub("\\..*", "", y)))
  cols = character()
  for(i in (1:length(clusters)-1)) cols[i+1] = paste0("Cluster_", i)
  colnames(all_features) = cols
  
  rows = character()
  for(i in 1:N ) rows[i] = paste0("Feature_", i)
  rownames(all_features) = rows
  
  return(all_features)
}


GSEA_calc_gene <- function(gene_list, # provided by Dr. Jiangyan Yu
                           DEG_list,
                           comparison = NULL,
                           species = "M", # H for human
                           genes_down = NULL,
                           genes_all = NULL,
                           GeneSets =c("GO","KEGG","DO","Hallmark","cannonicalPathways","Motifs","ImmunoSignatures"),
                           GOntology = "BP", #Alternative "MF" or "CC"
                           pCorrection = "bonferroni", # choose the p-value adjustment method
                           pvalueCutoff = 0.1, # set the unadj. or adj. p-value cutoff (depending on correction method)
                           qvalueCutoff = 0.2 # set the q-value cutoff (FDR corrected)
) {
  
  require(biomaRt)
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  require(clusterProfiler)
  require(DOSE)
  
  gsea_pathway_dir = ""
  
  # define universe
  universe <- as.character(gene_list)
  # change symbols to ENTREZ IDs (necessary for ClusterProfiler)
  if(species == "M") {
    universe = getLDS(attributes = c("mgi_symbol"), 
                      filters = "mgi_symbol", 
                      values = universe, 
                      mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/"), 
                      attributesL = c("hgnc_symbol"), 
                      martL = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/"), 
                      uniqueRows=T)[,2]
  }else{}
  
  universe_Entrez <- bitr(universe, 
                          fromType="SYMBOL", 
                          toType="ENTREZID", 
                          OrgDb="org.Hs.eg.db")$ENTREZID
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  # human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://www.ensembl.org")
  # mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://www.ensembl.org")
  
  
  genes_up = DEG_list 
  # genes_down = DEG_down_list 
  # genes_all = DEG_all_list
  
  if(species == "M") {
    genes_up <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genes_up, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)[,2]
    # genes_down <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genes_down, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)[,2]
    # genes_all <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genes_all, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)[,2]
  }
  
  cannonicalPathway_genes <- clusterProfiler::read.gmt(paste(gsea_pathway_dir,"c2.all.v7.5.entrez.gmt",sep=""))
  immuno_genes <- clusterProfiler::read.gmt(paste(gsea_pathway_dir,"c7.all.v7.5.entrez.gmt",sep=""))
  hallmark_genes <- clusterProfiler::read.gmt(paste(gsea_pathway_dir,"h.all.v7.5.entrez.gmt",sep=""))
  motifs <- clusterProfiler::read.gmt(paste(gsea_pathway_dir,"c3.all.v7.5.entrez.gmt",sep=""))
  
  # if(!is.null(group_by)){
  # SetIdents(object = tmp) <- group_by
  # } 
  
  # if(!is.null(genes_down)) {
  #   top_down <- genes_down
  #   entrez_down <- bitr(top_down, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  # }
  
  # if(!is.null(genes_all)) {
  #   top_all <- genes_all
  #   entrez_all <- bitr(top_all, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  # }
  
  top_up <- genes_up
  entrez_up <- bitr(top_up, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  
  OrgDb = org.Hs.eg.db
  
  results <- list()
  
  
  # GO enrichment
  if("GO" %in% GeneSets){
    print("Performing GO enrichment")
    results$GOup <- as.data.frame(enrichGO(gene = entrez_up,
                                           universe = universe_Entrez,
                                           OrgDb = OrgDb,
                                           ont = GOntology,
                                           pAdjustMethod = pCorrection,
                                           pvalueCutoff  = pvalueCutoff,
                                           qvalueCutoff  = qvalueCutoff,
                                           readable      = T))
    
    if(nrow(results$GOup)>0){results$GOup$Enrichment <- paste("GO enrichment for genes upregulated in comparison ",comparison,sep="")}
    
    # if(!is.null(genes_down)) {
    #   results$GOdown <- as.data.frame(enrichGO(gene = entrez_down,
    #                                            universe = universe_Entrez,
    #                                            OrgDb = OrgDb,
    #                                            ont = GOntology,
    #                                            pAdjustMethod = pCorrection,
    #                                            pvalueCutoff  = pvalueCutoff,
    #                                            qvalueCutoff  = qvalueCutoff,
    #                                            readable      = T))
    #   if(nrow(results$GOdown)>0){results$GOdown$Enrichment <- paste("GO enrichment for genes downregulated in comparison ",comparison,sep="")}
    # }
    
    #   if(!is.null(genes_all)) {
    #     results$GOall <- as.data.frame(enrichGO(gene = entrez_all,
    #                                             universe = universe_Entrez,
    #                                             OrgDb = OrgDb,
    #                                             ont = GOntology,
    #                                             pAdjustMethod = pCorrection,
    #                                             pvalueCutoff  = pvalueCutoff,
    #                                             qvalueCutoff  = qvalueCutoff,
    #                                             readable      = T))
    #     if(nrow(results$GOall)>0){results$GOall$Enrichment <- paste("GO enrichment for genes allregulated in comparison ",comparison,sep="")}
    #   }
  }
  
  # KEGG enrichment
  if("KEGG" %in% GeneSets){
    print("Performing KEGG enrichment")
    
    org = "hsa"
    
    results$KEGGup <- as.data.frame(enrichKEGG(gene = entrez_up, 
                                               organism = org,
                                               universe = universe_Entrez, 
                                               pAdjustMethod = pCorrection,
                                               pvalueCutoff  = pvalueCutoff,
                                               qvalueCutoff = qvalueCutoff))
    if(nrow(results$KEGGup)>0){results$KEGGup$Enrichment <- paste("KEGG enrichment for genes upregulated in comparison ",comparison,sep="")}
    
    # if(!is.null(genes_down)) {
    #   results$KEGGdown <- as.data.frame(enrichKEGG(gene = entrez_down, 
    #                                                organism = org,
    #                                                universe = universe_Entrez, 
    #                                                pAdjustMethod = pCorrection,
    #                                                pvalueCutoff  = pvalueCutoff,
    #                                                qvalueCutoff = qvalueCutoff))
    #   if(nrow(results$KEGGdown)>0){results$KEGGdown$Enrichment <- paste("KEGG enrichment for genes downregulated in comparison ",comparison,sep="")}
    # }
  }
  
  # DO enrichment
  if("DO" %in% GeneSets){
    print("Performing Disease Ontology enrichment")
    
    results$DOup <- as.data.frame(enrichDO(gene = entrez_up, 
                                           universe = universe_Entrez, 
                                           pAdjustMethod = pCorrection,
                                           pvalueCutoff  = pvalueCutoff,
                                           qvalueCutoff = qvalueCutoff,
                                           minGSSize     = 5,
                                           maxGSSize     = 500,
                                           readable=TRUE))
    if(nrow(results$DOup)>0){results$DOup$Enrichment <- paste("DO enrichment for genes upregulated in comparison ",comparison,sep="")}
    
    # if(!is.null(genes_down)) {
    #   results$DOdown <- as.data.frame(enrichDO(gene = entrez_down, 
    #                                            universe = universe_Entrez, 
    #                                            pAdjustMethod = pCorrection,
    #                                            pvalueCutoff  = pvalueCutoff,
    #                                            qvalueCutoff = qvalueCutoff,
    #                                            minGSSize     = 5,
    #                                            maxGSSize     = 500,
    #                                            readable=TRUE))
    #   if(nrow(results$DOdown)>0){results$DOdown$Enrichment <- paste("DO enrichment for genes downregulated in comparison ",comparison,sep="")}
    # }
  }
  
  # Hallmark enrichment
  if("Hallmark" %in% GeneSets){
    print("Performing Hallmark enrichment")
    
    results$HALLMARKup <- as.data.frame(enricher(entrez_up,
                                                 TERM2GENE=hallmark_genes,
                                                 universe = universe_Entrez,  
                                                 pAdjustMethod = pCorrection,
                                                 pvalueCutoff  = pvalueCutoff,
                                                 qvalueCutoff = qvalueCutoff))
    if(nrow(results$HALLMARKup)>0){results$HALLMARKup$Enrichment <- paste("HALLMARK enrichment for genes upregulated in comparison ",comparison,sep="")}
    
    # if(!is.null(genes_down)) {
    #   results$HALLMARKdown <- as.data.frame(enricher(entrez_down,
    #                                                  TERM2GENE=hallmark_genes,
    #                                                  universe = universe_Entrez,  
    #                                                  pAdjustMethod = pCorrection,
    #                                                  pvalueCutoff  = pvalueCutoff,
    #                                                  qvalueCutoff = qvalueCutoff))
    #   if(nrow(results$HALLMARKdown)>0){results$HALLMARKdown$Enrichment <- paste("HALLMARK enrichment for genes downregulated in comparison ",comparison,sep="")}
    # }
  }
  
  # Cannonical Pathway enrichment
  if("cannonicalPathways" %in% GeneSets){
    print("Performing Cannonical Pathway (C2) enrichment")
    
    results$cannonicalPathwaysup <- as.data.frame(enricher(entrez_up,
                                                           TERM2GENE=cannonicalPathway_genes,
                                                           universe = universe_Entrez,  
                                                           pAdjustMethod = pCorrection,
                                                           pvalueCutoff  = pvalueCutoff,
                                                           qvalueCutoff = qvalueCutoff))
    if(nrow(results$cannonicalPathwaysup)>0){results$cannonicalPathwaysup$Enrichment <- paste("Cannonical pathway enrichment for genes upregulated in comparison ",comparison,sep="")}
    
    # if(!is.null(genes_down)) {
    #   results$cannonicalPathwaysdown <- as.data.frame(enricher(entrez_down,
    #                                                            TERM2GENE=cannonicalPathway_genes,
    #                                                            universe = universe_Entrez,  
    #                                                            pAdjustMethod = pCorrection,
    #                                                            pvalueCutoff  = pvalueCutoff,
    #                                                            qvalueCutoff = qvalueCutoff))
    #   if(nrow(results$cannonicalPathwaysdown)>0){results$cannonicalPathwaysdown$Enrichment <- paste("Cannonical pathway enrichment for genes downregulated in comparison ",comparison,sep="")}
    # }
  }
  
  # Motif enrichment
  if("Motifs" %in% GeneSets){
    print("Performing Motif enrichment")
    
    results$Motifup <- as.data.frame(enricher(entrez_up,
                                              TERM2GENE=motifs,
                                              universe = universe_Entrez,  
                                              pAdjustMethod = pCorrection,
                                              pvalueCutoff  = pvalueCutoff,
                                              qvalueCutoff = qvalueCutoff))
    if(nrow(results$Motifup)>0){results$Motifup$Enrichment <- paste("TF binding motif enrichment for genes upregulated in comparison ",comparison,sep="")}
    
    # if(!is.null(genes_down)) {
    #   results$Motifdown <- as.data.frame(enricher(entrez_down,
    #                                               TERM2GENE=motifs,
    #                                               universe = universe_Entrez,  
    #                                               pAdjustMethod = pCorrection,
    #                                               pvalueCutoff  = pvalueCutoff,
    #                                               qvalueCutoff = qvalueCutoff))
    #   if(nrow(results$Motifdown)>0){results$Motifdown$Enrichment <- paste("TF binding motif enrichment for genes downregulated in comparison",comparison,sep="")}
    # }
  }
  
  # Immunosignatures enrichment
  if("ImmunoSignatures" %in% GeneSets){
    print("Performing immunesignature enrichment")
    
    results$ImmSigup <- as.data.frame(enricher(entrez_up,
                                               TERM2GENE=immuno_genes,
                                               universe = universe_Entrez,  
                                               pAdjustMethod = pCorrection,
                                               pvalueCutoff  = pvalueCutoff,
                                               qvalueCutoff = qvalueCutoff))
    if(nrow(results$ImmSigup)>0){results$ImmSigup$Enrichment <- paste("Immunosignature enrichment for genes upregulated in comparison ",comparison,sep="")}
    
    # if(!is.null(genes_down)) {
    #   results$ImmSigdown <- as.data.frame(enricher(entrez_down,
    #                                                TERM2GENE=immuno_genes,
    #                                                universe = universe_Entrez,  
    #                                                pAdjustMethod = pCorrection,
    #                                                pvalueCutoff  = pvalueCutoff,
    #                                                qvalueCutoff = qvalueCutoff))
    #   if(nrow(results$ImmSigdown)>0){results$ImmSigdown$Enrichment <- paste("Immunosignature enrichment for genes downregulated in comparison ",comparison,sep="")}
    # }
  }
  return(results)
}



dotplotGSEA <- function(x,   # provided by Dr. Jiangyan Yu
                        show=25,
                        font.size=10,
                        title.size=10,
                        title.width=100,
                        order="count"){
  if(nrow(x)<1){
    print("No enrichment found.")
  }else{
    x <- if(nrow(x)>show){x[c(1:show),]}else{x}
    if(order=="padj"){
      x <- x[order(x$Count,decreasing=FALSE),] 
      x <- x[order(x$p.adjust,decreasing=TRUE),]
      x$Description <- ifelse(nchar(x$Description)<500,
                              paste(substr(x$Description, 1, 500),"[...]",sep=""),
                              x$Description)
      x$Description <- factor(x$Description, levels = unique(x$Description))
    }
    if(order=="count"){
      x <- x[order(x$Count,decreasing=FALSE),]
      x$Description <- ifelse(nchar(x$Description)>500,
                              paste(substr(x$Description, 1, 500),"[...]",sep=""),
                              x$Description)
      x$Description <- factor(x$Description, levels = unique(x$Description))
      x$GeneRatio <- factor(x$GeneRatio, levels = unique(x$GeneRatio))
    }
    
    ggplot(x, aes(x = GeneRatio, y = Description, color = p.adjust)) +
      geom_point(aes(size = Count)) +
      scale_colour_gradientn(colours=c('red', 
                                       'orange', 
                                       'darkblue',
                                       'darkblue'),
                             limits=c(0,1),
                             values   = c(0,0.05,0.2,0.5,1),
                             breaks   = c(0.05,0.2,1),
                             labels = format(c(0.05,0.2,1))) +
      ylab(NULL) +
      ggtitle(paste(strwrap(unique(x$Enrichment), width=title.width), collapse = "\n"))+
      theme_classic() +
      scale_y_discrete(position = "right") +
      theme(text = element_text(size=font.size),
            plot.title = element_text(size=title.size)) 
  }
}


nichenetr_all_genes = function(seurat.object,        # provided by Dr. Jiangyan Yu
                               receiver.cells, 
                               sender.cells,
                               DE_table_receiver, 
                               label.receiver.cells, 
                               label.sender.cells, 
                               top_n_target_per_ligand, 
                               geneset_oi_avg_log2FC, 
                               best_ligand_top_n){
  require(nichenetr)
  function_dir = ""
  
  RNA.counts=seurat.object@assays$RNA@counts
  receiver.cells=receiver.cells
  sender.cells=sender.cells
  DE_table_receiver=DE_table_receiver
  label.receiver.cells=label.receiver.cells
  label.sender.cells=label.sender.cells
  top_n_target_per_ligand=top_n_target_per_ligand
  geneset_oi_avg_log2FC=geneset_oi_avg_log2FC
  best_ligand_top_n=best_ligand_top_n
  
  ### import nichenet database
  # Nichenetr packages
  ligand_target_matrix = readRDS(paste(function_dir,"nichenetr_ligand_target_matrix.rds", sep = ""))
  
  ##convert from human to mouse
  colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
  rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
  ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
  
  lr_network = readRDS(paste(function_dir,"nichenetr_lr_network.rds",sep = ""))
  
  lr_network$from_mouse <- lr_network$from %>% convert_human_to_mouse_symbols()
  lr_network$to_mouse <- lr_network$to %>% convert_human_to_mouse_symbols()
  
  weighted_networks = readRDS(paste(function_dir,"nichenetr_weighted_networks.rds",sep = ""))
  weighted_networks_m = weighted_networks
  weighted_networks_m$lr_sig = weighted_networks_m$lr_sig %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% tidyr::drop_na()
  weighted_networks_m$gr = weighted_networks_m$gr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% tidyr::drop_na()
  
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
  
  weighted_networks_lr = weighted_networks_lr %>% mutate(from_mouse = convert_human_to_mouse_symbols(from), to_mouse = convert_human_to_mouse_symbols(to)) %>% tidyr::drop_na()
  
  ligand_tf_matrix = readRDS(paste(function_dir,"nichenetr_ligand_tf_matrix.rds", sep = ""))
  colnames(ligand_tf_matrix) = ligand_tf_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
  rownames(ligand_tf_matrix) = ligand_tf_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
  ligand_tf_matrix = ligand_tf_matrix %>% .[!is.na(rownames(ligand_tf_matrix)), !is.na(colnames(ligand_tf_matrix))]
  
  sig_network = readRDS(paste(function_dir,"nichenetr_signaling_network.rds", sep = ""))
  sig_network$from <- sig_network$from %>% convert_human_to_mouse_symbols()
  sig_network$to <- sig_network$to %>% convert_human_to_mouse_symbols()
  
  gr_network = readRDS(paste(function_dir,"nichenetr_gr_network.rds", sep = ""))
  gr_network$from <- gr_network$from %>% convert_human_to_mouse_symbols()
  gr_network$to <- gr_network$to %>% convert_human_to_mouse_symbols()
  
  ### define genes of interests
  expressed_genes_receiver = RNA.counts %>% .[,receiver.cells] %>% apply(1,function(x){sum(x>0)/length(x)}) %>% .[. >=0.1] %>% names()
  
  receptors = lr_network %>% pull(to_mouse) %>% unique()
  expressed_receptors = intersect(receptors, expressed_genes_receiver)
  
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  ## define gene of interest
  # DE_table_receiver = FindMarkers(object = seurat.object, ident.1 = receiver.cells, min.pct = 0.10) %>% rownames_to_column("gene")
  
  geneset_oi = DE_table_receiver %>% dplyr::filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= geneset_oi_avg_log2FC) %>% pull(gene) %>% .[. %in% rownames(ligand_target_matrix)] %>% as.character()
  
  ligands = lr_network %>% pull(from_mouse) %>% unique()
  
  sender_genes_sender = RNA.counts %>% .[,sender.cells] %>% apply(1,function(x){sum(x>0)/length(x)}) %>% .[. >=0.1] %>% names()
  expressed_ligands_sender = intersect(ligands,sender_genes_sender)
  expressed_ligands = expressed_ligands_sender
  
  potential_ligands = lr_network %>% dplyr::filter(from_mouse %in% expressed_ligands & to_mouse %in% expressed_receptors) %>% pull(from_mouse) %>% unique()
  
  ligand_activities = predict_ligand_activities(
    geneset = geneset_oi, 
    background_expressed_genes = background_expressed_genes, 
    ligand_target_matrix = ligand_target_matrix, 
    potential_ligands = potential_ligands)
  
  ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(plyr::desc(pearson)))
  ligand_activities
  
  best_upstream_ligands = ligand_activities %>% 
    top_n(best_ligand_top_n, pearson) %>%
    arrange(-pearson) %>% 
    pull(test_ligand) %>%
    unique()
  
  # DotPlot(seurat.object, features = best_upstream_ligands %>% rev(), cols = "Blue", split.by = "celltype") + RotatedAxis()
  
  # show histogram of ligand activity scores
  ## define cut-off for the number of selected ligands
  p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
    geom_histogram(color="black", fill="darkorange")  + 
    # geom_density(alpha=.1, fill="orange") +
    geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
    labs(x="ligand activity (PCC)", y = "# ligands") +
    theme_classic()
  
  active_ligand_target_links_df = best_upstream_ligands %>% 
    lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = top_n_target_per_ligand) %>% 
    bind_rows() %>% drop_na()
  
  ### for visualisation
  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.1) ## old cutoff was 0.25, new is 0.33
  
  order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% 
    rev() %>% make.names()
  order_targets = active_ligand_target_links_df$target %>% 
    unique() %>% 
    intersect(rownames(active_ligand_target_links)) %>% 
    make.names()
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  
  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% 
    t()
  
  p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot(paste("Prioritized",label.sender.cells, "ligands"), paste("DE genes in", label.receiver.cells, "cells"), color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))
  
  lr_network_top = lr_network %>% 
    dplyr::filter(from_mouse %in% best_upstream_ligands & to_mouse %in% expressed_receptors) %>% 
    distinct(from_mouse,to_mouse)
  
  best_upstream_receptors = lr_network_top %>% pull(to_mouse) %>% unique()
  
  lr_network_top_df_large = weighted_networks_lr %>% 
    dplyr::filter(from_mouse %in% best_upstream_ligands & to_mouse %in% best_upstream_receptors)
  
  lr_network_top_df = lr_network_top_df_large %>% tidyr::spread("from_mouse","weight",fill = 0)
  ###remove from and to columns due to mouse orgnism
  lr_network_top_df = lr_network_top_df[,c(-1,-2)]
  lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to_mouse) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to_mouse)
  
  # perform hierarchical clustering to order the ligands and receptors
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  dist_ligands[is.na(dist_ligands)] = 0
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
  
  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
  
  p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
  
  print("Please mannually check the ligand_activities matrix to explore which ligands are biologically meaningful. Then set best_ligand_top_n accordingly.")
  
  return(list(p_hist_lig_activity=p_hist_lig_activity,
              p_ligand_target_network=p_ligand_target_network,
              p_ligand_receptor_network=p_ligand_receptor_network,
              order_receptors= order_receptors,
              order_ligands_receptor=order_ligands_receptor,
              ligand_activities = ligand_activities,
              active_ligand_target_links_df = active_ligand_target_links_df,
              lr_network_top_df_large = lr_network_top_df_large %>% dplyr::rename(ligand = from_mouse, receptor = to_mouse)))  
}

custom_LIANA_dotplot = function (liana_res, 
                                 source_groups = NULL, 
                                 target_groups = NULL,
                                 top_interactions = 25, 
                                 significance_threshold = Inf,
                                 scale_color.viridis = "magma",
                                 specificity = "natmi.edge_specificity", 
                                 magnitude = "sca.LRscore",
                                 y.label = "Interactions (Ligand => Receptor)", 
                                 size.label = "Interaction Specificity",
                                 colour.label = "Expression Magnitude",
                                 sort_for = c("magnitude", "specificity"),
                                 target_colors = RColorBrewer::brewer.pal(8, "Dark2"),
                                 source_color = "gray6",
                                 show_complex = TRUE,
                                 scale_color.inverse = FALSE,
                                 size_range = c(2, 10),
                                 compute_ChordPlots = FALSE,
                                 highlight.groups = NULL, 
                                 source_text.size = 12,
                                 target_text.size = 12,
                                 target_text.angle = 0) 
{  require(dplyr)
  require(ggplot2)
  require(tidyr)
  require(forcats)
  require(magrittr)
  
  require(RColorBrewer)
  
  if(scale_color.inverse == FALSE) scale_color.inverse = 1
  else scale_color.inverse = -1
  
  if (target_text.angle == 0) {
    target_text.hjust = NULL
    target_text.vjust = NULL
  } else if (target_text.angle == 45) {
    target_text.hjust = 1
    target_text.vjust = NULL
  } else {
    target_text.hjust = 1
    target_text.vjust = 0.5
  }
  
  if (show_complex) entities <- c("ligand.complex", "receptor.complex")
  else entities <- c("ligand", "receptor")
  
  if (!is.null(source_groups) & !is.null(target_groups)) {
    liana_mod = filter(liana_res, source %in% source_groups)
    liana_mod = filter(liana_mod, target %in% target_groups)
  } else if (!is.null(source_groups) & is.null(target_groups)) {
    liana_mod = filter(liana_res, source %in% source_groups)
  } else if (is.null(source_groups) & !is.null(target_groups)) {
    liana_mod = filter(liana_res, target %in% target_groups)
  } else
    liana_mod = liana_res
  
  liana_mod %<>% dplyr::rename(magnitude = !!magnitude) %>% dplyr::rename(specificity = !!specificity) %>%
    unite(all_of(entities), col = "interaction", sep = " => ", remove = FALSE) %>%
    unite(c("source", "target"), col = "source_target", remove = FALSE)
  
  if(sort_for == "magnitude") 
    liana_mod = arrange(liana_mod,dplyr::desc(magnitude),dplyr::desc(specificity))
  else liana_mod = arrange(liana_mod, dplyr::desc(specificity),dplyr::desc(magnitude))
  
  liana_mod = filter(liana_mod, aggregate_rank <= significance_threshold)
  
  if (!is.null(top_interactions)) {
    top_int <- liana_mod %>% distinct_at(entities) %>% head(top_interactions)
    liana_mod %<>% inner_join(top_int, by = entities)
  }
  
  liana_copy = liana_mod
  
  p = suppressWarnings(
    ggplot( liana_mod,
            aes(
              x = target,
              y = reorder(interaction,dplyr::desc(interaction)),
              colour = magnitude,
              size = specificity,
              group = target
            )
    ) +
      geom_point() + viridis::scale_color_viridis(option = scale_color.viridis,
                                                  direction = scale_color.inverse) +
      labs( label = NULL, subtitle = "Source", caption = "Target") +
      scale_size_continuous(range = size_range) + 
      facet_grid( . ~ source, space = "free", scales = "free", switch = "y") +
      theme_bw(base_size = 16) +
      theme(
        axis.text.x = element_text(
          color = target_colors
          # c(
          # "#E69F00",
          # "#56B4E9",
          # "#009E73",
          # "#F0E442",
          # "#0072B2",
          # "#D55E00",
          # "#CC79A7" )
          ,
          face = "bold",
          size = target_text.size,
          angle = target_text.angle,
          hjust = target_text.hjust,
          vjust = target_text.vjust
        ),
        axis.text.y = element_text(
          size = 12,
          vjust = 0.5,
          face = "bold",
          colour = "gray6"
        ),
        axis.title.y = element_text(size = 12, face = "bold"),
        
        strip.text.x = element_text(size = source_text.size, 
                                    color = source_color, 
                                    face = "bold"),
        legend.box.spacing = unit(1.5, "cm"),
        legend.position = "bottom",
        legend.text.align = 0,
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.box = "horizontal",
        legend.key.width = unit(2, "cm"),
        legend.key.height =  unit(0.6, "cm"),
        plot.subtitle = element_text(
          vjust = 0.5,
          hjust = 0.5,
          size = 12,
          face = "bold"
        ),
        plot.caption = element_text(
          vjust = 45,
          hjust = 0.5,
          size = 12,
          face = "bold"
        ),
        panel.grid = element_line(color = "grey90", linetype = "dotted"),
        
        panel.spacing = unit(0.1, "lines"),
        strip.background = element_rect(fill = NA)
      ) +
      guides(
        colour = guide_colourbar(title.position = "top", title.hjust = 0.5),
        size = guide_legend(title.position = "top", title.hjust = 0.5)
      ) +
      labs(
        y = y.label,
        colour = colour.label,
        size = size.label,
        x = NULL
      )
  )
  
  if (isTRUE(compute_ChordPlots)) {
    data <- liana_copy %>% dplyr::select(dplyr::all_of(c("source",
                                                         "target"))) %>% dplyr::group_by(.data$target, .data$source) %>%
      dplyr::summarise(value = dplyr::n()) %>% dplyr::rename(from = .data[["source"]],
                                                             to = .data[["target"]]) %>%
      dplyr::select(dplyr::all_of(c("from", "to", "value")))
    
    
    p.source_target <- SCpubr::do_ChordDiagramPlot(
      from_df = TRUE,
      df = data,
      padding_labels = 0,
      self.link = 2,
      highlight_group = highlight.groups,
      annotationTrack = c("grid"),
      link.border.color = NA,
      z_index = TRUE
    )
    data <-
      liana_copy %>% dplyr::select(dplyr::all_of(c("ligand.complex",
                                                   "receptor.complex"))) %>%
      dplyr::group_by(.data$ligand.complex, .data$receptor.complex) %>% dplyr::summarise(value = dplyr::n()) %>%
      dplyr::rename(from = .data[["ligand.complex"]],
                    to = .data[["receptor.complex"]]) %>%
      dplyr::select(dplyr::all_of(c("from", "to", "value")))
    
    ligands = unique(data$from)
    receptors = unique(data$to)
    # set.seed(12345)
    #colors = distinctive_colors(length(ligands), length(receptors))
    
    #colors_from = setNames(colors$Palette1, ligands)
    #colors_to = setNames(colors$Palette2, receptors)
    
    p.ligand_receptor <-
      SCpubr::do_ChordDiagramPlot(
        from_df = TRUE,
        #group = c("source", "target"),
        df = data,
        #colors.from = colors_from,
        #colors.to = colors_to,
        padding_labels = 0,
        self.link = 2,
        highlight_group = highlight.groups,
        annotationTrack = c("grid"),
        link.border.color = NA,
        z_index = TRUE
      )
    return(
      list(
        dotplot = p,
        chord_total_interactions = p.source_target,
        chord_ligand_receptor = p.ligand_receptor
      )
    )
    
  } else
    return(p)
}

custom_LIANA_top_interactions = function (liana_res, 
                                          top_N_interactions = 30, 
                                          source_groups = NULL, 
                                          target_groups = NULL, 
                                          specificity = "natmi.edge_specificity", 
                                          magnitude = "sca.LRscore",
                                          significance_threshold = Inf,
                                          show_complex = TRUE,
                                          sort_for = c("magnitude", "specificity")
)
  #
  
{ require(dplyr)
  require(tidyr)
  require(forcats)
  require(magrittr)
  
  if (show_complex) {
    entities <- c("ligand.complex", "receptor.complex")
  }
  else {
    entities <- c("ligand", "receptor")
  }
  
  if (!is.null(source_groups) & !is.null(target_groups)) {
    liana_mod = filter(liana_res, source %in% source_groups)
    liana_mod = filter(liana_mod, target %in% target_groups)
  } else if (!is.null(source_groups) & is.null(target_groups)) {
    liana_mod = filter(liana_res, source %in% source_groups)
  } else if (is.null(source_groups) & !is.null(target_groups)) {
    liana_mod = filter(liana_res, target %in% target_groups)
  } else
    liana_mod = liana_res
  
  
  liana_mod %<>% dplyr::rename(magnitude = !!magnitude) %>% 
    dplyr::rename(specificity = !!specificity) %>%
    unite(all_of(entities),
          col = "interaction",
          sep = " => ",
          remove = FALSE) %>%
    unite(c("source", "target"), col = "source_target",
          remove = FALSE)
  
  if(sort_for == "magnitude") 
    liana_mod = arrange(liana_mod,dplyr::desc(magnitude),dplyr::desc(specificity))
  else liana_mod = arrange(liana_mod, dplyr::desc(specificity),dplyr::desc(magnitude))
  
  liana_mod = filter(liana_mod, aggregate_rank <= significance_threshold)
  
  if (!is.null(top_N_interactions)) {
    top_int <-
      liana_mod %>% distinct_at(entities) %>% head(top_N_interactions)
    liana_mod %<>% inner_join(top_int, by = entities)
  }
  
  liana_mod %>%
    unite(all_of(entities), col = "interaction", sep = " => ", 
          remove = FALSE) %>%
    unite(c("source", "target"), col = "source_target",
          remove = FALSE) -> liana_mod
  
  return(liana_mod)
}
