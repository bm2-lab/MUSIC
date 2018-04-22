#################################       Data Preprocessing     ############################
# ********************** data importing
# if data format is 10X, for convenience, this function is better.
Input_preprocess_10X<-function(directory, sample_info_batch=NULL){
  require(Seurat)
  require(hash)
  require(stringr)
  perturb_seq <- Read10X(directory)
  perturb_seq <- as.matrix(perturb_seq)
  cbc_gbc <- read.table(paste(directory, "cbc_gbc_dict.tsv",sep = "/"), stringsAsFactors = FALSE)
  cbc_gbc <- unique(cbc_gbc)
  cbc_gbc_grna <- read.table(paste(directory, "cbc_gbc_dict_grna.tsv",sep = "/"), stringsAsFactors = FALSE)
  cbc_gbc_grna <- unique(cbc_gbc_grna)
  data_preprocess <- function(perturb_seq, cbc_gbc) {
    cell_KO_hash = hash()
    for (i in 1:nrow(cbc_gbc)) {
      cbc = cbc_gbc[i, 1]
      gbc = cbc_gbc[i, 2]
      if (has.key(cbc, cell_KO_hash)) {
        cell_KO_hash[cbc] = paste(cell_KO_hash[[cbc]],gbc, sep = ",")
      }
      else {
        cell_KO_hash[cbc] = gbc
      }
    }
    sample_info = c()
    j = 1
    k = 1
    nogbc_col = c()
    for (i in 1:ncol(perturb_seq)) {
      if (!is.null(cell_KO_hash[[colnames(perturb_seq)[i]]])) {
        sample_info[k] = cell_KO_hash[[colnames(perturb_seq)[i]]]
        k = k + 1
      }
      else {
        nogbc_col[j] <- i
        j = j + 1
      }
    }
    perturb_seq <- perturb_seq[, -nogbc_col]
    names(sample_info) = colnames(perturb_seq)
    for (i in 1:length(sample_info)) {
      sample_info_arr <- unlist(strsplit(sample_info[i],","))
      if (length(sample_info_arr) > 1) {
        sample_info_arr <- sort(sample_info_arr)
        begin=1
        if(!(sample_info_arr[1]=="CTRL")){
          sortedMultiKO=sample_info_arr[1]
          begin=2
        }
        else{
          sortedMultiKO=sample_info_arr[2]
          begin=3
        }
        if(length(sample_info_arr)>=3){
          for (j in begin:length(sample_info_arr)) {
            sortedMultiKO = paste(sortedMultiKO, sample_info_arr[j],sep = ",")
          }
        }
        sample_info[i] = sortedMultiKO
      }
    }
    perturb_seq <- perturb_seq[!str_detect(row.names(perturb_seq),"^MRP"), ]
    perturb_seq <- perturb_seq[!str_detect(row.names(perturb_seq),"^RP"), ]
    return(list(perturb_data = perturb_seq, sample_info = sample_info))
  }
  perturb_seq_list <- data_preprocess(perturb_seq, cbc_gbc)
  perturb_seq_grna_list <- data_preprocess(perturb_seq, cbc_gbc_grna)
  if(is.null(sample_info_batch)){
    perturb_list = list("expression_profile" = perturb_seq_list$perturb_data, "sample_info_gene" = perturb_seq_list$sample_info, "sample_info_sgRNA" = perturb_seq_grna_list$sample_info)
  }else{
    perturb_list = list("expression_profile" = perturb_seq_list$perturb_data, "sample_info_gene" = perturb_seq_list$sample_info, "sample_info_sgRNA" = perturb_seq_grna_list$sample_info,"sample_info_batch"=sample_info_batch)
  }
  return(perturb_list)
}
# no matter what the original data format is, users can handle the original data to this format.
Input_preprocess<-function(expression_profile,sample_info_gene,sample_info_sgRNA,sample_info_batch=NULL){
  require(stringr)
  options(warn=1)
  if(ncol(expression_profile)!=length(sample_info_gene)){
    warning("expression_profile and sample_info_gene have different length, please check and try again.")
  }
  if(ncol(expression_profile)!=length(sample_info_sgRNA)){
    warning("expression_profile and sample_info_sgrna have different length, please check and try again.")
  }
  if(length(sample_info_gene)!=length(sample_info_sgRNA)){
    warning("sample_info_gene and sample_info_sgrna have different length, please check and try again.")
  }
  if(!is.null(sample_info_batch)){
    if(ncol(expression_profile)!=length(sample_info_batch)){
      warning("expression_profile and sample_info_batch have different length, please check and try again.")
    }
    if(length(sample_info_batch)!=length(sample_info_gene)){
      warning("sample_info_batch and sample_info_gene have different length, please check and try again.")
    }
    if(length(sample_info_batch)!=length(sample_info_sgRNA)){
      warning("sample_info_batch and sample_info_sgRNA have different length, please check and try again.")
    }
  }
  expression_profile<-expression_profile[(!str_detect(row.names(expression_profile),"^MRP")) & (!str_detect(row.names(expression_profile),"^RP")),]#filter mitochondrial ribosomal protein and ribosomal protein
  if(is.null(sample_info_batch)){
    return(list("expression_profile"=expression_profile,"sample_info_gene"=sample_info_gene,"sample_info_sgRNA"=sample_info_sgRNA))
  }else{
    return(list("expression_profile"=expression_profile,"sample_info_gene"=sample_info_gene,"sample_info_sgRNA"=sample_info_sgRNA,"sample_info_batch"=sample_info_batch))
  }
}
# ********************** quality control of data
Cell_qc<-function(expression_profile,sample_info_gene,species="Hs",gene_low=500,gene_high=10000,mito_high=0.1,umi_low=1500,umi_high=Inf,plot=FALSE,plot_path="~/quality_control.pdf"){
  require(Seurat)
  require(stringr)
  sample_info_gene<-sample_info_gene[names(sample_info_gene) %in% colnames(expression_profile)]
  expression_profile <-expression_profile[,!duplicated(colnames(expression_profile))]
  KO <- new("seurat", raw.data = expression_profile)
  KO <- Setup(KO, min.cells = 0,min.gene=0,do.logNormalize = TRUE, total.expr = 1e4, project = "expression_profile")
  if(species=="Hs"){
    mito.genes <- grep("^MT-", rownames(KO@data), value = TRUE)
  }else if(species=="Mm"){
    mito.genes <- grep("^mt-", rownames(KO@data), value = TRUE)
  }else{
    stop("species should be 'Mm' or 'Hs'")
  }
  percent.mito <- colSums(as.matrix(expm1(KO@data[mito.genes, ])))/colSums(as.matrix(expm1(KO@data)))
  KO <- AddMetaData(KO, percent.mito, "percent.mito")
  KO <- SubsetData(KO, subset.name = "nUMI")
  KO <- SubsetData(KO, subset.name = "nGene")
  KO <- SubsetData(KO, subset.name = "percent.mito")
  if(species=="Hs"){
    KO@data<-KO@data[!str_detect(row.names(KO@data),"^MT-"),]

  }else if(species=="Mm"){
    KO@data<-KO@data[!str_detect(row.names(KO@data),"^mt-"),]
  }
  KO@data.info$QC<-"filter"
  for(i in 1:nrow(KO@data.info)){
    if(KO@data.info$nGene[i]>gene_low & KO@data.info$nGene[i]<gene_high & KO@data.info$percent.mito[i]<mito_high & KO@data.info$nUMI[i]>umi_low & KO@data.info$nUMI[i]<umi_high){
      KO@data.info$QC[i]<-"retain"
    }
  }
  if(plot==TRUE){
    pdf(file=plot_path)
    par(mfrow=c(1,3))
    hist(KO@data.info$nGene,breaks = 20,freq=FALSE,xlab = "Gene numbers",ylab = "Density",main = "Gene numbers distribution")
    lines(density(KO@data.info$nGene),col="red",lwd=1)
    hist(KO@data.info$nUMI,breaks = 20,freq=FALSE,xlab = "UMI numbers",ylab = "Density",main = "UMI numbers distribution")
    lines(density(KO@data.info$nUMI),col="red",lwd=1)
    hist(KO@data.info$percent.mito,breaks = 20,freq=FALSE,xlab = "Percent of mito",ylab = "Density",main = "Percent of mito distribution")
    lines(density(KO@data.info$percent.mito),col="red",lwd=1)
    VlnPlot(KO, c("nGene", "nUMI", "percent.mito"), nCol = 3,size.use=0.001,cols.use=c("palegreen","pink"),group.by = "QC")
    dev.off()
  }
  KO_filter<-as.matrix(KO@data[,which(KO@data.info$QC=="retain")])
  KO_data_qc<-KO_filter#have adopted total_expr=10000 and log
  sample_info_gene<-sample_info_gene[names(sample_info_gene) %in% colnames(KO_data_qc)]
  return(list("expression_profile_qc"=KO_data_qc,"sample_info_gene_qc"=sample_info_gene))
}
# ***********************batch effect if there is batch and you think it may affect the example
Batch_adjust<-function(expression_profile,sample_info_batch,sample_info_gene,plot=TRUE,plot_path="~/batch_check.pdf"){
  require(BatchQC)
  require(sva)
  gene_counts_sum<-apply(expression_profile,1,sum)
  expression_profile<-expression_profile[which(gene_counts_sum>0),]
  sample_info_batch<-sample_info_batch[names(sample_info_gene)]
  pdata <- data.frame(sample_info_batch, sample_info_gene)
  mod = model.matrix(~as.factor(sample_info_gene), data = pdata)
  batchTest<-combatPlot(expression_profile,sample_info_batch,mod)
  pValue<-batchTest$p.value
  batch_explain<-batchqc_explained_variation(expression_profile, sample_info_gene, sample_info_batch)
  pdf(file=plot_path)
  print(pValue)
  if(pValue<=0.05){
    if(plot==TRUE){
      par(mfrow=c(1,2),cex.axis=0.5)
      p1=boxplot(batch_explain$explained_variation)
      print(p1)
      # pca figure may take long time
      noPrint<-batchqc_pca(expression_profile,sample_info_batch)
    }
  }
  else{
    if(plot==TRUE){
      par(mfrow=c(2,2))
      p2=boxplot(batch_explain$explained_variation,cex.axis=0.5)
      print(p2)
      noPrint<-batchqc_pca(expression_profile,sample_info_batch,cex.axis=0.5)
    }
    expression_profile<-ComBat(expression_profile=expression_profile, sample_info_batch=sample_info_batch,mod=mod)
    batch_explain<-batchqc_explained_variation(expression_profile, sample_info_gene, sample_info_batch)
    if(plot==TRUE){
      p3=boxplot(batch_explain$explained_variation,cex.axis=0.5)
      print(p3)
      noPrint<-batchqc_pca(expression_profile,sample_info_batch,cex.axis=0.5)
    }
  }
  dev.off()
  return(expression_profile)
}
# ********************   sgrna efficiency assessment and filtering
Cell_filtering<-function(expression_profile,sample_info_gene,sample_info_sgRNA,nonzero=0.01,grna_cell_num=10,fold_change=0.5,cellNum=30,plot=FALSE,plot_path_zeroRatio="~/nonzero_ratio.pdf",plot_path_sgRNAefficiency="~/sgRNA_efficiency.pdf",plot_path_capturePhenotype="~/phenotype_capture.pdf",plot_path_KO_efficiency="~/KO_efficiency.pdf"){
  require(reshape2)
  require(ggplot2)
  require(stringr)
  require(hash)
  #building a hash from grna name to gene name
  grna_gene_hash<-hash()
  sample_info_sgRNA<-sample_info_sgRNA[intersect(names(sample_info_sgRNA),names(sample_info_gene))]
  sample_info_gene<-sample_info_gene[intersect(names(sample_info_sgRNA),names(sample_info_gene))]
  if(is.null(sample_info_sgRNA)){
    stop("sample_info_sgRNA is null, please check it.")
  }
  for(i in 1:length(sample_info_sgRNA)){
    grna_gene_hash[sample_info_sgRNA[i]]<-sample_info_gene[i]
  }
  sample_info_grna_KO<-sample_info_sgRNA[!str_detect(sample_info_sgRNA,"CTRL")]
  #turning combined sample label to single sample label
  sample_info_grna_KO_split<-c()
  for(i in 1:length(sample_info_grna_KO)){
    ko_split<-unlist(str_split(sample_info_grna_KO[i],","))
    names(ko_split)<-rep(names(sample_info_grna_KO)[i],length(ko_split))
    sample_info_grna_KO_split<-c(sample_info_grna_KO_split,ko_split)
  }
  sample_info_grna_CTRL<-sample_info_sgRNA[str_detect(sample_info_sgRNA,"CTRL")]
  sample_info_gene_KO<-sample_info_gene[sample_info_gene!="CTRL"]
  sample_info_gene_KO_split<-c()
  for(i in 1:length(sample_info_gene_KO)){
    ko_split<-unlist(str_split(sample_info_gene_KO[i],","))
    names(ko_split)<-rep(names(sample_info_gene_KO)[i],length(ko_split))
    sample_info_gene_KO_split<-c(sample_info_gene_KO_split,ko_split)
  }
  sample_info_gene_CTRL<-sample_info_gene[sample_info_gene=="CTRL"]
  gene_express_KO<-expression_profile[(row.names(expression_profile) %in% sample_info_gene_KO_split),(colnames(expression_profile) %in% names(sample_info_gene_KO_split))]
  gene_express_CTRL<-expression_profile[(row.names(expression_profile) %in% sample_info_gene_KO_split),(colnames(expression_profile) %in% names(sample_info_gene_CTRL))]
  ######step 1: Filtering out cells with an extremely large proportions of zero knockout expression in the control cells.
  print("Filtering cells by considering ratio of nonzero knockout expression value in control cells.")
  zeroRatio<-matrix(rep(NA,3*nrow(gene_express_CTRL)),nrow(gene_express_CTRL),3)
  colnames(zeroRatio)<-c("zero","nonzero","gene_order")
  row.names(zeroRatio)<-row.names(gene_express_CTRL)
  for(i in 1:nrow(gene_express_CTRL)){
    zero_ratio<-length(gene_express_CTRL[i,][which(gene_express_CTRL[i,]==0)])/length(gene_express_CTRL[i,])
    zeroRatio[i,1]<-zero_ratio
    zeroRatio[i,2]<-1-zero_ratio
    if(1-zero_ratio<=nonzero){
      gene_express_CTRL[i,]=NA
    }
  }
  zeroRatio<-zeroRatio[order(zeroRatio[,2]),]
  zeroRatio[,3]<-1:nrow(zeroRatio)
  zeroRatio<-melt(zeroRatio)
  zeroRatio<-zeroRatio[zeroRatio[,2]!="gene_order",]
  colnames(zeroRatio)<-c("gene","label","value")
  nonzeroRatio<-zeroRatio[which(zeroRatio[,2]=="nonzero"),]
  nonzeroRatio<-as.character(nonzeroRatio[which(nonzeroRatio[,"value"]<nonzero),"gene"])
  ######step 2: Reducing false positive rate of gene knockout. performing the filterings to make sure log2_fold_change is less than log2(fold_change)
  print("Filtering cells by considering sgRNA knockout efficiency.")
  gene_express_CTRL<-na.omit(gene_express_CTRL)
  gene_express_KO<-gene_express_KO[row.names(gene_express_CTRL),]
  gene_express_KO<-na.omit(gene_express_KO)
  gene_express_CTRL<-na.omit(gene_express_CTRL)
  gene_express_CTRL<-gene_express_CTRL[sort(row.names(gene_express_CTRL)),]
  grna_express_KO<-gene_express_KO
  grna<-sort(unique(sample_info_grna_KO_split))
  gene<-sort(row.names(gene_express_KO))
  #calculating log2_fold_change between case and control
  grna_gene_express<-matrix(rep(NA,length(grna)*4),length(grna),4)
  row.names(grna_gene_express)<-grna
  grna_gene_express_filtered<-matrix(rep(NA,length(grna)*4),length(grna),4)
  row.names(grna_gene_express_filtered)<-grna
  colnames(grna_gene_express)<-c("sample_num","ko_mean_express","ctrl_mean_express","log2_fold_change")
  colnames(grna_gene_express_filtered)<-c("sample_num_filtered","ko_mean_express_filtered","ctrl_mean_express_filtered","log2_fold_change")
  grna_list=list()
  grna_ctrl_list=list()
  k=1
  l=1
  m=1
  KO_gene_label=c()
  CTRL_gene_label=c()
  for(i in gene){
    #the name of sgRNA must like this "knockoutName_xxx"
    i2<-paste(i,"_",sep="")
    for(j in grna){
      if(str_detect(j,i2)){
        grna_express=grna_express_KO[i,names(sample_info_grna_KO_split[sample_info_grna_KO_split==j])]
        grna_list[[k]]=grna_express
        grna_ctrl_list[[k]]=as.numeric(gene_express_CTRL[i,])
        names(grna_list)[[k]]=j
        names(grna_ctrl_list)[[k]]=j
        k=k+1
        grna_gene_express[j,"sample_num"]=length(grna_express)
        KO_gene_label[l:(l+length(grna_express)-1)]=rep(i,each=length(grna_express))
        CTRL_gene_label[m:(m+length(gene_express_CTRL[i,])-1)]=rep(i,each=length(gene_express_CTRL[i,]))
        l=l+length(grna_express)
        m=m+length(gene_express_CTRL[i,])
        grna_gene_express[j,"ko_mean_express"]=mean(grna_express)
        grna_gene_express[j,"ctrl_mean_express"]=mean(gene_express_CTRL[i,])
        grna_gene_express[j,"log2_fold_change"]=log2(grna_gene_express[j,"ko_mean_express"]/grna_gene_express[j,"ctrl_mean_express"])
      }
    }
  }
  grna_gene_express=na.omit(grna_gene_express)
  grna_gene_express=grna_gene_express[grna_gene_express[,1]>=grna_cell_num,]
  grna_gene_fold<-as.matrix(grna_gene_express[,4])
  grna_gene_dataframe<-data.frame(row.names(grna_gene_fold),grna_gene_fold[,1])
  grna_gene_dataframe$label<-"less"
  grna_gene_dataframe[grna_gene_dataframe[,2]>0,3]<-"greater"
  colnames(grna_gene_dataframe)=c("grna_name","log2_fold_change","label")
  grna_gene_dataframe2<-grna_gene_dataframe[grna_gene_dataframe$log2_fold_change<1,]
  grna_gene_express_filtered<-grna_gene_express_filtered[row.names(grna_gene_dataframe2),]
  for(i in 1:nrow(grna_gene_dataframe2)){
    if(grna_gene_dataframe2$log2_fold_change[i]>log2(fold_change)){
      for(j in seq(1,0,by=-0.01)){
        grna_name<-as.character(grna_gene_dataframe2$grna_name[i])
        medi=quantile(grna_list[[grna_name]],j)
        grna_list[[grna_name]]=grna_list[[grna_name]][grna_list[[grna_name]]<=medi]
        fold_change_new<-mean(grna_list[[grna_name]])/mean(grna_ctrl_list[[grna_name]])
        log_fold_change_new=log2(fold_change_new)
        grna_gene_dataframe2$log2_fold_change[i]<-log_fold_change_new
        if(grna_gene_dataframe2$log2_fold_change[i]<=-1){
          grna_gene_express_filtered[i,"sample_num_filtered"]<-length(grna_list[[grna_name]])
          grna_gene_express_filtered[i,"ko_mean_express_filtered"]<-mean(grna_list[[grna_name]])
          grna_gene_express_filtered[i,"ctrl_mean_express_filtered"]<-mean(grna_ctrl_list[[grna_name]])
          grna_gene_express_filtered[i,"log2_fold_change"]<-log_fold_change_new
          break
        }
      }
    }
    else{
      grna_gene_express_filtered[i,"sample_num_filtered"]<-grna_gene_express[grna_gene_dataframe2$grna_name[i],"sample_num"]
      grna_gene_express_filtered[i,"ko_mean_express_filtered"]<-grna_gene_express[grna_gene_dataframe2$grna_name[i],"ko_mean_express"]
      grna_gene_express_filtered[i,"ctrl_mean_express_filtered"]<-grna_gene_express[grna_gene_dataframe2$grna_name[i],"ctrl_mean_express"]
      grna_gene_express_filtered[i,"log2_fold_change"]<-grna_gene_express[grna_gene_dataframe2$grna_name[i],"log2_fold_change"]
    }
  }
  grna_gene_express_filtered[is.infinite(grna_gene_express_filtered[,"log2_fold_change"]),"log2_fold_change"]<-mean(grna_gene_express_filtered[is.finite(grna_gene_express_filtered[,"log2_fold_change"]),"log2_fold_change"])
  #obtainning grna efficiency
  grna_efficiency<-grna_gene_express_filtered[,"log2_fold_change"]
  #calculating KO efficiency with grna efficiency
  print("Obtaining the overall knockout efficiency.")
  sample_ctrl<-colnames(gene_express_CTRL)
  sample_KO_filtered<-c()
  sample_info_gene_KO_filtered<-c()
  sample_info_grna_KO_filtered<-c()
  for(i in row.names(grna_gene_express_filtered)){
    sample_KO_filtered<-c(sample_KO_filtered,as.character(names(grna_list[[i]])))
  }
  sample_KO_filtered<-unique(sample_KO_filtered)
  sample_info_gene_KO_filtered<-sample_info_gene[sample_KO_filtered]
  sample_info_grna_KO_filtered<-sample_info_sgRNA[sample_KO_filtered]
  sample_gene_KO_table<-table(sample_info_gene_KO_filtered)
  sample_grna_KO_table<-table(sample_info_grna_KO_filtered)
  KO_efficiency<-rep(0,length(sample_gene_KO_table))
  names(KO_efficiency)<-names(sample_gene_KO_table)
  for(i in 1:length(sample_gene_KO_table)){
    ko_split<-unlist(str_split(names(sample_gene_KO_table)[i],","))
    ko_split<-sort(ko_split)
    for(j in 1:length(sample_grna_KO_table)){
      ko_grna_split<-unlist(str_split(names(sample_grna_KO_table)[j],","))
      ko_grna_split_gene<-c()
      for(k in ko_grna_split){
        ko_grna_split_gene<-c(ko_grna_split_gene,as.character(grna_gene_hash[[k]]))
      }
      ko_grna_split_gene<-sort(ko_grna_split_gene)
      if(identical(ko_split,ko_grna_split_gene)){
        KO_efficiency[i]<-KO_efficiency[i]+(sample_grna_KO_table[j]/sample_gene_KO_table[i])*abs(mean(grna_efficiency[ko_grna_split]))
      }
    }
  }
  KO_efficiency<-na.omit(KO_efficiency)
  sample_info_gene_filtered_qc_zr_se<-sample_info_gene[c(sample_KO_filtered,sample_ctrl)]
  expression_profile_filtered_qc_zr_se<-expression_profile[,names(sample_info_gene_filtered_qc_zr_se)]
  ######step 3:Filtering out knockout cells without sufficient cell number to capture the corresponding perturbation phenotype
  print("Filtering cells by considering phenotype capture.")
  sample_info_gene_filtered_qc_zr_se<-sample_info_gene_filtered_qc_zr_se[names(sample_info_gene_filtered_qc_zr_se) %in% colnames(expression_profile_filtered_qc_zr_se)]
  cellNum_eachKO<-table(sample_info_gene_filtered_qc_zr_se)
  cellNum_eachKO_filtered<-cellNum_eachKO[cellNum_eachKO>=cellNum]
  sample_info_gene_filtered_qc_zr_se_pc<-sample_info_gene_filtered_qc_zr_se[sample_info_gene_filtered_qc_zr_se %in% names(cellNum_eachKO_filtered)]
  expression_profile_filtered_qc_zr_se_pc<-expression_profile_filtered_qc_zr_se[,intersect(colnames(expression_profile_filtered_qc_zr_se),names(sample_info_gene_filtered_qc_zr_se_pc))]
  gene_counts_sum<-apply(expression_profile_filtered_qc_zr_se_pc,1,sum)
  expression_profile_filtered_qc_zr_se_pc<-expression_profile_filtered_qc_zr_se_pc[which(gene_counts_sum>0),]
  #
  if(plot==TRUE){
    #plot nonzero ratio in control cells
    zeroRatio_new<-zeroRatio[zeroRatio$label=="nonzero",-2]
    zeroRatio_new2<-zeroRatio_new[order(-zeroRatio_new[,2]),]
    colnames(zeroRatio_new2)<-c("knock","ratio")
    pdf(file=plot_path_zeroRatio)
    p=ggplot(zeroRatio_new2,aes(x=knock,y=ratio))+geom_bar(stat = "identity", position="stack",width = 0.5)+theme(axis.text.y=element_text(size=5.5,face="bold"))+theme_bw()+theme(axis.text.x=element_text(angle=90,size=8))+geom_hline(aes(yintercept=nonzero),col="red",linetype="dashed")+scale_y_continuous(breaks=c(0,nonzero,0.25,0.5,0.75,1))+xlab("Knockout")+ylab("Ratio of nonzero value")
    print(p)
    dev.off()
    #plot grna efficiency
    pdf(file=plot_path_sgRNAefficiency)
    p=ggplot(grna_gene_dataframe,aes(grna_name,log2_fold_change,fill=label))+geom_bar(stat="identity",position = "identity")+theme(axis.text.x=element_text(angle=90,size=6))+geom_hline(aes(yintercept=1),col="red",linetype="dashed")+geom_hline(aes(yintercept=-1),col="red",linetype="dashed")+xlab("sgRNA")+ylab("log2 fold change")
    print(p)
    dev.off()
    #plot cells for phenotype capture
    phenotype<-table(sample_info_gene_filtered_qc_zr_se)
    phenotype_plot<-as.matrix(phenotype)
    phenotype_plot<-as.data.frame(phenotype_plot)
    phenotype_plot$knockout<-row.names(phenotype_plot)
    phenotype_plot2<-phenotype_plot[which(phenotype_plot[,1]>0),]
    colnames(phenotype_plot2)<-c("cellNum","knockout")
    pdf(plot_path_capturePhenotype)
    p=ggplot(phenotype_plot2,aes(x=knockout,y=cellNum,group=1))+geom_line()+geom_point(size=3)+xlab("Knockout")+ylab("Cell number")+theme_bw()+theme(axis.text.x=element_text(angle=90,size=8))+geom_hline(aes(yintercept=30),col="red",linetype="dashed")+scale_y_continuous(breaks=c(0,cellNum,100,200,300,400,500))+
      geom_text(aes(label = cellNum,vjust = -0.8, hjust = 0.5),size=3,show.legend=FALSE)
    print(p)
    dev.off()
    #plot knockout efficiency
    correction_KO<-round(KO_efficiency,2)
    correction_KO_matrix<-as.matrix(correction_KO)
    correction_KO_frame<-as.data.frame(correction_KO_matrix)
    correction_KO_frame$knockout<-row.names(correction_KO_frame)
    colnames(correction_KO_frame)<-c("KO_efficiency","knockout")
    pdf(file=plot_path_KO_efficiency)
    p=ggplot(correction_KO_frame,aes(x=knockout,y=KO_efficiency))+geom_bar(stat = "identity", position="stack",width = 0.5)+theme(axis.text.y=element_text(size=5.5,face="bold"))+theme_bw()+theme(axis.text.x=element_text(angle=90,size=8))+xlab("Knockout")+ylab("Overall knockout efficiency")+
      geom_text(aes(label = KO_efficiency,vjust = -0.8, hjust = 0.5),show.legend=FALSE)
    print(p)
    dev.off()
  }
  Filtering_And_KOefficiency<-list("expression_profile_qc_zr_se"=expression_profile_filtered_qc_zr_se,"sample_info_gene_qc_zr_se"=sample_info_gene_filtered_qc_zr_se,"expression_profile_qc_zr_se_pc"=expression_profile_filtered_qc_zr_se_pc,"sample_info_gene_qc_zr_se_pc"=sample_info_gene_filtered_qc_zr_se_pc,"KO_efficiency"=KO_efficiency,"sgRNA_efficiency"=grna_efficiency,"nonzeroRatio"=nonzeroRatio)
  return(Filtering_And_KOefficiency)
}
# ********************   plotting all information component off cells
Plot_filtering_overview<-function(sample_info_gene_original=c(),sample_info_gene_qc=c(),nonzeroRatio=c(),sample_info_gene_qc_zr_se=c(),sample_info_gene_qc_zr_se_pc=c(),plot_path="~/overview_of_cell_filterings.pdf"){
  require(ggplot2)
  original<-table(sample_info_gene_original)
  filter1_qc<-table(sample_info_gene_qc)
  filter2_nonzeroRatio<-nonzeroRatio
  filter3_sgrna<-table(sample_info_gene_qc_zr_se)
  filter4_phenotype<-table(sample_info_gene_qc_zr_se_pc)
  original<-sort(original)
  filter1<-original
  for(n in names(original)){
    if(n %in% names(filter1_qc)){
      filter1[n]=filter1_qc[n]
    }
    else{
      filter1[n]=0
    }
  }
  filter2<-filter1
  filter2[filter2_nonzeroRatio]=0
  filter3=filter2
  for(n in names(filter2)){
    if(n %in% names(filter3_sgrna)){
      filter3[n]=filter3_sgrna[n]
    }
    else{
      filter3[n]=0
    }
  }
  filter4=filter3
  for(n in names(filter3)){
    if(n %in% names(filter4_phenotype)){
      filter4[n]=filter4_phenotype[n]
    }
    else{
      filter4[n]=0
    }
  }
  filter_m<-matrix(c(original,filter1,filter2,filter3,filter4),5,byrow = TRUE)
  colnames(filter_m)<-names(original)
  row.names(filter_m)<-c("original_number","filter1_qualityControl","filter2_zeroRatio","filter3_sgRNA","filter4_phenotypeCapture")
  filter_m_median<-matrix(rep(0,5*ncol(filter_m)),5)
  filter_m_median[1,]<-filter_m[1,]-filter_m[2,]
  filter_m_median[2,]<-filter_m[2,]-filter_m[3,]
  filter_m_median[3,]<-filter_m[3,]-filter_m[4,]
  filter_m_median[4,]<-filter_m[4,]-filter_m[5,]
  filter_m_median[5,]<-filter_m[5,]
  colnames(filter_m_median)<-colnames(filter_m)
  row.names(filter_m_median)<-c("Cells filtered at Step 1","Cells filtered at Step 2","Cells filtered at Step 3","Cells filtered at Step 4","Cells retained for subsequent analysis")
  filter_median_melt<-melt(filter_m_median)
  colnames(filter_median_melt)<-c("Cell_components","Knockout","Cell_number")
  pdf(plot_path)
  p<-ggplot(filter_median_melt,aes(x=Knockout,weight=Cell_number,fill=Cell_components))+geom_bar(position="stack")+xlab("Knockout")+ylab("Cell number")+theme(axis.text.x=element_text(angle=90,size=8))
  plot(p)
  dev.off()
  return(filter_median_melt)
}
################################     Model building          ##############################
# ************************   obtaining high dispersion different genes
Get_high_var_genes<-function(expression_profile,sample_info_gene,x.low.cutoff=0.0125,x.high.cutoff=5,y.cutoff=1,do.spike=FALSE,num.bin=30,plot=FALSE,plot_path="~/get_high_var_genes.pdf"){
  logVarDivMean=function(x) return(log(var(exp(x)-1)/mean(exp(x)-1)))
  expMean=function(x) return(log(mean(exp(x)-1)+1))
  data_norm<-function(xy,num.bin){
    seperate<-seq(0,max(xy[,1]),length.out=num.bin+1)
    for(i in 2:length(seperate)){
      xy[which(xy[,1]>seperate[i-1] & xy[,1]<=seperate[i]),3]<-(xy[which(xy[,1]>seperate[i-1] & xy[,1]<=seperate[i]),2]-mean(xy[which(xy[,1]>seperate[i-1] & xy[,1]<=seperate[i]),2]))/sd(xy[which(xy[,1]>seperate[i-1] & xy[,1]<=seperate[i]),2])
    }
    return(xy)
  }
  label=c()
  j=1
  for(i in 1:length(sample_info_gene)){
    if(sample_info_gene[i]=="CTRL"){
      label[j]="ctrl"
      j=j+1
    }
    else{
      label[j]="ko"
      j=j+1
    }
  }
  if(do.spike){
    require(DESeq2)
    geneTypes <- factor( c( endoGene="NNNN", ERCC="ERCC" )[substr( rownames(expression_profile), 1, 4 ) ] )
    countsERCC <- expression_profile[ which( geneTypes=="ERCC" ), ]
    countsEndo<-expression_profile[-which(geneTypes=="ERCC"),]
    lengthsEndo <- expression_profile[ -which( geneTypes=="ERCC" ), 1 ]
    lengthsERCC <- expression_profile[ which( geneTypes=="ERCC" ), 1 ]
    sfERCC <- estimateSizeFactorsForMatrix(countsERCC)
    sfEndo<-sfERCC
    nCountsEndo<-t(t(countsEndo)/sfEndo)
    expression_profile=nCountsEndo
  }
  nCountsEndo_ko<-expression_profile[,which(label=="ko")]
  nCountsEndo_ctrl<-expression_profile[,which(label=="ctrl")]
  data.x_ko<-apply(nCountsEndo_ko,1,expMean)
  data.y_ko<-apply(nCountsEndo_ko,1,logVarDivMean)
  data.y_ko[is.nan(data.y_ko)]<-0
  data.y_ctrl<-apply(nCountsEndo_ctrl,1,logVarDivMean)
  data.y_ctrl[is.nan(data.y_ctrl)]<-0
  data.y_diff<-abs(data.y_ko-data.y_ctrl)

  diff_xy<-cbind(data.x_ko,data.y_diff,0,1)
  row.names(diff_xy)<-row.names(nCountsEndo_ko)
  colnames(diff_xy)<-c("data.x","data.y","data.norm.y","vargene")
  diff_xy<-data_norm(diff_xy,num.bin)
  diff_xy[is.na(diff_xy[,3]),3]<-0
  diff_xy[which(diff_xy[,1]>x.low.cutoff & diff_xy[,1]<x.high.cutoff & diff_xy[,3]>y.cutoff),4]<-2
  if(plot==TRUE){
    pdf(file=plot_path)
    plot(diff_xy[,1],diff_xy[,3],type="p",xlab="Average expression",pch=16,ylab="Dispersion difference",col=diff_xy[,4])
    dev.off()
  }
  singleCellCRISPRscreen_vargene<-expression_profile[which(diff_xy[,4]==2),]
  return(singleCellCRISPRscreen_vargene)
}
# ************************   selecting the optimal topic number automatically
Get_topics<-function(expression_profile_var_gene,sample_info_gene,topic_number_min=3,topic_number_max=8,alpha=0.5,seed_num=2017,burnin=1000,thin=100,iter=1000,plot=FALSE,plot_path="~/select_topic_number.pdf"){
  require(slam)
  require(topicmodels)
  require(ggplot2)
  #adjust data for topic model
  label=c()
  j=1
  for(i in 1:length(sample_info_gene)){
    if(sample_info_gene[i]=="CTRL"){
      label[j]="ctrl"
      j=j+1
    }
    else{
      label[j]="ko"
      j=j+1
    }
  }
  expression_profile_var_gene<-expression_profile_var_gene+1
  nCountsEndo_ctrl<-expression_profile_var_gene[,which(label=="ctrl")]
  for(i in 1:nrow(expression_profile_var_gene)){
    expression_profile_var_gene[i,]=(expression_profile_var_gene[i,]-mean(nCountsEndo_ctrl[i,]))/mean(nCountsEndo_ctrl[i,])
  }
  expression_profile_var_gene<-round((expression_profile_var_gene+abs(min(expression_profile_var_gene)))*10)
  #
  gene_counts_sum<-apply(expression_profile_var_gene,1,sum)
  expression_profile_var_gene<-expression_profile_var_gene[which(gene_counts_sum>0),]
  control=list(seed=seed_num,burnin=burnin,thin=thin,iter=iter)
  dtm<-as.simple_triplet_matrix(t(expression_profile_var_gene))
  topic_model_list=list()
  i=1
  specificity_score<-c()
  purity_score<-c()
  combination_score<-c()
  topic_name=c()
  for(k in topic_number_min:topic_number_max){
    print(paste("The scope of topic number you choose is between",topic_number_min,"and",topic_number_max,",now the evaluating topic number is",k,sep=" "))
    topic_model=LDA(dtm,k=k,method="Gibbs",control=control)
    topic_model_list[[i]]=topic_model
    topic_model<-topic_model_list[[i]]@gamma
    topic_num<-ncol(topic_model)
    col_varDivSq<-apply(topic_model,2,var)/apply(topic_model,2,mean)^2
    row_var<-apply(topic_model,1,var)
    specificity_score[i]=log(mean(col_varDivSq))
    purity_score[i]=log(mean(row_var))
    combination_score[i]=alpha*specificity_score[i]+(1-alpha)*purity_score[i]
    topic_name[i]=topic_num
    i<-i+1
  }
  if(plot==TRUE){
    pdf(file=plot_path)
    Topic_number<-factor(topic_name,levels = topic_name)
    selectTopic_dataFrame<-data.frame(Topic_number,Score=combination_score)
    p=ggplot(selectTopic_dataFrame,aes(x=Topic_number,y=Score,group=1))+theme(axis.text.x=element_text(size=10))+geom_line()+geom_point(size=6)+xlab("Topic number")+ylab("Score")
    print(p)
    dev.off()
  }
  m=order(combination_score,decreasing = TRUE)[1]
  return(topic_model_list[[m]])
}
# **********************   plotting heatmap between cells and topics
Plot_cell_topic<-function(model,plot_path="~/distribution_cell_in_topics.pdf"){
  require(gplots)
  pdf(file=plot_path)
  topic_cell<-model@gamma
  quan_up<-quantile(topic_cell,0.98)
  quan_down<-quantile(topic_cell,0.2)
  for(i in 1:nrow(topic_cell)){
    for(j in 1:ncol(topic_cell)){
      if(topic_cell[i,j]>=quan_up){
        topic_cell[i,j]<-quan_up
      }
      if(topic_cell[i,j]<quan_down){
        topic_cell[i,j]<-quan_down
      }
    }
  }
  colnames(topic_cell)<-paste('Topic',1:ncol(topic_cell),sep=" ")
  heatmap.2(topic_cell,col=bluered,Colv=FALSE,dendrogram="row",labRow="",labCol=colnames(topic_cell),cexCol = 1,adjCol=c(NA,1),xlab="Topic",ylab="Cell",srtCol = 45, key=TRUE,trace="none", breaks=seq.int(from = min(topic_cell), to = max(topic_cell), length.out = 50), main = 'Distribution of cells for each topic',hclustfun=hclust)
  dev.off()
}
# *********************   annotating each topic's functions for Hs(homo sapiens) or Mm(mus musculus)
Topic_func_anno <-function(model,species="Hs",topGene_percent=0.2,topNum=5,FDR=0.1,plot=TRUE,plot_path="~/topic_annotation.pdf"){
  require(clusterProfiler)
  require(reshape2)
  require(Biostrings)
  require(dplyr)
  if(species=="Hs"){
    require(org.Hs.eg.db)
    organism="hsa"
    Alia_entrez<-org.Hs.egALIAS2EG
    OrgDb<-"org.Hs.eg.db"
  }else if(species=="Mm"){
    require(org.Mm.eg.db)
    organism="mmu"
    Alia_entrez<-org.Mm.egALIAS2EG
    OrgDb<-"org.Mm.eg.db"
  }else{
    stop("species should be 'Hs' or 'Mm'")
  }
  my_beta <- model@beta
  my_geneName<-model@terms
  colnames(my_beta) <- my_geneName
  rownames(my_beta) <- paste('Topic',1:nrow(my_beta),sep=" ")
  x<-Alia_entrez
  my_geneName<-my_geneName[my_geneName %in% mappedkeys(x)]#filter some unrecognized alia name
  my_beta<-my_beta[,my_geneName]
  my_geneID<-as.list(x[my_geneName])
  col_index=c()
  geneID=c()
  i=1
  for(k in 1:length(my_geneID)){
    if(length(my_geneID[[k]])==1){
      col_index[i]=k
      geneID[i]=my_geneID[[k]]
      i=i+1
    }
  }
  my_beta<-my_beta[,col_index]
  colnames(my_beta)<-geneID
  my_beta <- exp(my_beta)
  data_melt <- melt(my_beta)
  colnames(data_melt)<-c("topics","gene","value")
  data_list<-split(data_melt,data_melt$topics)
  topic_id=list()
  for(i in 1:length(data_list)){
    cutoff<-quantile(data_list[[i]]$value,1-topGene_percent)
    entr_character<-as.character(data_list[[i]][data_list[[i]]$value>=cutoff,"gene"])
    topic_id[i]=list(entr_character)
  }
  topic_id_name=c()
  for(i in 1:length(topic_id)){
    topic_id_name[i]=paste("Topic",i,sep=" ")
  }
  names(topic_id)<-topic_id_name
  Compare_go <- compareCluster(geneCluster=topic_id, fun="enrichGO",OrgDb=OrgDb,ont="BP",minGSSize=1,qvalueCutoff=1,pvalueCutoff=1)
  #choose top-ranked go terms
  enrich_result<-Compare_go@compareClusterResult
  enrich_result$GeneModelCount<-rep(0,nrow(enrich_result))
  my_beta_toOne<-my_beta
  for(i in 1:nrow(my_beta_toOne)){
    my_beta_toOne[i,]=my_beta_toOne[i,]/max(my_beta_toOne[i,])
  }
  for(i in 1:nrow(enrich_result)){
    geneIDArr<-unlist(strsplit(enrich_result$geneID[i],"/"))
    geneModelRatioCount=0
    for(j in 1:length(geneIDArr)){
      geneModelRatioCount=geneModelRatioCount+my_beta_toOne[enrich_result$Cluster[i],geneIDArr[j]]
    }
    enrich_result$GeneModelCount[i]=geneModelRatioCount
    enrich_result$GoGeneAll[i]<-as.numeric(unlist(strsplit(enrich_result$BgRatio[i],"/"))[1])
  }
  enrich_result$GeneModelCount=enrich_result$GeneModelCount*(mean(enrich_result$Count)/mean(enrich_result$GeneModelCount))
  enrich_result$GeneChooseAll<-as.numeric(unlist(strsplit(enrich_result$GeneRatio[1],"/"))[2])
  enrich_result$BgGeneAll<-as.numeric(unlist(strsplit(enrich_result$BgRatio[1],"/"))[2])
  enrich_result$ratio<-enrich_result$GeneModelCount/enrich_result$GeneChooseAll
  enrich_result$p.value_new<-phyper(enrich_result$GeneModelCount-1,enrich_result$GoGeneAll,enrich_result$BgGeneAll-enrich_result$GoGeneAll,enrich_result$GeneChooseAll,lower.tail = FALSE)
  enrich_result$FDR<-p.adjust(enrich_result$p.value_new,method = "fdr")
  enrich_result<-enrich_result[order(enrich_result$Cluster,enrich_result$FDR),]
  ex<-enrich_result %>% group_by(Cluster) %>% summarise(enrich_count=length(Cluster))
  ex<-as.data.frame(ex)
  ex$enrich_count<-cumsum(ex$enrich_count)
  enrich_number=topNum
  enrich_number2=topNum
  enrich_choose_index=c()
  enrich_choose_index=c(1:enrich_number2)
  for(i in 1:(nrow(ex)-1)){
    enrich_choose_index[(enrich_number2+1):(enrich_number2+enrich_number)]=c((ex$enrich_count[i]+1):(ex$enrich_count[i]+enrich_number))
    enrich_number2=enrich_number2+enrich_number
  }
  enrich_result<-enrich_result[enrich_choose_index,]
  topic_annotation_result<-enrich_result[,c("Cluster","Description","ratio","FDR")]
  topic_annotation_result<-topic_annotation_result[topic_annotation_result$FDR<=FDR,]
  #plot the top-ranked go terms of each topic
  if(plot==TRUE){
    require(ggplot2)
    topic_annotation_result<-na.omit(topic_annotation_result)
    topic_annotation_result$Description<-factor(as.character(topic_annotation_result$Description),levels=rev(as.character(topic_annotation_result$Description)))
    pdf(file=plot_path,width=14,height=10)
    p=ggplot(topic_annotation_result,aes(Cluster,Description)) +
      geom_point(aes(colour=FDR),size=6) +
      scale_color_gradient(low = "red", high = "blue")+
      theme_bw() +
      theme(axis.text.x=element_text(angle=45,size=14))+
      theme(axis.text.y=element_text(size=14))+
      theme(axis.title =element_text(size = 0))
    print(p)
    dev.off()
  }
  return(topic_annotation_result)
}
# ********************   evaluating off target effect
Get_offtarget<-function(offTarget_results,expression_profile,sample_info_gene,sample_info_sgRNA){
  offTargetGene<-offTarget_results$offtarget[,c("name","gene")]
  require(hash)
  grna_gene_hash<-hash()
  sample_info_sgRNA<-sample_info_sgRNA[names(sample_info_gene)]
  for(i in 1:length(sample_info_sgRNA)){
    grna_gene_hash[sample_info_sgRNA[i]]<-sample_info_gene[i]
  }
  offTargetGene$target<-NA
  for(i in 1:nrow(offTargetGene)){
    if(has.key(as.character(offTargetGene$name[i]),grna_gene_hash)){
      offTargetGene$target[i]<-grna_gene_hash[[as.character(offTargetGene$name[i])]]
    }
  }
  offTargetGene<-na.omit(offTargetGene)
  offTargetGene<-offTargetGene[offTargetGene$gene!="",]
  offTargetGene<-unique(offTargetGene)
  offTargetGene<-offTargetGene[offTargetGene$target!=offTargetGene$gene,]
  offTargetGene$off<-NA
  control_data<-expression_profile[,names(sample_info_gene[which(sample_info_gene=="CTRL")])]
  for(i in 1:nrow(offTargetGene)){
    offTarget_data<-expression_profile[,names(sample_info_sgRNA[which(sample_info_sgRNA==as.character(offTargetGene$name[i]))])]
    if(as.character(offTargetGene[i,"gene"]) %in% row.names(offTarget_data)){
      offgene<-offTarget_data[as.character(offTargetGene[i,"gene"]),]
    }else{
      offgene<-NA
    }
    if(as.character(offTargetGene[i,"target"]) %in% row.names(offTarget_data)){
      ongene<-offTarget_data[as.character(offTargetGene[i,"target"]),]
    }else{
      ongene<-NA
    }
    if(!is.na(offgene) && !is.na(ongene)){
      r_off<-cor(offgene,ongene)
      offgene_ctrl<-control_data[offTargetGene[i,"gene"],]
      ongene_ctrl<-control_data[offTargetGene[i,"target"],]
      r_ctrl_off<-cor(offgene_ctrl,ongene_ctrl)
      if(!is.na(r_off) && !is.na(r_ctrl_off) && (r_off-r_ctrl_off)/r_ctrl_off>0.05){
        offTargetGene$off[i]="yes"
      }
    }else{
      offTargetGene$off[i]=NA
    }
  }
  offTargetGene<-na.omit(offTargetGene)
  offTargetGene_hash<-hash()
  offTargetGene<-offTargetGene[,-1]
  offTargetGene<-unique(offTargetGene)
  if(nrow(offTargetGene)>0){
    for(i in 1:nrow(offTargetGene)){
      if(offTargetGene$off[i]=="yes"){
        if(has.key(offTargetGene$target[i],offTargetGene_hash)){
          offTargetGene_hash[offTargetGene$target[i]]<-paste(offTargetGene_hash[[offTargetGene$target[i]]],offTargetGene$gene,sep=",")
        }else{
          offTargetGene_hash[offTargetGene$target[i]]<-offTargetGene$gene[i]
        }
      }
    }
  }
  return(offTargetGene_hash)
}
##############################   perturbation effect prioritzing   ##################################################
# *******************   calculating topics distribution for each cell
Get_distribution_diff<-function(model,sample_info_gene,KO_efficiency,index=-4.35e-05){
  require(reshape2)
  require(dplyr)
  require(entropy)
  pmatrix<-model@gamma
  topicNum<-ncol(pmatrix)
  topicName<-paste('Topic',topicNum,sep='')
  rownames(pmatrix)<-model@documents
  colnames(pmatrix)<-paste('Topic',1:topicNum,sep='')
  p.matrix<-data.frame(pmatrix,samples=rownames(pmatrix),knockout=sample_info_gene)
  p.matrix <- melt(p.matrix,id=c('samples','knockout'))
  p.step1=p.matrix %>% group_by(knockout,variable) %>% summarise(number=sum(value))
  total_number=sum(p.step1$number)
  p.step2=p.step1 %>% group_by(knockout) %>% summarise(cellNum=sum(number))
  p.step1=merge(p.step1,p.step2,by='knockout')
  p.step3<-(p.step1$number)/(p.step1$cellNum)
  p.step4<-data.frame(p.step1,ratio=p.step3)
  p.step4$ctrlNum<-p.step4[which(p.step4$knockout=="CTRL"),"cellNum"]
  p.step4$percent<-apply(p.step4[,c("cellNum","ctrlNum")],1,min)/(p.step4$ctrlNum+p.step4$cellNum)
  p.step4$background<-index*log(p.step4$percent)
  p.step4$KO_eff<-KO_efficiency[as.character(p.step4$knockout)]
  p.step4$ctrl_ratio<-p.step4[which(p.step4$knockout=="CTRL"),"ratio"]
  p.step4$diff_index<-(p.step4$ratio-p.step4$ctrl_ratio)/p.step4$ctrl_ratio
  p.step4$diff_index_adjust<-(p.step4$ratio-p.step4$ctrl_ratio)/(p.step4$ctrl_ratio*p.step4$background*p.step4$KO_eff)
  p.step4<-na.omit(p.step4)
  return(p.step4)
}
# *******************  calculating overall perturbation effect ranking list
Rank_overall<-function(distri_diff,offTarget_hash=hash()){
  require(dplyr)
  require(entropy)
  require(hash)
  KO_offTarget_hash<-hash()
  topic_number=length(distri_diff$variable[!duplicated(distri_diff$variable)])
  #diff scores get from the sum of abs(diff_index) of each topic
  decision_sum<-distri_diff %>% group_by(knockout) %>% summarise(decision_sum=round(sum(abs(diff_index_adjust)),2))
  decision_sum<-as.data.frame(decision_sum)
  decision_max<-distri_diff %>% group_by(knockout) %>% summarise(decision_max=round(max(abs(diff_index_adjust)),2))
  decision_max<-as.data.frame(decision_max)
  decision_topicName<-distri_diff[round(abs(distri_diff[,"diff_index_adjust"]),2) %in% decision_max$decision_max,c("knockout","variable")]
  colnames(decision_topicName)<-c("knockout","mainly_influenced_topic")
  KO_cellNum<-distri_diff %>% group_by(knockout) %>% summarise(KO_cellNum=sum(number))
  KO_cellNum<-as.data.frame(KO_cellNum)
  JSD<-distri_diff %>% group_by(knockout) %>% summarise(JSD=0.5*(KL.plugin(ctrl_ratio,0.5*(ratio+ctrl_ratio))+KL.plugin(0.5*(ratio+ctrl_ratio),ctrl_ratio)))
  JSD<-as.data.frame(JSD)
  background<-distri_diff %>% group_by(knockout) %>% summarise(background=unique(background))
  background<-as.data.frame(background)
  KO_eff<-distri_diff %>% group_by(knockout) %>% summarise(KO_eff=round(unique(KO_eff),3))
  KO_eff<-as.data.frame(KO_eff)
  decision<-merge(KO_eff,merge(background,merge(JSD,merge(KO_cellNum,merge(decision_max,merge(decision_sum,decision_topicName))))))
  decision$perturb_percent<-round(decision$decision_max/decision$decision_sum,4)
  decision$JSD_adjust<-decision$JSD/(decision$KO_eff*decision$background)
  #KO gene score and rank
  decision_sort<-decision[order(decision$JSD_adjust,decreasing = TRUE),]
  decision_sort$ranking=1:nrow(decision_sort)
  row.names(decision_sort)=1:nrow(decision_sort)
  decision_sort$off_target<-"none"
  for(i in 1:nrow(decision_sort)){
    ko_gene_arr<-unlist(split(decision_sort$knockout,","))
    off_target_arr<-c()
    k=1
    for(j in ko_gene_arr){
      if(has.key(j,offTarget_hash)){
        off_target_arr[k]=offTarget_hash[[j]]
        k=k+1
      }
    }
    off_target=off_target_arr[1]
    if(length(off_target_arr)>1){
      for(k in 2:length(off_target_arr)){
        off_target<-paste(off_target,off_target_arr[k],sep=",")
      }
    }
    if(!is.null(off_target)){
      KO_offTarget_hash[decision_sort$off_target[i]]=off_target
      decision_sort$off_target[i]=off_target
    }
  }
  rankOverall_result_detail<-decision_sort
  rankOverall_result_summary<-decision_sort[,c("knockout","ranking","mainly_influenced_topic")]
  return(list("rank_overall_result_summary"=rankOverall_result_summary,"rank_overall_result_detail"=rankOverall_result_detail))
}
#********************   calculating topic-specific ranking list by considering efficiency and specificity
Rank_topic_specific<-function(rank_overall_result_detail,alpha=0.5){
  rankTopicSpecific_result_detail<-rank_overall_result_detail[order(rank_overall_result_detail$mainly_influenced_topic,-rank_overall_result_detail$decision_max,decreasing = FALSE),c("mainly_influenced_topic","knockout","decision_max","perturb_percent","ranking","off_target")]
  topic_times<-table(as.character(rankTopicSpecific_result_detail$mainly_influenced_topic))
  rank<-c()
  for(i in topic_times){
    rank<-c(rank,1:i)
  }
  rankTopicSpecific_result_detail$rank<-rank
  rankTopicSpecific_result_detail$decision_max_uniform<-(rankTopicSpecific_result_detail$decision_max--min(rankTopicSpecific_result_detail$decision_max))/(max(rankTopicSpecific_result_detail$decision_max)-min(rankTopicSpecific_result_detail$decision_max))
  rankTopicSpecific_result_detail$perturb_percent_uniform<-(rankTopicSpecific_result_detail$perturb_percent--min(rankTopicSpecific_result_detail$perturb_percent))/(max(rankTopicSpecific_result_detail$perturb_percent)-min(rankTopicSpecific_result_detail$perturb_percent))
  rankTopicSpecific_result_detail$recommand_score<-alpha*rankTopicSpecific_result_detail$decision_max_uniform+(1-alpha)*rankTopicSpecific_result_detail$perturb_percent_uniform
  rankTopicSpecific_result_detail<-rankTopicSpecific_result_detail[order(rankTopicSpecific_result_detail$mainly_influenced_topic,-rankTopicSpecific_result_detail$recommand_score,decreasing=FALSE),]
  topic_times<-table(as.character(rankTopicSpecific_result_detail$mainly_influenced_topic))
  rank<-c()
  for(i in topic_times){
    rank<-c(rank,1:i)
  }
  rankTopicSpecific_result_detail$rank<-rank
  rankTopicSpecific_result_detail<-rankTopicSpecific_result_detail[,c("mainly_influenced_topic","knockout","decision_max_uniform","perturb_percent_uniform","recommand_score","ranking","off_target")]
  rankTopicSpecific_result_summary<-rankTopicSpecific_result_detail[,c("mainly_influenced_topic","knockout","ranking")]
  return(list("rank_topic_specific_result_detail"=rankTopicSpecific_result_detail,"rank_topic_specific_result_summary"=rankTopicSpecific_result_summary))
}
#********************   compare analysis if there are different condition
Rank_diff<-function(rankOverall_result_summary_condition1,rankOverall_result_summary_condition2,difference_threshold=0.5,plot=TRUE,plot_path="~/rank_diff.pdf"){
  require(ggplot2)
  time1_rankings<-as.character(rankOverall_result_summary_condition1$knockout)
  sortKnock<-sort(time1_rankings)
  ranking<-1:length(sortKnock)
  names(ranking)<-time1_rankings
  ranking_time1<-ranking[sortKnock]

  time2_rankings<-as.character(rankOverall_result_summary_condition2$knockout)
  sortKnock<-sort(time2_rankings)
  ranking<-1:length(sortKnock)
  names(ranking)<-time2_rankings
  ranking_time2<-ranking[sortKnock]

  ranking_time1_2<-ranking_time1[intersect(names(ranking_time2),names(ranking_time1))]
  ranking_time2_2<-ranking_time2[intersect(names(ranking_time2),names(ranking_time1))]
  ranking_time1_3<-1:length(ranking_time1_2)
  ex<-sort(ranking_time1_2)
  names(ranking_time1_3)<-names(ex)
  ranking_time1_ok<-ranking_time1_3[names(ranking_time2_2)]
  ranking_time2_3<-1:length(ranking_time2_2)
  ex<-sort(ranking_time2_2)
  names(ranking_time2_3)<-names(ex)
  ranking_time2_ok<-ranking_time2_3[names(ranking_time2_2)]
  ranking_difference<-ranking_time2_ok-ranking_time1_ok
  threshold_ranking<-round(length(ranking_difference)*difference_threshold)
  perturbation<-names(ranking_time1_ok)
  DataFrame<-data.frame(perturbation,ranking_difference,threshold_ranking)
  if(plot==TRUE){
    pdf(plot_path)
    p=ggplot(DataFrame,aes(x=perturbation,y=ranking_difference,group=1))+geom_line()+geom_point(size=3)+xlab("Knockout")+ylab("Ranking difference")+theme_bw()+theme(axis.text.x=element_text(angle=90,size=9))+geom_hline(aes(yintercept=threshold_ranking),col="red",linetype="dashed")+geom_hline(aes(yintercept=threshold_ranking*-1),col="red",linetype="dashed")
    print(p)
    dev.off()
  }
  difference_result<-DataFrame[which(abs(DataFrame$ranking_difference)>=threshold_ranking),]
  return(difference_result)
}
