#takes a directory location where you want to save a file, a list of plot objects to save and optional info on the quality of the life to be made
batch_plot_save<-function(directory,obj_list,name.prefix,width=1020,height=1020,units="px",type="png"){
  
  for(p in 1:length(obj_list)){
    switch(type,
           'png'={png(filename = paste0(directory,"/",p,".png"),width=width,height=height,units = units)
             plot(obj_list[[p]])
             dev.off()
           },
           'tiff'={
             ggsave(plot=obj_list[[p]],filename = paste0(directory,"/",name.prefix,p,".tiff"),device = "tiff",width = width,height = height,units=units)
           },
           'emf'={
             emf(file = paste0(directory,"/",p,".emf"),width=width,height=height)
             plot(obj_list[[p]])
             dev.off()
           },
           "eps"={
             ggsave(plot=obj_list[[p]],filename = paste0(directory,"/",name.prefix,obj_list[[p]]$labels$title,".eps"),device = "eps",width = width,height = height,units=units)
           },
           
    )
  } 
  
}


#takes a seurat object that has two or more integrated objects and returns a graph of the relative ratios of object a to object b in each cluster 
#grouping is required for integrated objects with more than two component objects 
#grouping should be a vector of two vectors that are the names of the orig.ident entries that should be combined for the graph 
integration_graph<-function(obj,data_source="orig.ident",values="proportion",label_size=0.25,min.segment.length=0.1,label=TRUE){
  comp_obj<-table(obj@meta.data[[data_source]])#get the names and frequency of the component objects of the integrated object
  if(length(comp_obj)==1){
    stop("object only has one annotated original identity in metadata")
  }
  object_1<-names(comp_obj)[1]
  object_2<-names(comp_obj)[2]
  ratio<-matrix(nrow=0,ncol=3)
  colnames(ratio)<-c(object_1,object_2,"cluster_name")
  for (cluster in levels(obj@active.ident)) {#pull out each cluster and calculate seperately
    obj_cluster<-subset(obj,idents=cluster)#TODO rewrite to use a fetch data request to get the needed subset so the function runs faster
    num_cells<-table(obj_cluster@meta.data[[data_source]])#pull out cell numbers from meta-data
    switch(values,
           "proportion"={
             ratio<-rbind(ratio,list((num_cells[object_1][[1]]/comp_obj[[1]]),(num_cells[object_2][[1]]/comp_obj[[2]]),as.character(cluster)))
             reference_line<-geom_abline(slope = 1 ,intercept = 0,color='green')
           },#make the table row, the proportion of values from each sample
           
           "number"={
             ratio<-rbind(ratio,list(num_cells[object_1][[1]],num_cells[object_2][[1]],as.character(cluster)))
             reference_line<-geom_abline(slope = comp_obj[[2]]/comp_obj[[1]],intercept = 0,color='green')
           }
    )
  }
  ratio<-data.frame(ratio)
  ratio[,1]<-unlist(ratio[,1])
  ratio[,2]<-unlist(ratio[,2])
  ratio$cluster_name<-unlist(ratio$cluster_name)
  rownames(ratio)<-ratio$cluster_name
  #now plot, currently dots maybe some kind of radial graph would be better given that the determinant factor for exclusion is slope
  #green line is the expected slope of an "evenly divided cluster based on the total number of cells in each group
  coexpress_plot<-ggplot(ratio,aes_string(x=object_1,y=object_2,color=factor(ratio$cluster_name,levels = ratio$cluster_name)))+geom_point()+reference_line# +ggtitle(paste(names(comp_obj[2]),"vs",names(comp_obj[1])))  
  if(label){
    coexpress_plot=coexpress_plot+geom_label_repel(label = rownames(ratio),label.size=label_size,min.segment.length = min.segment.length,max.overlaps = 30)
  }
  return(coexpress_plot)
}

###takes a seurat object and returns a plot of the ratios between source identities for each cluster
#plot has clusters on the x axis and the y axis is the ratio of one source identity to the other 
integrated_cluster_ratio_graph<-function(obj,data_source="orig.ident",reverse=FALSE){
  comp_obj<-table(obj@meta.data[[data_source]])#get the names and frequency of the component objects of the integrated object
  if(length(comp_obj)==1){
    stop("object only has one annotated original identity in metadata")
  }
  object_1<-names(comp_obj)[1]
  object_2<-names(comp_obj)[2]
  if(reverse){#lets the user reverse the proportion of the graph
    object_1<-names(comp_obj)[2]
    object_2<-names(comp_obj)[1] 
  }
  ratio<-matrix(nrow=0,ncol=2)
  colnames(ratio)<-c("proportion","cluster_name")
  for (cluster in levels(obj@active.ident)) {#pull out each cluster and calculate seperately
    obj_cluster<-subset(obj,idents=cluster)#this would be much faster if it wasn't subsetting
    num_cells<-table(obj_cluster@meta.data[[data_source]])#pull out cell numbers from meta-data
    ratio<-rbind(ratio,list((num_cells[object_2][[1]]/num_cells[object_1][[1]]),as.character(cluster)))#make the table row, the proportion of values from each sample
  }
  
  ratio<-as.data.frame(ratio)
  ratio[,1]<-unlist(ratio[,1])
  ratio$cluster_name<-unlist(ratio$cluster_name)#still not clear to me why ggplot doesn't like these when they are not unlisted
  rownames(ratio)<-ratio$cluster_name
  coexpress_plot<-ggplot(ratio,aes_string(x=factor(ratio$cluster_name,levels = ratio$cluster_name),y="proportion",color=factor(ratio$cluster_name,levels = ratio$cluster_name)))+geom_point()+geom_hline(yintercept = 1,color='green')+
    ggtitle(paste(object_2,"/",object_1))+ylim(0,NA)+ylab(paste("proportion",object_2,'vs',object_1))#color=scale_color_manual(values = cluster_name)
  return(coexpress_plot)
}


coexpress_bar<-function(object,idents,genes){
  #check that the given idents and genes are in the object
  data<-data.frame()
  for(ident in idents){
    ident.cells<-WhichCells(object = object,idents =ident)#gets the total cells for this object, length is the number of cells in the onject and the list of cells is useeful for fetchdata later    
    gene.ident.cells<-FetchData(object = object,vars=genes,cells=ident.cells)#fetch the expression data for the requested genes 
    gene.colocalize<-gene.ident.cells[apply(gene.ident.cells,1,function(row) all(row !=0)),]#takes the subset of the cells in the object with no zero values in the pulled genes
    gene.percent<-dim(gene.colocalize)[1]/dim(gene.ident.cells)[1]*100# get the percent of cells with requested genes
    data<-rbind(data,data.frame(cells=gene.percent,cluster=ident)) 
    
    coexpress.plot<-ggplot(data = data,aes(x=factor(cluster,levels = idents),y=cells))+
      geom_bar(stat = "identity",position = "dodge",color="black",width=0.5)+
      scale_fill_manual(values=c("#3A3B3C","grey"))+ggtitle("Coexpression")+
      theme_cowplot()+
      theme(axis.title.x = element_blank(),title = element_text(hjust = 0.5))+
      ylab("% Of Total Cells")
  }
  return(coexpress.plot)
}

quality_frame<-function(object){
  #takes a seurat object and returns a table of the average, and median, ncounts, n features and mitochondrial precentage for each cluster
  identities<-levels(object)
  output_frame<-data.frame('identity','# of cells','mean_count','median_count','mean_feature','median_feature','mean_mitochondrial%','median_mitochondrial%')
  for(i in identities){
    temp<-subset(object,cells = WhichCells(object,idents=i))
    new_row<-c(i,dim(temp@meta.data)[1],round(mean(temp@meta.data[["nCount_RNA"]])),round(median(temp@meta.data[["nCount_RNA"]])),round(mean(temp@meta.data[["nFeature_RNA"]])),round(median(temp@meta.data[["nFeature_RNA"]])),round(mean(temp@meta.data[["percent_mito"]]),digits = 2),round(median(temp@meta.data[["percent_mito"]]),digits = 2))
    output_frame<-rbind(output_frame,new_row)
  }
  return(output_frame)
}