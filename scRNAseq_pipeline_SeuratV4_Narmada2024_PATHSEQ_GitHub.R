#####################################################################################
# scRNA-seq pipeline
# Detection of Hepatitis B Virus mRNA from single cell RNA sequencing data without prior knowledge
# Nicolaas Van Renne
######################################################################################

#set working directory
	setwd("C://your_working_directory/") #adapt to your working directory
	
#load libraries
	library(Seurat)
	library(harmony)
	library(ggsci)
	library(ggplot2)
	library(patchwork)

# settings
	feat.min.cells = 3        # 3 works well;  delete features (gene) expressed in less than this number of cells (for each sample).
	set.mito = 20             # set maximal mitochondrial content
	min.nCount_RNA = 1500     # 1500 works well; minimum RNA counts per cell; cells with less than this are removed 
	max.nCount_RNA = 99999    # 99999 works well; maximum RNA counts per cell; cells with more than this are removed 
	min.nFeature_RNA = 200    # 200 generally ok; minimim amount of RNA features (i.e. genes) per cell. Will remove cells with low number of genes
	n.variablefeatures = 2000 # 2000 is ok for harmony integration 

	visualization.dimensions = 30    #  uses the first x dimensions for the dimensionality reduction plot
	PCA.dimensions = 30              #  uses the first x dimensions for PCA reduction	
	
	FindCluster.resolution = 1.0     # 2.0 for harmony is a good start for large data sets. 0.4 is a good start for small data sets. Increase for more clusters, decrease for less clusters

	seed = 18101983 # set seed for reproducibility
	set.seed(seed)

		
# load single cell data 
	#note: the directory structure for sample CHB20 should be "C://your_working_directory/scRNAseq_data/Narmada2024/CHB20" but feel free to adapt
	#Create sample name vector and data directory vector
	sample.name.vector  <- c("CHB20","CHB14","CHB13","CHB9","CHB8","CHB7","CHB5","FC9","FC6","CHB4","FC5","FC4","FC3","FC2","FC1","CHB31","CHB30","CHB28","CHB24","CHB23","CHB3","CHB2") #this will be orig.ident in meta.data
	data.dir.vector     <- c("CHB20","CHB14","CHB13","CHB9","CHB8","CHB7","CHB5","FC9","FC6","CHB4","FC5","FC4","FC3","FC2","FC1","CHB31","CHB30","CHB28","CHB24","CHB23","CHB3","CHB2") # directories where the data is to be found
	root.dir.path 		<- "scRNAseq_data/Narmada2024" #path to main directory that contains data 
	SO.list <- list()
		
	#create seurat object list with selected samples 
		for(i in seq_along(sample.name.vector)){
		  print(paste("processing sample",i,":",sample.name.vector[i]))
		  SO.list[[i]] <- Read10X(data.dir = eval(paste0(root.dir.path,"/",data.dir.vector[i])))
		  SO.list[[i]] <- CreateSeuratObject(counts = SO.list[[i]], project = sample.name.vector[i], min.cells = feat.min.cells, min.features = min.nFeature_RNA) 
		  SO.list[[i]] <- RenameCells(SO.list[[i]], add.cell.id = data.dir.vector[i])
		}

	#pre-process PATHSEQ reads
		key <- rbind(
			c("sample001","SRR28562212","CHB20") ,
			c("sample002","SRR28562213","CHB14") ,
			c("sample003","SRR28562214","CHB13") ,
			c("sample004","SRR28562215","CHB9" ) ,
			c("sample005","SRR28562216","CHB8" ) ,
			c("sample006","SRR28562227","CHB7" ) ,
			c("sample007","SRR28562238","CHB5" ) ,
			c("sample008","SRR28562246","FC9"  ) ,
			c("sample009","SRR28562248","FC6"  ) ,
			c("sample010","SRR28562249","CHB4" ) ,
			c("sample011","SRR28562250","FC5"  ) ,
			c("sample012","SRR28562251","FC4"  ) ,
			c("sample013","SRR28562252","FC3"  ) ,
			c("sample014","SRR28562253","FC2"  ) ,
			c("sample015","SRR28562254","FC1"  ) ,
			c("sample016","SRR28562255","CHB31") ,
			c("sample017","SRR28562256","CHB30") ,
			c("sample018","SRR28562257","CHB28") ,
			c("sample019","SRR28562258","CHB24") ,
			c("sample020","SRR28562259","CHB23") ,
			c("sample021","SRR28562260","CHB3" ) ,
			c("sample022","SRR28562261","CHB2" ) ,
			c("MISSING_SRA","SRR28562247","FC8")
		)
			
			
		#load csv (output of pathseq)
			umi_table_filename <- list()
			umi_table_filename[[01]] = "pathseq_output/sample001.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[01] <- "CHB20"
			umi_table_filename[[02]] = "pathseq_output/sample002.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[02] <- "CHB14"
			umi_table_filename[[03]] = "pathseq_output/sample003.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[03] <- "CHB13"
			umi_table_filename[[04]] = "pathseq_output/sample004.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[04] <- "CHB9"
			umi_table_filename[[05]] = "pathseq_output/sample005.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[05] <- "CHB8"
			umi_table_filename[[06]] = "pathseq_output/sample006.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[06] <- "CHB7"
			umi_table_filename[[07]] = "pathseq_output/sample007.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[07] <- "CHB5"
			umi_table_filename[[08]] = "pathseq_output/sample008.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[08] <- "FC9"
			umi_table_filename[[09]] = "pathseq_output/sample009.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[09] <- "FC6"
			umi_table_filename[[10]] = "pathseq_output/sample010.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[10] <- "CHB4"
			umi_table_filename[[11]] = "pathseq_output/sample011.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[11] <- "FC5"
			umi_table_filename[[12]] = "pathseq_output/sample012.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[12] <- "FC4"
			umi_table_filename[[13]] = "pathseq_output/sample013.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[13] <- "FC3"
			umi_table_filename[[14]] = "pathseq_output/sample014.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[14] <- "FC2"
			umi_table_filename[[15]] = "pathseq_output/sample015.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[15] <- "FC1"
			umi_table_filename[[16]] = "pathseq_output/sample016.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[16] <- "CHB31"
			umi_table_filename[[17]] = "pathseq_output/sample017.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[17] <- "CHB30"
			umi_table_filename[[18]] = "pathseq_output/sample018.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[18] <- "CHB28"
			umi_table_filename[[19]] = "pathseq_output/sample019.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[19] <- "CHB24"
			umi_table_filename[[20]] = "pathseq_output/sample020.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[20] <- "CHB23"
			umi_table_filename[[21]] = "pathseq_output/sample021.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[21] <- "CHB3"
			umi_table_filename[[22]] = "pathseq_output/sample022.gex.filtered_matrix.genus.csv"  ; names(umi_table_filename)[22] <- "CHB2"
		
			umi_table <- list()
			total_list <- list()
			for(i in 1:length(umi_table_filename)){
				umi_table[[i]] <- read.csv(umi_table_filename[[i]],sep=',',header=TRUE,row.names = 1)
				total_list[[i]] <- rowSums(umi_table[[i]])
			}
			#umi_table[umi_table==0] <- NA
	
		#adapt umi_table
			#change barcodes so they match the seurat object
				for(i in 1:length(umi_table_filename)){
					print(paste("should refer to same sample:", key[i,1], key[i,3], sample.name.vector[i], umi_table_filename[[i]]))
					rownames(umi_table[[i]]) <- gsub(key[i,1],key[i,3],rownames(umi_table[[i]]))
				}	
				
			#make sure the separate matrices contain the same column names (microbes) in the same order	
				all_columns <- unique(unlist(lapply(umi_table, colnames)))
		
				umi_list <- lapply(umi_table, function(mat) {
				missing_cols <- setdiff(all_columns, colnames(mat))
				if (length(missing_cols) > 0) {
					mat <- cbind(mat, matrix(0, nrow = nrow(mat), ncol = length(missing_cols)))
					colnames(mat)[(ncol(mat) - length(missing_cols) + 1):ncol(mat)] <- missing_cols
				}
				mat <- mat[, all_columns]  # Reorder columns to ensure consistency
				return(mat)
				})
		
				# QC: Check if all matrices in umi_list have the same column names in the same order
					are_colnames_identical <- all(sapply(umi_list, function(mat) {
					identical(colnames(mat), colnames(umi_list[[1]]))
					}))
					
					if (are_colnames_identical) {
					message("OK: All matrices have the same column names in the same order.")
					} else {
					message("error: Column names differ across matrices.")
					}
			
			#merge matrices into a single matrix
				merged_matrix <- do.call(rbind, umi_list)		
				all_microbes <- colnames(merged_matrix)
				
	#add the adapted matrix to the raw readcounts in the seurat object
		updated_counts_matrix <- list()
		for(i in 1:length(sample.name.vector)){
			print(paste("add microbe reads to sample",i,":",sample.name.vector[i]))
			#transpose umi_table
				transposed_umi_table <- t(umi_list[[i]])
			
			# Retrieve the read count matrix from the Seurat object
				counts_matrix <- SO.list[[i]][["RNA"]]@counts
			
			# Identify missing cells in the transposed UMI table
				missing_cells <- setdiff(colnames(counts_matrix), colnames(transposed_umi_table))
			
			# Add missing cells (columns) with zeros to the transposed UMI table
				if (length(missing_cells) > 0) {
				zero_matrix_cells <- matrix(0, nrow = nrow(transposed_umi_table), ncol = length(missing_cells))
				colnames(zero_matrix_cells) <- missing_cells
				transposed_umi_table <- cbind(transposed_umi_table, zero_matrix_cells)
				}

			# Add missing cells (columns) with zeros to counts_matrix 
				missing_columns_in_counts <- setdiff(colnames(transposed_umi_table), colnames(counts_matrix))
				if(length(missing_columns_in_counts) > 0){
				zero_matrix <- matrix(0, nrow = nrow(counts_matrix), ncol = length(missing_columns_in_counts))
				colnames(zero_matrix) <- missing_columns_in_counts
				counts_matrix <- cbind(counts_matrix, zero_matrix)
				}
			
			# ensure the columns in both matrices are in the same order
				counts_matrix <- counts_matrix[, colnames(transposed_umi_table)]
			
			#QC: check if colnames (cell barcodes) match
				if(!all(colnames(transposed_umi_table)==colnames(counts_matrix))){print("ERROR: barcodes do not match")}

			#add microbe count matrix to seurat count matrix 
				updated_counts_matrix[[i]] <- rbind(counts_matrix, transposed_umi_table)
		}


	#create seurat object list with selected samples 
		for(i in seq_along(sample.name.vector)){
		  print(paste("processing sample",i,":",sample.name.vector[i]))
		  SO.list[[i]] <- CreateSeuratObject(counts = updated_counts_matrix[[i]], project = sample.name.vector[i], min.cells = feat.min.cells, min.features = min.nFeature_RNA) 
		  SO.list[[i]][["percent.mt"]]    <- PercentageFeatureSet(SO.list[[i]], pattern = "^MT-")
		  SO.list[[i]] <- subset(SO.list[[i]], subset = nFeature_RNA > min.nFeature_RNA & percent.mt < set.mito & nCount_RNA > min.nCount_RNA & nCount_RNA < max.nCount_RNA)			
		}
	
				
	#add meta.data (you can add as many meta data as you want.
		SO.list[[1]][["SLoss"]] <- "NoSLoss" 	#	NonSLoss	CHB20
		SO.list[[2]][["SLoss"]] <- "NoSLoss"	#	NonSLoss	CHB14
		SO.list[[3]][["SLoss"]] <- "NoSLoss"	#	NonSLoss	CHB13
		SO.list[[4]][["SLoss"]] <- "NoSLoss"	#	NonSLoss	CHB9
		SO.list[[5]][["SLoss"]] <- "NoSLoss"	#	NonSLoss	CHB8
		SO.list[[6]][["SLoss"]] <- "NoSLoss"	#	NonSLoss	CHB7
		SO.list[[7]][["SLoss"]] <- "NoSLoss"	#	NonSLoss	CHB5
		SO.list[[8]][["SLoss"]] <- "SLoss"		#	Sloss	FC9
		SO.list[[9]][["SLoss"]] <- "SLoss"		#	Sloss	FC6
		SO.list[[10]][["SLoss"]] <- "NoSLoss"	#	NonSLoss	CHB4
		SO.list[[11]][["SLoss"]] <- "SLoss"		#	Sloss	FC5
		SO.list[[12]][["SLoss"]] <- "SLoss"		#	Sloss	FC4
		SO.list[[13]][["SLoss"]] <- "SLoss"		#	Sloss	FC3
		SO.list[[14]][["SLoss"]] <- "SLoss"		#	Sloss	FC2
		SO.list[[15]][["SLoss"]] <- "SLoss"		#	Sloss	FC1
		SO.list[[16]][["SLoss"]] <- "NoSLoss"	#	NonSLoss	CHB31
		SO.list[[17]][["SLoss"]] <- "NoSLoss"	#	NonSLoss	CHB30
		SO.list[[18]][["SLoss"]] <- "NoSLoss"	#	NonSLoss	CHB28
		SO.list[[19]][["SLoss"]] <- "NoSLoss"	#	NonSLoss	CHB24
		SO.list[[20]][["SLoss"]] <- "NoSLoss"	#	NonSLoss	CHB23
		SO.list[[21]][["SLoss"]] <- "NoSLoss"	#	NonSLoss	CHB3
		SO.list[[22]][["SLoss"]] <- "NoSLoss"	#	NonSLoss	CHB2
		

	#integrate using Harmony 
		#if one sample, run the non-integrated pipeline	
		if(length(SO.list)==1){		
		  SO.integrated <- SO.list[[1]]
		  rm(SO.list)
		  SO.integrated <- NormalizeData(SO.integrated)
		  SO.integrated <- FindVariableFeatures(SO.integrated, nfeatures = n.variablefeatures)
		  SO.integrated <- ScaleData(SO.integrated)
		  SO.integrated <- RunPCA(SO.integrated)
		  SO.integrated <- RunUMAP(SO.integrated, dims = 1:visualization.dimensions)
		  SO.integrated <- FindNeighbors(SO.integrated, dims = 1:PCA.dimensions) 
		  SO.integrated <- FindClusters(SO.integrated, resolution = FindCluster.resolution) 
		}else{
		  #if more than one sample integrate using Harmony 
		  SO.integrated <- merge(SO.list[[1]], y = SO.list[-1])
		  rm(SO.list)
		  SO.integrated <- NormalizeData(SO.integrated)
		  SO.integrated <- FindVariableFeatures(SO.integrated, nfeatures = n.variablefeatures)
		  SO.integrated <- ScaleData(SO.integrated)
		  SO.integrated <- RunPCA(SO.integrated)
		  SO.integrated <- RunHarmony(SO.integrated, group.by.vars = "orig.ident")
		  SO.integrated <- RunUMAP(SO.integrated, reduction = "harmony", dims = 1:visualization.dimensions)
		  SO.integrated <- FindNeighbors(SO.integrated, reduction = "harmony", dims = 1:PCA.dimensions) 
		  SO.integrated <- FindClusters(SO.integrated, resolution = FindCluster.resolution) 
		}
	
##############################################################
# Find markers for annotating cell clusters
##############################################################
	# settings
	min.pct.expression = 0.25 #standard setting: 0.25
	min.logfc = 0.25 #0.25 is standard
	p.val.cutoff <- (1/10^3) #(1/10^3) is standard, use (1/10^0) to ignore
	DefaultAssay(SO.integrated) <- "RNA" # "RNA" or "SCSIGN" is also possible if single cell signatures are calculated, or "ADT" for citeSeq
	
	# run algorithm
	library(future)
	plan("multicore", workers = 4)
	
	cluster.names <- unique(Idents(SO.integrated))[order(unique(Idents(SO.integrated)))];
	markers.fm.list <- list()
	for (i in 1:length(cluster.names)) {
	  print(paste0("calculating markers for cluster ",cluster.names[i],". Total: ",length(cluster.names)," clusters"))
	  markers.fm.list[[i]] <- FindMarkers(SO.integrated, ident.1 = cluster.names[i], min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)
	  markers.fm.list[[i]] <- markers.fm.list[[i]][which(markers.fm.list[[i]]$p_val_adj < (p.val.cutoff) ),] #select genes with p_val_adj > p.val.cutoff setting
	}
	
	#visualise
	for (i in 1:length(markers.fm.list)){
	  print(DotPlot(SO.integrated, cols = "RdYlBu",features = rev( rownames(markers.fm.list[[i]])[1:40]) ) + RotatedAxis() )
	  print(paste0(cluster.names[i]))
	  readline(prompt="Press [enter] to continue")
	}	
	
#############################################################################################
# Manually annotate clusters
# 
# Here you can computationally isolate clusters, explore their identity and  
# 	manually save the annotations, and then inject them back in the original seurat object
#
#############################################################################################
	# set annotations 
		annotations.list <- list()
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =0 )  ;names(annotations.list)[length(annotations.list)] <-   "T cytotox"       #celltype 0     
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =1 )  ;names(annotations.list)[length(annotations.list)] <-   "MAIT/ILC (1)"     
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =2 )  ;names(annotations.list)[length(annotations.list)] <-   "NK CD56bright"     
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =3 )  ;names(annotations.list)[length(annotations.list)] <-   "NK CD56dim"     
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =4 )  ;names(annotations.list)[length(annotations.list)] <-   "Tn/Tcm"     
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =5 )  ;names(annotations.list)[length(annotations.list)] <-   "Texh"  #celltype 5         
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =6 )  ;names(annotations.list)[length(annotations.list)] <-   "neutrophil"     
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =7 )  ;names(annotations.list)[length(annotations.list)] <-   "endothelial (1)"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =8 )  ;names(annotations.list)[length(annotations.list)] <-   "mono CL"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =9 )  ;names(annotations.list)[length(annotations.list)] <-   "B cell"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =10)  ;names(annotations.list)[length(annotations.list)] <-   "cDC2 (1)"   #celltype 10
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =11)  ;names(annotations.list)[length(annotations.list)] <-   "hepatocyte"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =12)  ;names(annotations.list)[length(annotations.list)] <-   "plasmocyte"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =13)  ;names(annotations.list)[length(annotations.list)] <-   "T mito hi"                 
		#annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =14)  ;names(annotations.list)[length(annotations.list)] <-   "KC (2)"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =15)  ;names(annotations.list)[length(annotations.list)] <-   "mono NCL"  #celltype 15               
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =16)  ;names(annotations.list)[length(annotations.list)] <-   "endothelial (2)"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =17)  ;names(annotations.list)[length(annotations.list)] <-   "cholangiocyte"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =18)  ;names(annotations.list)[length(annotations.list)] <-   "cycling"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =19)  ;names(annotations.list)[length(annotations.list)] <-   "hep/T/NK doublet"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =20)  ;names(annotations.list)[length(annotations.list)] <-   "KC (1)"      #celltype 20                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =21)  ;names(annotations.list)[length(annotations.list)] <-   "pDC"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =22)  ;names(annotations.list)[length(annotations.list)] <-   "cDC1 (1)"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =23)  ;names(annotations.list)[length(annotations.list)] <-   "endothelial (3)"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =24)  ;names(annotations.list)[length(annotations.list)] <-   "stromal"                 
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =25)  ;names(annotations.list)[length(annotations.list)] <-   "MAIT/ILC (2)"    #celltype 25                    
		annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =26)  ;names(annotations.list)[length(annotations.list)] <-   "mast cell"                 

		
		#partition cluster KC (2) low quality cells 
			SO.backup <- SO.integrated 
			SO.integrated <- subset(SO.integrated, idents=14)
		
			FindCluster.resolution = 0.2
			SO.integrated <- FindNeighbors(SO.integrated, reduction = "pca", dims = 1:30)
			SO.integrated <- FindClusters(SO.integrated, resolution = FindCluster.resolution) #hoe hoger res, hoe meer clusters
				
				DimPlot(SO.integrated)
				dotplot.markers <- c("XCR1","CLEC9A","CD1C","CD1D","CD1E","CLEC10A","MARCO","CD5L","SLC40A1","HMOX1")
				p=DotPlot(SO.integrated, cols = "RdYlBu",features = dotplot.markers) + RotatedAxis() & theme(axis.title.x = element_blank(),axis.title.y = element_blank())
				p
			
			#add annotations of the cluster KC (1) ribo high
			annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =0 )  ;names(annotations.list)[length(annotations.list)] <-   "KC (2)"       #celltype 0     
			annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =1 )  ;names(annotations.list)[length(annotations.list)] <-   "cDC2 (2)"     
			annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =2 )  ;names(annotations.list)[length(annotations.list)] <-   "KC (2)"     
			annotations.list[[length(annotations.list)+1]] <- WhichCells(SO.integrated, idents =3 )  ;names(annotations.list)[length(annotations.list)] <-   "cDC1 (2)"     
	
		#reset seurat object
			SO.integrated <- SO.backup

	# Inject new annotations in Seurat object
		for(i in seq_along(annotations.list)){
		  Idents(SO.integrated, cells = annotations.list[[i]])  <- names(annotations.list)[[i]]
		}
		
	# save in meta.data$annotation
		SO.integrated[["annotation"]] <- Idents(SO.integrated)
		DimPlot(SO.integrated, label=T, label.size=5)

	# reorder idents and create DimPlots 
		Idents(SO.integrated) <- factor(Idents(SO.integrated), levels = c(
			"hepatocyte"               ,
			"hep/T/NK doublet"  	   ,
			"cholangiocyte"  		   ,	
			"endothelial (1)"          ,
			"endothelial (2)"		   ,
			"endothelial (3)" 		   ,
			"stromal"  		           ,
			"Tn/Tcm"   		           ,
			"T cytotox"                ,
			"Texh"                     ,
			"T mito hi"                ,				
			"MAIT/ILC (1)"             ,
			"MAIT/ILC (2)"    		   ,
			"NK CD56bright"            ,
			"NK CD56dim"               ,
			"B cell"                   ,
			"plasmocyte" 		       ,
			"mono CL"                  ,
			"mono NCL"                 ,
			"cDC1 (1)"                 ,
			"cDC1 (2)"                 ,
			"cDC2 (1)"                 ,
			"cDC2 (2)"                 ,
			"KC (1)"                   ,
			"KC (2)"                   ,
			"pDC"               	   ,
			"neutrophil"               ,
			"mast cell"                ,   
			"cycling"                  )
		)
		DimPlot(SO.integrated, label=T, label.size=5)
		
###########################
# Reproduce all plots
###########################
	#create color vector
		library(ggsci)
			n.clusters <- length(unique(Idents(SO.integrated)))
			if(n.clusters <= 10){
			color.vector <- pal_npg("nrc")(n.clusters)
			}else{
			#retrieve color codes
			ggplotColours <- function(n = 6, h = c(0, 360) + 15){
				if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
				hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
			}
			color.vector <- c(pal_npg("nrc")(10),ggplotColours(n= (n.clusters-10)))
			}
			color.vector.alpha2 <- alpha(color.vector,0.2)	
			color.vector.alpha5 <- alpha(color.vector,0.5)	 
		
	#UMAP representation plots
		p=DimPlot(SO.integrated, cols= color.vector, pt.size=1)
		p

		p=DimPlot(SO.integrated, cols= color.vector.alpha5, pt.size=0.5) + theme_void() + theme(legend.position = 'none')
		p
		

	#DotPlots
		dotplot.markers <- c("TTR","CYP2E1","SOX9","KRT19","VWF","MECOM","CLEC1B","CLEC4G","CD34","CXCL12","MYL9","COL1A1","CD2","CD3D","CCR7","TCF7","LTB","CD8A","CD8B","TIGIT","PDCD1","CTLA4","IL23R","SLC4A10","CCR6","NCAM1","XCL1","XCL2","CD160","GNLY","GZMB","MS4A1","BANK1","MZB1","TNFRSF17","CD14","VCAN","FCGR3A","XCR1","CLEC9A","CD1C","CLEC10A","MARCO","CD5L","CLEC4C","SCT","G0S2","CSF3R","CPA3","KIT","MKI67","TOP2A")
		p=DotPlot(SO.integrated, cols = "RdYlBu",features = dotplot.markers) + RotatedAxis() & theme(axis.title.x = element_blank(),axis.title.y = element_blank())
		p

	#HBV plots
		p=VlnPlot(SO.integrated, "Orthohepadnavirus", split.by="SLoss")
		p
		 
		p=DotPlot(SO.integrated, cols = "Spectral",features = all_microbes) + RotatedAxis() & theme(axis.title.x = element_blank(),axis.title.y = element_blank())
		p

		#generate plots in patients without SLoss
			SO.NoSLoss <- subset(SO.integrated, SLoss=="NoSLoss")
			p=DotPlot(SO.NoSLoss,features = 'Orthohepadnavirus', cols="Spectral") + RotatedAxis() & theme(axis.title.x = element_blank(),axis.title.y = element_blank())
			p
				
			p=FeaturePlot(SO.NoSLoss, features = 'Orthohepadnavirus', order=T) &   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
			p


	#Cluster proportion table
		clustertable <- table(Idents(SO.integrated), SO.integrated@meta.data$orig.ident) #color by dataset (orig.ident)
		plot.orig.ident <- ggplot(as.data.frame(clustertable), aes(fill=Var2, y=Freq, x=Var1)) + 
		  geom_bar(position="fill", stat="identity")
		
		denominator <-  as.numeric(rowSums(clustertable))	
		clustertable.normalized <- clustertable
		for(i in 1:nrow(clustertable)){
		  clustertable.normalized[i,] <- clustertable[i,]/denominator[i]
		}
		
		p=ggplot(as.data.frame(clustertable.normalized), aes(fill=Var2, y=Freq, x=Var1)) + 
		  geom_bar(position="fill", stat="identity") + 
		  theme_classic() +
		  xlab("") +
		  ylab("relative contribution") +
		  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.25, size = 20, color="black")) +
		  theme(axis.text.y = element_text(size = 20, color="black")) +
		  theme(axis.title.y = element_text(size = 20)) +
		  scale_y_continuous(expand = c(0, 0), limits=c(0,1.05), breaks=c(0,0.25,0.5,0.75,1)) +
		  theme(legend.text=element_text(size=20)) +
		  theme(legend.title=element_text(size=20)) +
		  guides(fill=guide_legend("donor sample"))
		p		
	
	#low quality myeloid cells ridgeplots
		SO.DCKC <- subset(SO.integrated, idents=c("cDC1 (1)","cDC1 (2)","cDC2 (1)","cDC2 (2)","KC (1)","KC (2)"))
		Idents(SO.DCKC) <- factor(Idents(SO.DCKC), levels = c("cDC1 (1)","cDC2 (1)","KC (1)","cDC1 (2)","cDC2 (2)","KC (2)"))
			
		p=RidgePlot(SO.DCKC, features = c("nFeature_RNA"), cols=color.vector[c(20,22,24,21,23,25)])+ theme_classic() +	guides(fill = guide_legend(reverse = TRUE))
		p
		
		p=RidgePlot(SO.DCKC, features = c("nCount_RNA"),cols=color.vector[c(20,22,24,21,23,25)])+theme_classic()
		p

#session info
 sessionInfo()
		#	R version 4.0.3 (2020-10-10)
		#	Platform: x86_64-w64-mingw32/x64 (64-bit)
		#	Running under: Windows 10 x64 (build 19045)
		#	
		#	Matrix products: default
		#	
		#	locale:
		#	[1] LC_COLLATE=Dutch_Belgium.1252  LC_CTYPE=Dutch_Belgium.1252    LC_MONETARY=Dutch_Belgium.1252 LC_NUMERIC=C                   LC_TIME=Dutch_Belgium.1252    
		#	
		#	attached base packages:
		#	[1] stats     graphics  grDevices utils     datasets  methods   base     
		#	
		#	other attached packages:
		#	[1] patchwork_1.1.1    ggplot2_3.4.4      ggsci_2.9          harmony_0.1.0      Rcpp_1.0.7         SeuratObject_4.0.4 Seurat_4.0.2      
		#	
		#	loaded via a namespace (and not attached):
		#	[1] nlme_3.1-149          matrixStats_0.58.0    spatstat.sparse_3.0-2 RcppAnnoy_0.0.18      RColorBrewer_1.1-2    httr_1.4.2            sctransform_0.3.2     tools_4.0.3           utf8_1.2.1           
		#	[10] R6_2.5.0              irlba_2.3.3           rpart_4.1-15          KernSmooth_2.23-17    uwot_0.1.10           mgcv_1.8-33           lazyeval_0.2.2        colorspace_2.0-0      withr_2.5.0          
		#	[19] tidyselect_1.2.0      gridExtra_2.3         compiler_4.0.3        cli_3.6.2             plotly_4.9.3          labeling_0.4.2        scales_1.2.1          lmtest_0.9-38         spatstat.data_3.0-1  
		#	[28] ggridges_0.5.3        pbapply_1.4-3         goftest_1.2-2         stringr_1.5.1         digest_0.6.27         spatstat.utils_3.0-3  pkgconfig_2.0.3       htmltools_0.5.2       parallelly_1.32.0    
		#	[37] fastmap_1.1.0         htmlwidgets_1.5.3     rlang_1.1.3           shiny_1.6.0           farver_2.1.0          generics_0.1.3        zoo_1.8-12            jsonlite_1.7.2        ica_1.0-2            
		#	[46] dplyr_1.1.1           magrittr_2.0.1        Matrix_1.3-4          munsell_0.5.0         fansi_0.4.2           abind_1.4-5           reticulate_1.18       lifecycle_1.0.3       stringi_1.5.3        
		#	[55] MASS_7.3-58.3         Rtsne_0.15            plyr_1.8.6            grid_4.0.3            parallel_4.0.3        listenv_0.8.0         promises_1.2.0.1      ggrepel_0.9.1         crayon_1.4.1         
		#	[64] miniUI_0.1.1.1        deldir_1.0-6          lattice_0.20-41       cowplot_1.1.1         splines_4.0.3         tensor_1.5            pillar_1.9.0          igraph_1.2.6          spatstat.geom_3.2-1  
		#	[73] future.apply_1.7.0    reshape2_1.4.4        codetools_0.2-16      leiden_0.3.7          glue_1.6.2            data.table_1.14.0     png_0.1-7             vctrs_0.6.5           httpuv_1.5.5         
		#	[82] gtable_0.3.0          RANN_2.6.1            purrr_1.0.2           spatstat.core_2.0-0   polyclip_1.10-0       tidyr_1.3.1           scattermore_0.7       future_1.26.1         mime_0.10            
		#	[91] xtable_1.8-4          RSpectra_0.16-0       later_1.1.0.1         survival_3.2-7        viridisLite_0.4.2     tibble_3.2.1          cluster_2.1.0         globals_0.15.1        fitdistrplus_1.1-3   
		#	[100] ellipsis_0.3.2        ROCR_1.0-11          

