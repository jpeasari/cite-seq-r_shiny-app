library(shiny)
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(writexl)
library(readxl)
library(data.table)
library(magrittr)
library(tidyverse)


dirs_group01 <- list.dirs(path = 'data/Group01/', recursive = F, full.names = F)
dirs_group02 <- list.dirs(path = 'data/Group02/', recursive = F, full.names = F)


skus <- list('Group1' = as.list(dirs_group01),'Group2' = as.list(dirs_group02))

sample_names <- list.dirs(path = 'all_data/', recursive = F, full.names = F)



###################################################
testing <- function(sample_name) {
  dir <- getwd()
  setwd(dir)
  
  sample_name <- "dana"
  
  sample <- paste("output",sample_name,sep = "/")
  
  sample <- as.character(sample)
  
  
  ## Creating a new directory
  dir_1 <- paste(sample,"Output",sep="_")
  dir_2 <- paste(sample,"Plots",sep="_")
  dir_3 <- paste(sample,"RData",sep="_")
  
  
  if (!dir.exists(dir_1)){
    dir.create(dir_1)
  }else{
    print("dir exists")
  }
  
  if (!dir.exists(dir_2)){
    dir.create(dir_2)
  }else{
    print("dir exists")
  }
  
  if (!dir.exists(dir_3)){
    dir.create(dir_3)
  }else{
    print("dir exists")
  }
}

######################
run_citeseq <- function(dirs_group01,dirs_group02)  {
print("Started CITE SEQ Analysis...................................")
  dir <- getwd()
  setwd(dir)
  
  sample_name <- "dana"
  
  sample <- paste("output",sample_name,sep = "/")
  
  sample <- as.character(sample)
  
  
  ## Creating a new directory
  dir_1 <- paste(sample,"Output",sep="_")
  dir_2 <- paste(sample,"Plots",sep="_")
  dir_3 <- paste(sample,"RData",sep="_")
  
  
  if (!dir.exists(dir_1)){
    dir.create(dir_1)
  }else{
    print("dir exists")
  }
  
  if (!dir.exists(dir_2)){
    dir.create(dir_2)
  }else{
    print("dir exists")
  }
  
  if (!dir.exists(dir_3)){
    dir.create(dir_3)
  }else{
    print("dir exists")
  }
sample_name <- "dana"
flag_x <- 0
flag_y <- 0

name_x <- "R"
name_y <- "NR"

list_x <- list() # To store seurat objects
list_y <- list()

list_all_x <- list() # TO store sample names
list_all_y <- list()


######################

#### Function to create Seurat object for all the samples
makeObject <- function(dataRead,sampleName){
  
  al1 <- CreateSeuratObject(counts = dataRead$`Gene Expression`,project = sampleName)
  al1[["ADT"]] <- CreateAssayObject(counts = dataRead$`Antibody Capture`)
  al1 <- PercentageFeatureSet(al1, "^MT-", col.name = "percent_mito")
  al1 <- PercentageFeatureSet(al1, "^RP[SL]", col.name = "percent_ribo")
  al1 <- PercentageFeatureSet(al1, "^HB[^(P)]", col.name = "percent_hb")
  al1 <- PercentageFeatureSet(al1, "PECAM1|PF4", col.name = "percent_plat")
  selected_c <- WhichCells(al1, expression = nFeature_RNA > 200 & nCount_RNA < 10000 & percent_mito < 20 & percent_ribo > 5)
  al1 <- subset(al1, cells = selected_c)
  
}


##########



for (x in dirs_group01) {
  list_all_x[[flag_x + 1]] <- paste0(name_x,"_",as.character(flag_x+1))
  namex <- paste0(name_x,as.character(flag_x+1))
  namex <- Read10X(data.dir = paste0("data/Group01/", x))
  namexx <- paste0(name_x,"_",as.character(flag_x+1))
  namexx <- makeObject(namex,paste0(name_x,"_",as.character(flag_x+1)))
  list_x[[flag_x + 1]] <- namexx
  flag_x <- flag_x + 1
  
}
flag_x <- 0



for (y in dirs_group02) {
  list_all_y[[flag_y + 1]] <- paste0(name_y,"_",as.character(flag_y+1))
  namey <- paste0(name_y,as.character(flag_y+1))
  namey <- Read10X(data.dir = paste0("data/Group02/", y))
  nameyy <- paste0(name_y,"_",as.character(flag_y+1))
  nameyy <- makeObject(namey,paste0(name_y,"_",as.character(flag_y+1)))
  list_y[[flag_y + 1]] <- nameyy
  flag_y <- flag_y + 1
  
  
}
flag_y <- 0



print("Read 10X genome files......sucess")
all_samples_names <- unlist(c(list_all_x,list_all_y))
all_sample_objects <- unlist(c(list_x,list_y))
rest_samples <- all_sample_objects[2:length(all_sample_objects)]

group01_objects <- unlist(list_all_x)
group02_objects <- unlist(list_all_y)


merged_seurat <- merge(all_sample_objects[[1]], y = rest_samples,add.cell.ids = all_samples_names,project = 'DanaFarber')
#View(merged_seurat@meta.data)
# create a sample empty column in merged_seurat
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# split sample column into three more columns separated by _
###                                   sample  No  Bar code
###   NR_1_AAACCTGAGATGTGTA-1 becomes   NR      1   AAACCTGAGATGTGTA-1
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('sample', 'No', 'Barcode'), 
                                    sep = '_')
#############################################################################################
# calculate mitochondrial percentage
## We don't need to run below two lines again. We already ran  ran this  step at the begining in line 35
##merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')
##merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 800 & nFeature_RNA > 500 & mitoPercent < 20)
#############################################################################################

merged_seurat_filtered <- merged_seurat

# perform standard workflow steps to figure out if we see any batch effects --------
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)


merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)


tiff(paste(dir_2,'dimplot_00.tiff',sep="/"),units="in", width=7, height=6, res=300)
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'sample', label = TRUE, repel = TRUE,raster=FALSE)
p1
dev.off()

# perform integration to correct for batch effects ------
obj.list <- SplitObject(merged_seurat_filtered, split.by = 'sample')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 1000)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features,dims = 1:30,reduction = "rpca",reference = c(1, length(obj.list)),k.anchor = 5)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors, dims = 1:20)
seurat.integrated_copy <- seurat.integrated
seurat.integrated <- seurat.integrated_copy

DefaultAssay(seurat.integrated) <- "integrated"


# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated,npcs = 50)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:20)
seurat.integrated <- FindNeighbors(seurat.integrated, reduction = "pca", dims = 1:20)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 1)

seurat.integrated_1 <- seurat.integrated

DefaultAssay(seurat.integrated_1) <- "RNA"
##DefaultAssay(seurat.integrated_1) <- "ADT" antibody-derived tags (ADT)

path <- paste(sample_name,".rds",sep = "")

saveRDS(seurat.integrated_1, file = paste(dir_3,path,sep="/"))

print("CITE SEQ Analysis Completed...................................")
return("Finished !!!!!!!!!!!!!!!1")
}

########################################################################################################################


user_function <- function (dirs_group01) {
  
  
  for (x in dirs_group01) {
    y <- x + 1
    print(paste(x,"-hello"))
    
  }
  
  
}




ui <- fluidPage(

  

  titlePanel("CITE SEQ Analysis using Seurat R Library...."),
  textInput("sample_name", "Enter sample name text to create directoy:",value = "dana"),
  actionButton("runScript", "Run: To Create Directory",position = "left",style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
  tableOutput("result"),
  
  br(),
  
  textInput("directory","Enter directory name:",""),
  tableOutput("dir"),
  
  textInput("group1","Enter Group 1 Samples seperated by a comma:",""),
  tableOutput("grp1"),
  
  textInput("group2","Enter Group 2 Samples seperated by a comma:",""),
  tableOutput("grp2"),
  
  br(),
  
  textAreaInput("function_code", "Define a function in R code:", ""),
  
  
  actionButton("run_button", "Run: To Create CITE SEQ Object:",position = "left",style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
  tableOutput("result_text"),
  

  sidebarLayout(position = "right",
                
                sidebarPanel (
                  
                  selectInput("sku", "Check Samples in each group:", skus),
                  textOutput("skuchoice"),
                  
              
                  selectInput('selectfile3','Select Group1 samples from list',choice = list.dirs(path = 'all_data/', recursive = F, full.names = F),multiple = TRUE),
                  textOutput('fileselected3'),

                  
                  selectInput('selectfile4','Select Group2 samples from list',choice = list.dirs(path = 'all_data/', recursive = F, full.names = F),multiple = TRUE),
                  textOutput('fileselected4'),
                  
                  
                            ),
                
                mainPanel("output stuff here", tableOutput('table'))
                
            )
  
          )


server <- function(input, output,session) {
  
  output$skuchoice <- renderText({paste("You choose:",input$sku)})
  
  output$fileselected3 <- renderText({paste("You choose:",input$selectfile3," ")})
  output$fileselected4 <- renderText({paste("You choose:",input$selectfile4,"  ")})
  
  

  
  output$dir <- renderPrint({
    cat("Directory name:\n")
    cat(input$directory)
  })
  
  output$grp1 <- renderPrint({
    x <- as.character(unlist(strsplit(input$group1,",")))
    cat("Samples in Group 1:\n")
    cat(x)
    })
  
  output$grp2 <- renderPrint({
    y <- as.character(unlist(strsplit(input$group2,",")))
    cat("Samples in Group 2:\n")
    cat(y)
    })
  
#########################################################################################################################################
  mylist <- reactiveVal() 
  observe({ mylist(list(sample_name = input$sample_name))})
  

  output$result <- renderTable({
    input$runScript 
    isolate(as.data.frame(testing(input$sample_name)))

  })
  
  values_vector_group1 <- reactive({
    grp1_input <- input$group1
    values_grp1 <- as.character(unlist(strsplit(grp1_input, ",")))
    values_grp1
  })
  
  values_vector_group2 <- reactive({
    grp2_input <- input$group2
    values_grp2 <- as.character(unlist(strsplit(grp2_input, ",")))
    values_grp2
  })
  
  run_function <- eventReactive(input$run_button, {
    function_code <- input$function_code
    values_grp1 <- values_vector_group1()
    values_grp2 <- values_vector_group2()
    
  
    run_citeseq <- eval(parse(text = function_code))
    
    # Call the user-defined function on the vector
    result <- run_citeseq(values_grp1,values_grp2)
    result
  })
  
  
  output$result_text <- renderPrint({
    result <- run_function()
    if (!is.null(result)) {
      result
    } else {
      "Function result is NULL. Please check your function code and input."
    }
  })
  
  output$result_text1 <- renderPrint({
    run_function()
  })
  

}

  
  
  
  
  
  
  
  
  
  
  
  
  
  # mylist2 <- reactiveVal() 
  # observe({ mylist2(list(sample_name = input$sample_name,dirs_group01=input$dirs_group01,dirs_group02=input$dirs_group02))})
  # 
  # output$result2 <- renderTable({
  # 
  #   input$runScript2 
  #   isolate(as.data.frame(run_citeseq(input$sample_name,input$dirs_group01,input$dirs_group02)))
  # })
  # output$table1 <- renderDataTable({sample <- input$run
  #       if (sample == "start") {
  #       source("multiple_samples_cite_seq.R")}
  # })
  


shinyApp(ui, server)
