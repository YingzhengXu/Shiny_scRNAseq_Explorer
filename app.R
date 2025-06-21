library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(dplyr)
library(Seurat)
library(stringr)
library(ggrepel)
library(ggplot2)
library(viridis)

fgsea_sets<-readRDS('fgsea_sets.RDS')

ui <- dashboardPage(
  dashboardHeader(title = "scRNAseq Explorer"),
  dashboardSidebar(
    sidebarMenu(id = "sidebarid",
                ############
                #CLUSTER#
                menuItem("CLUSTER", tabName = 'clusterpage'),
                conditionalPanel(
                  'input.sidebarid == "clusterpage"',
                  sliderInput("width", label = "Plot width:", min = 300,max = 5000,value = 1000, step = 50),
                  sliderInput("height", label = "Plot height:", min = 300,max = 5000,value = 1000, step = 50),
                  selectizeInput('visclu',"Visulization type: ", c('tsne', 'umap', 'pca'), selected = "tSNE",multiple = FALSE),
                  checkboxInput('quickplot','Speed plot'),
                  
                  selectizeInput("groupclu", "Group by: ", choice=c('RNA_snn_res.0.1','RNA_snn_res.0.2','RNA_snn_res.0.3','RNA_snn_res.0.4','RNA_snn_res.0.5',
                                                                    'RNA_snn_res.0.6','RNA_snn_res.0.7','RNA_snn_res.0.8', 'RNA_snn_res.0.9','RNA_snn_res.1',
                                                                    'RNA_snn_res.1.1','RNA_snn_res.1.2',
                                                                    'singleR','hash.ID', 'HTO_maxID'),
                                 options = list(
                                   maxItems=2,
                                   delimiter = ',',
                                   create = I("function(input, callback){
                                              return {
                                                value: input,
                                                text: input
                                               };
                                            }")),
                                 selected = NULL, multiple = T),
                  selectizeInput("splitclu", "Split by: ", choice=c('RNA_snn_res.0.1','RNA_snn_res.0.2','RNA_snn_res.0.3','RNA_snn_res.0.4','RNA_snn_res.0.5',
                                                                    'RNA_snn_res.0.6','RNA_snn_res.0.7','RNA_snn_res.0.8', 'RNA_snn_res.0.9','RNA_snn_res.1',
                                                                    'RNA_snn_res.1.1','RNA_snn_res.1.2',
                                                                    'singleR','hash.ID', 'HTO_maxID'),
                                 options = list(
                                   maxItems=2,
                                   delimiter = ',',
                                   create = I("function(input, callback){
                                              return {
                                                value: input,
                                                text: input
                                               };
                                            }")),
                                 selected = NULL, multiple = T),
                  downloadButton('downloadcluster', 'Download Plot')
                  
                  
                ),
                ############
                #EXPRESSION#
                menuItem("EXPRESSION", tabName = 'expressionpage'),
                conditionalPanel(
                  'input.sidebarid == "expressionpage"',
                  sliderInput("width2", label = "Plot width:", min = 300,max = 5000,value = 1000, step = 50),
                  sliderInput("height2", label = "Plot height:", min = 300,max = 5000,value = 1000, step = 50),
                  checkboxInput('quickplot2','Speed plot'),
                  
                  radioButtons('visexp',label = NULL, selected = 'tsne', choices = c('tsne','umap','pca')),
                  radioButtons('assaytype',label = NULL, selected = 'RNA', choices = c('RNA')),
                  
                  selectizeInput("gene", "Genes:", choice=NULL, 
                                 options = list(
                                   maxItems=20,
                                   delimiter = ',',
                                   create = I("function(input, callback){
                                              return {
                                                value: input,
                                                text: input
                                               };
                                            }")),
                                 selected = NULL, multiple = T),
                  selectizeInput("splitexp", "Split Plot by: ", choice=c('RNA_snn_res.0.1','RNA_snn_res.0.2','RNA_snn_res.0.3','RNA_snn_res.0.4','RNA_snn_res.0.5',
                                                                         'RNA_snn_res.0.6','RNA_snn_res.0.7','RNA_snn_res.0.8', 'RNA_snn_res.0.9','RNA_snn_res.1',
                                                                         'RNA_snn_res.1.1','RNA_snn_res.1.2',
                                                                         'singleR','hash.ID', 'HTO_maxID'),
                                 options = list(
                                   maxItems=2,
                                   delimiter = ',',
                                   create = I("function(input, callback){
                                              return {
                                                value: input,
                                                text: input
                                               };
                                            }")),
                                 selected = NULL, multiple = T),
                  downloadButton('downloadexpression', 'Download Plot')
                ),
                ############
                #Differential expression table#
                #Pathway tab#
                menuItem("PATHWAY", tabName = 'pathwaypage'),
                conditionalPanel(
                  'input.sidebarid == "pathwaypage"',
                  sliderInput("height.path", label = "Plot height:", min = 300,max = 5000,value = 1000, step = 50),
                  
                  selectizeInput("pathwayinput", "Select pathway:", choice=NULL,
                                 options = list(
                                   maxItems=1,
                                   delimiter = ',',
                                   create = I("function(input, callback){
                                              return {
                                                value: input,
                                                text: input
                                               };
                                            }")),
                                 selected = NULL, multiple = F),
                  tags$style(    type = 'text/css',
                                 ".selectize-input { word-wrap : break-word;}
                         .selectize-dropdown {word-wrap : break-word;} "),
                  selectizeInput("splitpath", "Split by: ", choice=c('RNA_snn_res.0.1','RNA_snn_res.0.2','RNA_snn_res.0.3','RNA_snn_res.0.4',
                                                                     'RNA_snn_res.0.5','RNA_snn_res.0.6','RNA_snn_res.0.7','RNA_snn_res.0.8',
                                                                     'singleR','hash.ID', 'HTO_maxID'),
                                 options = list(
                                   maxItems=2,
                                   delimiter = ',',
                                   create = I("function(input, callback){
                                              return {
                                                value: input,
                                                text: input
                                               };
                                            }")),
                                 selected = NULL, multiple = T),
                  selectizeInput("groupheat", "Group heatmap by: ", choice=c('RNA_snn_res.0.1','RNA_snn_res.0.2','RNA_snn_res.0.3','RNA_snn_res.0.4',
                                                                             'RNA_snn_res.0.5', 'RNA_snn_res.0.6','RNA_snn_res.0.7','RNA_snn_res.0.8', 
                                                                             'singleR','hash.ID', 'HTO_maxID'),
                                 options = list(
                                   maxItems=1,
                                   delimiter = ',',
                                   create = I("function(input, callback){
                                              return {
                                                value: input,
                                                text: input
                                               };
                                            }")),
                                 selected = NULL, multiple = T),
                  actionButton("runHeat", "Do Heatmap"),
                  downloadButton('downloadpath', 'Download Plot')
                ),
                
                ############
                # User DE page#
                menuItem("CUSTOMIZE DE", tabName = 'deuserpage'),
                conditionalPanel(
                  'input.sidebarid == "deuserpage"',
                  selectizeInput("makemeta", "1:Construct metadata ", choice=c('RNA_snn_res.0.1','RNA_snn_res.0.2','RNA_snn_res.0.3','RNA_snn_res.0.4',
                                                                               'RNA_snn_res.0.5','RNA_snn_res.0.6','RNA_snn_res.0.7','RNA_snn_res.0.8',
                                                                               'singleR','hash.ID', 'HTO_maxID'),
                                 options = list(
                                   delimiter = ',',
                                   create = I("function(input, callback){
                                              return {
                                                value: input,
                                                text: input
                                               };
                                            }")),
                                 selected = NULL, multiple = T),
                  selectizeInput("deinput", "2:Select 2 inputs to compare:", choice=NULL,
                                 options = list(
                                   maxItems=2,
                                   delimiter = ',',
                                   create = I("function(input, callback){
                                              return {
                                                value: input,
                                                text: input
                                               };
                                            }")),
                                 selected = NULL, multiple = T),
                  
                  div(style="display:inline-block; float:left",actionButton("runDE", "Do DE")),
                  div(style="display:inline-block; float:left", actionButton("runVol", "Do Volcano")),
                  div(style="display:inline-block; float:left", downloadButton('downloaduserde', 'DE Table')),
                  div(style="display:inline-block; float:right", downloadButton('downloadvoc', 'Volcano'))
                )
    )
  ),
  
  # body ----
  dashboardBody(
    tabItems(
      tabItem(tabName = "clusterpage",
              uiOutput('cluster.ui')
      ),
      tabItem(tabName = "expressionpage",
              uiOutput("expression.ui")
      ),
      tabItem(tabName = "deuserpage",
              column(width = 4,
                     dataTableOutput('deusertable'),style="height:1000px; overflow-y: scroll;overflow-x: scroll;"),
              column(width = 8,
                     uiOutput('volcano.ui')
              )
      ),
      tabItem(tabName = "pathwaypage",
              column(width = 7,
                     uiOutput('pathway.ui')),
              column(width = 5,
                     uiOutput('heatmap.ui')
              )
      )
    )
  )
)

# server -----------------------------------------------------------------------

server <- function(input, output, session){

  seu<-withProgress(message = 'Loading data, be patient...',value = 0.7,{
    seu<-readRDS('sample_seuratobj.rds')
  })
  
  # cluster page
  plotWidth <- reactive(as.numeric(input$width))
  plotHeight <- reactive(as.numeric(input$height))
  plotInput<-reactive({
    if(input$quickplot){
      if(!is.null(input$groupclu) & is.null(input$splitclu)){
        seu$multi<-paste(seu@meta.data[[input$groupclu[1]]], seu@meta.data[[input$groupclu[2]]], sep = '.')
        DimPlot(seu, reduction = as.character(input$visclu), group.by = 'multi', label = T
        )
      }else if(!is.null(input$splitclu) & is.null(input$groupclu)){
        seu$multi<-paste(seu@meta.data[[input$splitclu[1]]], seu@meta.data[[input$splitclu[2]]], sep = '.')
        DimPlot(seu, reduction = as.character(input$visclu), split.by = 'multi', label = T
        )
      }
    }else{
      if(!is.null(input$groupclu) & is.null(input$splitclu)){
        seu$multi<-paste(seu@meta.data[[input$groupclu[1]]], seu@meta.data[[input$groupclu[2]]], sep = '.')
        DimPlot(seu, reduction = as.character(input$visclu), group.by = 'multi', raster = F, label = T, pt.size = 0.5
        )
      }else if(!is.null(input$splitclu) & is.null(input$groupclu)){
        seu$multi<-paste(seu@meta.data[[input$splitclu[1]]], seu@meta.data[[input$splitclu[2]]], sep = '.')
        DimPlot(seu, reduction = as.character(input$visclu), split.by = 'multi', raster = F, label = T, pt.size = 0.5
        )
      }
    }
  })
  
  observe({
    output$cluster <- renderPlot({
      print(plotInput())
    })
    output$cluster.ui <- renderUI({
      plotOutput("cluster", width = plotWidth(), height = plotHeight()) %>% withSpinner(color="#0dc5c1")
    })
    output$downloadcluster <- downloadHandler(
      filename = function() { paste('scatterplot', '.png', sep='') },
      content = function(file) {
        png(file, width = plotWidth(), height = plotHeight())
        print(plotInput())
        dev.off()
      })
  })
  
  
  # expression page
  updateSelectizeInput(session, 'gene', choices = rownames(seu@assays[['RNA']]), server = TRUE, selected =rownames(seu@assays[['RNA']])[1])
  plotInput2<-reactive({
    if(!is.null(input$splitexp)){
      seu$multi<-paste(seu@meta.data[[input$splitexp[1]]], seu@meta.data[[input$splitexp[2]]], sep = '_')
      
      if(input$quickplot2){
        if(input$assaytype == 'RNA'){
          DefaultAssay(seu)<-'RNA'
          FeaturePlot(seu, c(input$gene), reduction = input$visexp,raster=T, split.by = 'multi')
        }
      }else{
        if(input$assaytype == 'RNA'){
          DefaultAssay(seu)<-'RNA'
          FeaturePlot(seu, c(input$gene), reduction = input$visexp,raster=F, split.by = 'multi')
        }
      }
    }else{
      if(input$quickplot2){
        if(input$assaytype == 'RNA'){
          DefaultAssay(seu)<-'RNA'
          FeaturePlot(seu, c(input$gene), reduction = input$visexp,raster=T)
        }
      }else{
        if(input$assaytype == 'RNA'){
          DefaultAssay(seu)<-'RNA'
          FeaturePlot(seu, c(input$gene), reduction = input$visexp,raster=F)
        }
      }
    }
  })
  
  plotWidth2 <- reactive(as.numeric(input$width2))
  plotHeight2 <- reactive(as.numeric(input$height2))
  observe({
    output$expression<-renderPlot({
      plotInput2()
    })
    output$expression.ui <- renderUI({
      plotOutput("expression", width = plotWidth2(), height = plotHeight2())%>% withSpinner(color="#0dc5c1")
    })
    output$downloadexpression <- downloadHandler(
      filename = function() { paste('expression', '.png', sep='') },
      content = function(file) {
        png(file, width = plotWidth2(), height = plotHeight2())
        print(plotInput2())
        dev.off()
      })
  })
  
  # DE user page
  observe({
    if(!is.null(input$makemeta)){
      list.meta<-list()
      for(i in 1:length(input$makemeta)){
        list.meta[[i]]<-seu@meta.data[[input$makemeta[i]]]
      }
      seu$DEmeta<-do.call(paste, c(list.meta, sep='_'))
      updateSelectizeInput(session, 'deinput', choices = names(table(seu$DEmeta)), server = TRUE, selected =NULL)
      
      #do FindMarker
      markers <- eventReactive(input$runDE, {
        Idents(seu)<-'DEmeta'
        temp<-withProgress(message = paste('Calculating:', input$deinput[1], ' vs ', input$deinput[2], sep = " "), value = 0.5,{
          FindMarkers(seu, ident.1 = input$deinput[1], ident.2 = input$deinput[2])
        })
        temp$`Genes`<-rownames(temp)
        return(temp)
      })
      output$deusertable <- renderDataTable({
        markers()
      })
      output$downloaduserde <- downloadHandler(
        filename = function() { paste('Customized_DE', '.csv', sep='') },
        content = function(file) {
          write.csv(markers(), file)
        })
      # make Volcano plot
      volconoPlot<-eventReactive(input$runVol, {
        de<-markers()
        de$diffexpressed <- "NO"
        de$diffexpressed[de$avg_log2FC > 1 & de$p_val < 0.05] <- "UP"
        de$diffexpressed[de$avg_log2FC < -1 & de$p_val < 0.05] <- "DOWN"
        de$delabel <- NA
        de$delabel[de$diffexpressed != "NO"] <- de$Genes[de$diffexpressed != "NO"]
        plot<-withProgress(message = "Rendering Volcano Plot...", value = 0.5,{
          ggplot(data=de, aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed, label=delabel)) +
            geom_point() +
            theme_minimal() +
            geom_text_repel() +
            scale_color_manual(values=c("blue", "black", "red")) + theme_bw()
        })
        return(plot)
      })
      
      observe({
        output$volcano <- renderPlot({
          print(volconoPlot())
        })
        output$volcano.ui <- renderUI({
          plotOutput("volcano", height = 1000) %>% withSpinner(color="#0dc5c1")
        })
        output$downloadvoc <- downloadHandler(
          filename = function() { paste('volcanoplot', '.png', sep='') },
          content = function(file) {
            png(file)
            print(volconoPlot())
            dev.off()
          })
        
      })
      
    }
  })
  
  # Pathway page
  updateSelectizeInput(session, 'pathwayinput', choices = names(fgsea_sets), server = TRUE, selected =names(fgsea_sets)[1])
  #Pathway scatter plot
  # get pathway genes
  plotHeight.path <- reactive(as.numeric(input$height.path))
  
  pathway.df<-as.data.frame(seu@reductions[["tsne"]]@cell.embeddings)
  pathway.df$RNA_snn_res.0.1<-seu$RNA_snn_res.0.1
  pathway.df$RNA_snn_res.0.2<-seu$RNA_snn_res.0.2
  pathway.df$RNA_snn_res.0.3<-seu$RNA_snn_res.0.3
  pathway.df$RNA_snn_res.0.4<-seu$RNA_snn_res.0.4
  pathway.df$RNA_snn_res.0.5<-seu$RNA_snn_res.0.5
  pathway.df$RNA_snn_res.0.6<-seu$RNA_snn_res.0.6
  pathway.df$RNA_snn_res.0.7<-seu$RNA_snn_res.0.7
  pathway.df$RNA_snn_res.0.8<-seu$RNA_snn_res.0.8
  pathway.df$hash.ID<-seu$hash.ID
  pathway.df$HTO_maxID<-seu$HTO_maxID
  pathway.df$singleR<-seu$singleR

  temp.genes<-reactive({
    path.genes<-fgsea_sets[[input$pathwayinput]]
    path.genes<-intersect(path.genes, rownames(seu))
    temp.genes<-as.data.frame(seu@assays[["RNA"]]@data[path.genes,])
    return(temp.genes)
  })
  PlotInput.path<-reactive({
    if(!is.null(input$pathwayinput)){
      sums<-colSums(temp.genes())
      if(!is.null(input$splitpath)){
        if(length(input$splitpath)<2){
          plot<-ggplot(pathway.df, aes(x=tSNE_1, y=tSNE_2, color=sums))+
            geom_point(size=0.5)+
            facet_wrap(~get(input$splitpath[1]))+
            labs(color=NULL)+
            scale_colour_viridis()+ theme_bw()
          return(plot)
        }
        if(length(input$splitpath)>1){
          plot<-ggplot(pathway.df, aes(x=tSNE_1, y=tSNE_2, color=sums))+
            geom_point(size=0.5)+
            facet_grid(get(input$splitpath[1])~get(input$splitpath[2]))+
            labs(color=NULL)+
            scale_colour_viridis()+ theme_bw()
          return(plot)
        }
        
      }else{
        plot<-ggplot(pathway.df, aes(x=tSNE_1, y=tSNE_2, color=sums))+
          geom_point(size=0.5)+
          labs(color=NULL)+
          scale_colour_viridis()+ theme_bw()
        return(plot)
      }
    }
  })
  
  
  #Heatmap
  PlotInput.heat<-eventReactive(input$runHeat, {
    if(!is.null(input$pathwayinput)){
      DefaultAssay(seu)<-'RNA'
      DoHeatmap(subset(seu, downsample = 5000), features = fgsea_sets[[input$pathwayinput]], group.by = input$groupheat)+
        scale_fill_gradientn(colors = c('navy','white','red'))+guides(color=FALSE)
    }
  })
  
  observe({
    output$pathway <- renderPlot({
      print(PlotInput.path())
    })
    output$heatmap <- renderPlot({
      PlotInput.heat()
    })
    output$pathway.ui <- renderUI({
      plotOutput("pathway", height = plotHeight.path()) %>% withSpinner(color="#0dc5c1")
    })
    output$heatmap.ui <- renderUI({
      plotOutput("heatmap", height = plotHeight.path()) %>% withSpinner(color="#0dc5c1")
    })
    
    output$downloadpath <- downloadHandler(
      filename = function() { paste('pathwayplot', '.pdf', sep='') },
      content = function(file) {
        pdf(file)
        print(PlotInput.path())
        print(PlotInput.heat())
        dev.off()
      })
  })
}

shinyApp(ui, server)



