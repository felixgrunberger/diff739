library("shiny")
library("semantic.dashboard")
library("plotly")
library("ggplot2")
library("tidyverse")
library("heatmaply")
library("LaCroixColoR")
library("DESeq2")
library("XVector")
library("data.table")
library("ggiraph")


#setwd("/Users/felix/Documents/R/differential_0739/shinyapps/diff739")
old_new <- fread("genes_gff_Table_old_new.tsv")
load("copper_resIHW_table")
load("deletion_resIHW_table")
load("full_table_ids_table")
load("dds_table")
#load("mat_copper_table")
#load("mat_deletion_table")
##load("mat_all_table")

Sys.setlocale(locale="en_US.UTF-8")

ui <- dashboardPage(title = "exploratory analysis of RNA seq data",theme = "simplex",
    
  dashboardHeader(
    
  ),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Compare counts", tabName = "single_genes")#,
      #menuItem("Heatmaps", tabName = "heatmap_boxes"#,
      #         menuSubItem("Copper shock",tabName = "copper_boxes"),
      #         menuSubItem("Deltion 739", tabName = "deletion_boxes"),
      ##         menuSubItem("all data", tabName = "all_boxes"))
    )
  ),
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
      ),
    tabItems(
      tabItem(tabName = "single_genes",
              fluidRow(
                box(title = "Select gene of interest", status = "primary", 
                    collapsible = F, width = 5,
                           selectizeInput(inputId = "gene",selected = "PF0739",
                                          label = "",
                                          choices = full_table_ids$old_lt, 
                                          options = list(maxOptions = 3000)
                    )
                    ),
                box(title = "Name in new genome", status = "primary", collapsible = F, width = 5,
                    column(3, 
                           label = "",
                           textOutput("new_tag")
                           )
                    ),
                box(title = "function", status = "primary", collapsible = F, width = 5,
                    
                    column(3,
                           label = "",
                           textOutput("name")
                    )
                    )
                ),
              fluidRow(
    
                box(title = "Show statistics", status = "primary", collapsible = T, width = 15,
                           tableOutput('table')
                    )
                ),
              fluidRow(
                box(title = "Show point-plot", status = "primary", collapsible = T, width = 15,
                           plotlyOutput("plot_any")
                )
                ),
              fluidRow(
                box(title = "Show box-plot", status = "primary", collapsible = T, width = 15,
                    plotlyOutput("box_plot")
                )
              
              ),
              fluidRow(
                box(title = "Show genomic region", status = "primary", collapsible = T, width = 15,
                    ggiraphOutput("region_plot")
                )
                
              )
      )
  )
)
)

server <- shinyServer(function(input, output, session) {
  
  output$new_tag = renderText({
    print(full_table_ids$gene[full_table_ids$old_lt == input$gene])
  })
  
  output$name = renderText({
    print(full_table_ids$product[full_table_ids$old_lt == input$gene])
  })
  
  # table
  output$table <- renderTable({
    copper_test <- copper_resIHW %>%
      as.data.frame() %>%
      tibble::rownames_to_column("genes") %>%
      left_join(old_new,c("genes" = "id")) %>%
      dplyr::filter(old_lt == input$gene) %>%
      as.tibble() %>%
      mutate('adjusted p value' = padj,
             condition = "copper shock vs 52",
             'type of regulation' = ifelse(log2FoldChange < 0, "down", "up")) %>%
      dplyr::select(log2FoldChange,'adjusted p value', 'type of regulation',condition)
    
    deletion_test <- deletion_resIHW %>%
      as.data.frame() %>%
      tibble::rownames_to_column("genes") %>%
      left_join(old_new,c("genes" = "id")) %>%
      dplyr::filter(old_lt == input$gene) %>%
      as.tibble() %>%
      mutate('adjusted p value' = padj,
             condition = "deletion vs 70",
             'type of regulation' = ifelse(log2FoldChange < 0, "down", "up")) %>%
      dplyr::select(log2FoldChange,'adjusted p value', 'type of regulation',condition)
    
    test <- rbind(copper_test, deletion_test)
    test
    
  }
  )
  
  
  #plot
  
  output$plot_any = renderPlotly({
    
    shiny::validate(
      need(input$gene, "What gene are you interested in?")
    )
    
    plot_everything <- counts(dds) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("genes") %>%
      left_join(old_new,c("genes" = "id")) %>%
      dplyr::filter(old_lt == input$gene) %>%
      t() %>%
      as.tibble() %>%
      dplyr::slice(2:13) %>%
      mutate(value = V1,
             dataset = substring(colnames(counts(dds)),76,82),
             strain = substring(dataset,1,2),
             replicate = substring(dataset, 4,4),
             strain = ifelse(substring(dataset,6,7) != "cu", strain, paste(strain,"_cu", sep = "")),
             replicate = substring(replicate, 1, 1)) %>%
      mutate(value = as.numeric(value),
             strain = as.factor(strain),
             replicate = as.factor(replicate))
    
    counts_plot <- ggplot(data = plot_everything,aes(x = strain,y = value, color = replicate)) +
      geom_jitter(width = 0.1, size = 5) + 
      xlab("strain") + 
      ylab("counts") +
      ggtitle(input$gene) +
      theme_light() +
      scale_color_manual(values = rev(lacroix_palette("PeachPear", n = 3,type = "continuous")))
    
    ggplotly(counts_plot)
    
  })
  
  ## BOX PLOT
  output$box_plot = renderPlotly({
    
    shiny::validate(
      need(input$gene, "What gene are you interested in?")
    )
    
    plot_everything <- counts(dds) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("genes") %>%
      left_join(old_new,c("genes" = "id")) %>%
      dplyr::filter(old_lt == input$gene) %>%
      t() %>%
      as.tibble() %>%
      dplyr::slice(2:13) %>%
      mutate(value = V1,
             dataset = substring(colnames(counts(dds)),76,82),
             strain = substring(dataset,1,2),
             replicate = substring(dataset, 4,4),
             strain = ifelse(substring(dataset,6,7) != "cu", strain, paste(strain,"_cu", sep = "")),
             replicate = substring(replicate, 1, 1)) %>%
      mutate(value = as.numeric(value),
             strain = as.factor(strain),
             replicate = as.factor(replicate))
    
    box_plot_data <- ggplot(data = plot_everything, aes(x = strain, y = value, color = strain)) +
      geom_boxplot() +
      xlab("strain") + 
      ylab("counts") +
      ggtitle(input$gene) +
      theme_light() +
      scale_color_manual(values = rev(lacroix_palette("PeachPear", n = 4,type = "continuous")))
    ggplotly(box_plot_data)
    
  })
  
  output$region_plot <- renderggiraph({
    shiny::validate(
      need(input$gene, "What gene are you interested in?")
    )
    which_row <- counts(dds) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("genes") %>%
      left_join(old_new,c("genes" = "id")) %>%
      summarise(what = as.numeric(which(old_lt == input$gene))) %>%
      as.data.frame()
    
    plot_everything1 <- counts(dds) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("genes") %>%
      left_join(old_new,c("genes" = "id")) %>%
      dplyr::slice((which_row$what-5):(which_row$what+5)) %>%
      dplyr::slice(rep(1:n(), each = 12)) %>% 
      mutate(strain = rep(unlist(list(rep(factor(52),3),rep(factor("52_cu"),3),rep(factor(70),3),rep(factor(74),3))),11),
             replicate = as.factor(rep(c(1,2,3),44)),
             n = as.factor(rep(1:12,11)),
             length = (end - start)/1000) 
    
    colnames(plot_everything1)[2:13]  <- c("A","B","C","D","E", "F", "G", "H", "I", "J", "K", "L")
    plot_everything2 <- plot_everything1 %>%
      mutate(value = ifelse(n == 1, A,
                            ifelse(n == 2, B, 
                                   ifelse(n == 3, C,
                                          ifelse(n == 4, D, 
                                                 ifelse(n ==5, E,
                                                        ifelse(n == 6, F, 
                                                               ifelse(n == 7, G,
                                                                      ifelse(n == 8, H, 
                                                                             ifelse(n == 9, I,
                                                                                    ifelse(n == 10 , J,
                                                                                           ifelse(n ==11, K, 
                                                                                                  L)))))))))))) %>%
      rowwise() %>%
      mutate(value = value/length) %>%
      select(old_lt, start, end, strand, strain, replicate,value, genes)
    
    plot_everything2$mean_value <- rep(c(aggregate(plot_everything2$value, by=list((seq(nrow(plot_everything2))-1) %/% 3), FUN=mean)$x),each = 3)
    
    
    coordinate_data <- plot_everything2 %>%
      dplyr::slice(rep(1:n(), each = 2)) 
    coordinate_data$coordinate <- coordinate_data$end
    coordinate_data[seq(1, nrow(coordinate_data), 2), ]$coordinate <- coordinate_data[seq(1, nrow(coordinate_data), 2), ]$start
    
    ## reflecting strand position
    coordinate_data_strand <- coordinate_data %>%
      mutate(new_start = ifelse(strand == "+", start,end),
             new_end = ifelse(strand == "+", end, start))
    
    coordinate_annotation <- left_join(coordinate_data_strand,full_table_ids, by = ("old_lt")) %>%
      mutate(tooltip_label = sprintf("gene: %s \n strain: %s \n counts: %s \n product: %s", old_lt, strain, mean_value, product))
    
    coord_plot <- ggplot(data = coordinate_annotation, aes(x = coordinate, y = mean_value, color = strain)) +
      geom_segment_interactive(aes(tooltip = tooltip_label, x = new_start, y = mean_value, xend = new_end, yend = mean_value),
                   lineend = "butt", linejoin = "mitre", 
                   size = 1, arrow = arrow(length = unit(0.1, "inches"))) +
      geom_text(aes(x=new_start-(new_start-new_end)/2, y=-max(mean_value)/10,  
                    label = old_lt), angle = 45, color = "black") +
      scale_y_continuous(limits = c(-max(coordinate_data_strand$mean_value)/7, max(coordinate_data_strand$mean_value))) +
      xlab("genomic position") + 
      ylab("RPK (reads per kilobase)") +
      ggtitle(paste("genomic region around ", input$gene, sep = "")) +
      theme_light() +
      scale_color_manual(values = rev(lacroix_palette("PeachPear", n = 4,type = "continuous")))
    
    ggiraph(print(coord_plot),width_svg = 10,zoom_max = 5)
  })
  
  ## heatmap copper
  #output$copper_heatmap <- d3heatmap::renderD3heatmap({
  #  d3heatmap(mat_copper, 
  #            colors = rev(lacroix_palette("PeachPear", n = 2000,type = "continuous")))
  ##})
  
  ## heatmap deletion
  #output$deletion_heatmap <- renderD3heatmap({
  #  d3heatmap(mat_deletion, 
  #            colors = rev(lacroix_palette("PeachPear", n = 2000,type = "continuous")))
  #})
  #
  ## heatmap all 4
  #output$all_heatmap <- renderD3heatmap({
  #  d3heatmap(mat, 
  #            colors = rev(lacroix_palette("PeachPear", n = 2000,type = "continuous")))
  #})
  
})

shinyApp(ui, server)
