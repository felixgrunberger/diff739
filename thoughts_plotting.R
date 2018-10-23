input <- 0
input$gene <- "PF0740"
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




box_plot <- ggplot(data = plot_everything, aes(x = strain, y = value, color = strain)) +
  geom_boxplot() +
  xlab("strain") + 
  ylab("counts") +
  ggtitle(input$gene) +
  theme_light() +
  scale_color_manual(values = rev(lacroix_palette("PeachPear", n = 4,type = "continuous")))
ggplotly(box_plot)




### now lets create something interesting 
# therefore i need a annotation file ? --> stored in old_new
# get surrounding 10 (?! proteins)
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
         n = as.factor(rep(1:12,11))) 
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

coord_plot <- ggplot(data = coordinate_data_strand, aes(x = coordinate, y = mean_value, color = strain)) +
  geom_segment(aes(x = new_start, y = mean_value, xend = new_end, yend = mean_value),
               lineend = "butt", linejoin = "mitre", 
               size = 1, arrow = arrow(length = unit(0.1, "inches"))) +
  #geom_point(aes(x = new_end, y = mean_value),shape = 18) +
  geom_text(aes(x=new_start-(new_start-new_end)/2, y=-max(mean_value)/10,  
                label = old_lt), angle = 45, color = "black") +
  scale_y_continuous(limits = c(-max(coordinate_data_strand$mean_value)/7, max(coordinate_data_strand$mean_value))) +
  xlab("genomic position") + 
  ylab("counts") +
  ggtitle(paste("genomic region around ", input$gene, sep = "")) +
  theme_light() +
  scale_color_manual(values = rev(lacroix_palette("PeachPear", n = 4,type = "continuous")))


ggiraph(print(coord_plot))

a <- left_join(coordinate_data_strand, full_table_ids, by = c("old_lt"))
a$product



### new heatmaps



