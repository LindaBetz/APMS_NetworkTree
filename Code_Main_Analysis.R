
tiff(filename = "Main_Network_Plots/women.tiff",
    width = 900,
    height = 900)
qgraph(
  cor(data_missings_removed %>% filter(sex == "female") %>% .[, 1:6],
      method = "pearson"),
  graph = "pcor",
  vsize = 10,
  edge.width = 1.75,
  label.cex = 1.75,
  cut = 0,
  borders = T,
  border.width  = 1,
  width = 1,
  minimum = 0,
  maximum = 0.4,
  edge.labels = TRUE,
  edge.label.color = "black",
  edge.label.position = 0.56,
  edge.label.cex = 1.45,
  theme = "colorblind",
  layout = "circle",
  labels = c("1", "2", "3", "4", "5", "6"),
  color =  c(rep("#f7f5f2", 4),
             rep("#c5ceed", 2))
)
dev.off()

tiff(filename = "Main_Network_Plots/men.tiff",
    width = 900,
    height = 900)
qgraph(
  cor(data_missings_removed %>% filter(sex == "male") %>% .[, 1:6],
      method = "pearson"),
  graph = "pcor",
  vsize = 10,
  edge.width = 1.75,
  label.cex = 1.75,
  cut = 0,
  borders = T,
  border.width  = 1,
  width = 1,
  minimum = 0,
  maximum = 0.4,
  edge.labels = TRUE,
  edge.label.color = "black",
  edge.label.position = 0.56,
  edge.label.cex = 1.45,
  theme = "colorblind",
  layout = "circle",
  labels = c("1", "2", "3", "4", "5", "6"),
  color =  c(rep("#f7f5f2", 4),
             rep("#c5ceed", 2))
)
dev.off()
