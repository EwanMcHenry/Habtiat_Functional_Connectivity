## save individual hex plots together as one pdf ----
# ggsave(
#   filename = paste0(func.conect.path, 
#                     "\\analysis outputs\\.maps\\",
#   "Treescapes_", paste(years.considered, collapse = "_"),"_hex.ECA.pdf"), 
#   plot = marrangeGrob(eca.hexmap, nrow=1, ncol=1, top = NULL), 
#   width = stand.plot.width , height = stand.plot.height
# )

## plotly for each individual map ----  
# eca.hexmap.plotly <- lapply(seq_along(all.hexgrids), FUN = function(i) {
#   plotly_build(ggplotly(eca.hexmap[[i]], tooltip = "text", 
#                         dynamicTicks = T) %>%
#                  config(displayModeBar = FALSE) %>% layout(hoverlabel = list(align = "left")))
#   # 
#   # all.hexgrids[[i]] = append(all.hexgrids[[i]], print(eca.hexmap) ) %>%
#   #   append(., print(eca.hexmap.plotly) )
# })
