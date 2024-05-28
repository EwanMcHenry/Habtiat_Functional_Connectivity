# Ewan McHenry


# exploration of species dispersal data
# libraries
library(tidyverse)
library(plotly)
library(ggpubr)#publr


# data
dsp_df = read.csv("Data\\species_dispersal_data_woodland.csv")

str(dsp_df)
#lengh non na in each col
sapply(dsp_df, function(x) sum(!is.na(x)))

# plot movement data WITH SPECIES LABELLED
# hover text using ggplotly()
# point distribution jittered verticle



library(tidyverse)
library(hrbrthemes)
library(viridis)

# create a dataset
data <- data.frame(
  name=c( rep("A",500), rep("B",500), rep("B",500), rep("C",20), rep('D', 100)  ),
  value=c( rnorm(500, 10, 5), rnorm(500, 13, 1), rnorm(500, 18, 1), rnorm(20, 25, 4), rnorm(100, 12, 1) )
)

# Interactive plot
plot <- dsp_df %>% 
  filter(woodland.tree.associated == 1) %>%
  group_by(species) %>%
  summarise(movement = median(movement, na.rm = TRUE),
            mean_movement =  median(mean_movement, na.rm = TRUE),
            class = class[1],
            common_name = common_name[1]) %>%
  ggplot( aes(y=class, x=mean_movement, fill=class)) +
  geom_boxplot(outlier.colour = NULL) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(size=0.4, alpha=0.9, 
              aes(text = map(
                paste0("<b>", common_name, "</b>",
                       "<br>", mean_movement, " m"), 
                             HTML),
                  colour = class)) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Mean species movement") +
  xlab("Distance (m)") +
  ylab("Species class")


  ggplotly(plot, tooltip = "text", dynamicTicks = T) %>% 
  config(displayModeBar = F)
  




#dotplot mean_movement on x, y is a jitter centered on 0
dsp_df %>% 
  filter(woodland.specialist == 1) %>%
  group_by(species) %>%
  summarise(movement = median(movement, na.rm = TRUE),
            mean_movement =  median(mean_movement, na.rm = TRUE)) %>%
  ggplot()+
  geom_point(aes(x = mean_movement,
                 y = jitter(rep(0, length(mean_movement)), amount = 1),
                 col = class,
                 text = species),
             show.legend = T) + #no legend
#make ggplotly
  ggplotly(plot, tooltip = "text", dynamicTicks = F) %>% 
  config(displayModeBar = F)
  # log transform the axis
  # scale_x_log10()+
  # add hover text
  ggplotly()

dsp_df %>% 
  #subset woodland species
  filter(woodland.specialist == 1) %>%
  # average movement by species
  group_by(species) %>%
  summarise(movement = median(movement, na.rm = TRUE),
            mean_movement =  median(mean_movement, na.rm = TRUE)) %>%
              ) %>%
  
  ggplot(aes(x = mean_movement))+
  geom_point(aes(col = species,
                 text = species),
             show.legend = FALSE)+ #no legend
  geom_text(aes( label = common_name), nudge_x = 0.1, nudge_y = 0.1) +
  # log transform the axis
  scale_x_log10()+




dsp_df %>% 
  ggplot(aes(x = movement , y = movement ))+
  geom_point(aes(col = species))+
  geom_text(aes( label = species), nudge_x = 0.1, nudge_y = 0.1)
