# expert opinion 

library(tidyverse)
library(ggpubr)
library(cowplot)
library(scales)
library(DT)
library(kableExtra)
library(rmarkdown)



draw_label_theme <- function(label, theme = NULL, element = "text", ...) {
  if (is.null(theme)) {
    theme <- ggplot2::theme_get()
  }
  if (!element %in% names(theme)) {
    stop("Element must be a valid ggplot theme element name")
  }}
  
# working directories ----
maindrive = "D:\\Users\\Ewan McHenry\\OneDrive - the Woodland Trust"
#maindrive = "C:\\Users\\emc2\\OneDrive - The Woodland Trust"
gis.wd = paste0( maindrive, "\\GIS")
func.conect.path = paste0(gis.wd, "\\Connectivity\\Functional connectivity\\functional conectivity metric dev")

setwd(func.conect.path)
# data ----
eo = read.csv("WT expert opinion01_03.22.csv")

# curation ----

eo$edge.low = eo$edge_value - ((1-eo$edge_certainty)/2)
eo$edge.hi = eo$edge_value + ((1-eo$edge_certainty)/2)
eo$aw.low = eo$nonancient_value - ((1-eo$nonancient_certainty)/2)
eo$aw.hi = eo$nonancient_value + ((1-eo$nonancient_certainty)/2)


eo$edge.hi [which(eo$edge.low<0)] = eo$edge.hi[which(eo$edge.low<0)] - eo$edge.low[which(eo$edge.low<0)]
eo$edge.low[which(eo$edge.low<0)] = 0

eo$edge.low[which(eo$edge.hi>1)] = eo$edge.low[which(eo$edge.hi>1)] + 1-eo$edge.hi[which(eo$edge.hi>1)]
eo$edge.hi [which(eo$edge.hi>1)] = 1

eo$aw.hi [which(eo$aw.low<0)] = eo$aw.hi[which(eo$aw.low<0)] - eo$aw.low[which(eo$aw.low<0)]
eo$aw.low[which(eo$aw.low<0)] = 0

eo$id = 1:dim(eo)[1]


# plot results for each round for each expert (their guess in red) ----
this.round = 1
plots = list(NA)
for (this.round in unique(eo$round)){
  for (you in  unique(eo$expert)){
    gg.edge = eo[eo$round == this.round,] %>% 
      ggplot()+
      geom_point(aes(x = edge_value,y = id ))+
      geom_errorbarh(aes(xmin = edge.low, xmax = edge.hi,y = id ))+
      geom_point(data = eo[eo$round == this.round & eo$expert == you,], aes(x = edge_value,y = id ), col = "red") +
      scale_x_continuous(lim = c(0,1), breaks = pretty_breaks(3))+
      ggtitle(paste("Detrimental edge",eo$expert[eo$expert == you & eo$round == this.round] , "round", this.round))+
      labs(x = "value", y = "Expert ID")+
      theme_pubr()
    
    gg.aw = eo[eo$round == this.round,] %>% 
      ggplot()+
      geom_point(aes(x = nonancient_value,y = id ))+
      geom_errorbarh(aes(xmin = aw.low, xmax = aw.hi,y = id ))+
      geom_point(data = eo[eo$round == this.round & eo$expert == you,], aes(x = nonancient_value,y = id ), col = "red") +
      scale_x_continuous(lim = c(0,1), breaks = pretty_breaks(3))+
      ggtitle("Non-ancient woodland")+
      labs(x = "value", y = "Expert ID")+
      theme_pubr()
    
    ggsave( plot_grid(gg.edge, gg.aw), width = 220, height = 120, units = "mm",
            file = paste0("expert opinion\\",eo$round[eo$expert == you & eo$round == this.round] , eo$expert[eo$expert == you & eo$round == this.round],".pdf"))
    # plots =  plot_grid(gg.edge, gg.aw)
    # cat(paste0("Hi, expert opinion results are in! Could you please take 2 min to review and possibly change the value and certainty of your guesess (tell me if you want to stick with what you have).\n",
    #            "Dots show the guesses (yours is in red) and wider whiskers mean less certainty\n\n",
    #            eo$expert[eo$expert == you], ", you said you were: \n", 100*eo$edge_certainty[eo$expert == you], "% certain woodland affected by detrimental edge impacts are ", 
    #            100*eo$edge_value[eo$expert == you],"% the value of non-edge\nAnd... \n", 100*eo$nonancient_certainty[eo$expert == you], "% certain that non-ancient woodland has ", 
    #            100*eo$nonancient_value[eo$expert == you],"% the value of ancient\n\n"))
    
  }
  
}
# plt all rounds ----
plot_grid(eo[eo$round == 1,] %>% 
            ggplot()+
            geom_point(aes(x = edge_value,y = id ))+
            geom_errorbarh(aes(xmin = edge.low, xmax = edge.hi,y = id ))+
            scale_x_continuous(lim = c(0,1), breaks = pretty_breaks(3))+
            ggtitle(paste("Detrimental edge round 1"))+
            labs(x = "value", y = "Expert ID")+
            theme_pubr(),
          eo[eo$round == 1,] %>% 
            ggplot()+
            geom_point(aes(x = nonancient_value,y = id ))+
            geom_errorbarh(aes(xmin = aw.low, xmax = aw.hi,y = id ))+
            scale_x_continuous(lim = c(0,1), breaks = pretty_breaks(3))+
            ggtitle("Non-ancient woodland value round 1")+
            labs(x = "value", y = "Expert ID")+
            theme_pubr(),
          eo[eo$round == 2,] %>% 
            ggplot()+
            geom_point(aes(x = edge_value,y = id ))+
            geom_errorbarh(aes(xmin = edge.low, xmax = edge.hi,y = id ))+
            scale_x_continuous(lim = c(0,1), breaks = pretty_breaks(3))+
            ggtitle(paste("Detrimental edge round 2"))+
            labs(x = "value", y = "Expert ID")+
            theme_pubr(),
          eo[eo$round == 2,] %>% 
            ggplot()+
            geom_point(aes(x = nonancient_value,y = id ))+
            geom_errorbarh(aes(xmin = aw.low, xmax = aw.hi,y = id ))+
            scale_x_continuous(lim = c(0,1), breaks = pretty_breaks(3))+
            ggtitle("Non-ancient woodland value round 2")+
            labs(x = "value", y = "Expert ID")+
            theme_pubr()
          )

# summariese results ----
delphi.summary = eo[,c(1,3:6)] %>% 
  group_by(round) %>% 
  summarise(nonAW.val = median(nonancient_value, na.rm = T),
            AW.cert = median(nonancient_certainty, na.rm = T),
            AW.sd = sd(nonancient_value, na.rm = T),
            edge.val = median(edge_value, na.rm = T),
            edge.cert = median(edge_certainty, na.rm = T),
            edge.sd = sd(edge_value, na.rm = T) )

delphi.point.est.df = data.frame(' ' = c("Favourable core", "Unfavourable edge"),
                           'Non-ancient woodland' = c(1,
                                                      1 * delphi.summary[dim(delphi.summary)[1], "edge.val" ]) %>% 
                             as.numeric(),
                           'Ancient woodland' = c(1/delphi.summary[dim(delphi.summary)[1], "nonAW.val" ],
                                                  1/delphi.summary[dim(delphi.summary)[1], "nonAW.val" ] * delphi.summary[dim(delphi.summary)[1], "edge.val" ]) %>% 
                             as.numeric()
)

# fancy table ----
colnames(delphi.point.est.df) = c(" ", "Non-ancient woodland", "Ancient woodland" )
delphi.point.est.df  %>%
  mutate_if(is.numeric, round, digits = 2) %>% 
  kable(.,"html",  align = c("l", "c", "c"),
        caption = "Relative habiatat quality model parameters") %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)  %>% 
  kable_classic(full_width = F, html_font = "Cambria") %>% 
  footnote( general = "This parametarisation is  an imperfect attempt to account for varying habitat quality in functional connectivity calculations. It does not relate to the value of ancient woodland (which is considered irreplaceable).") 

# save objects ----
write.csv(delphi.point.est.df, 
     file = paste0(func.conect.path, "\\delphi point estimate summary.csv"))  



  
