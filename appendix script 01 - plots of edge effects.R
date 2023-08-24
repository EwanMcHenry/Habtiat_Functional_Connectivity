
# info for development of functional connectivity metric,

##------ Thu Mar 10 11:07:10 2022 ------##
# plots of expert opinion  etc
timestart = Sys.time()

# working directories ----
maindrive = "S:\\Users\\Ewan McHenry\\OneDrive - the Woodland Trust"
maindrive = "C:\\Users\\emc2\\OneDrive - The Woodland Trust"

ts.wd = paste0(maindrive , "\\Treescapes analysis")
gis.wd = paste0( maindrive, "\\GIS")

# read subscript 01 - load libraries and functions ----
source(paste0(gis.wd, 
              "\\Connectivity\\Functional connectivity\\functional conectivity metric dev\\subscript01- Connectivity metric - loading libraries and functions 01.R"))

# load data
eycott = read.csv(paste0(gis.wd, "\\Connectivity\\Functional connectivity\\functional conectivity metric dev\\hab costs and edge effects Eycott 2011.csv"))


eycott$hab.short = eycott$hab 
eycott$hab.short[c(4:7, 10, 12, )]
  
  factor(eycott$hab, levels = eycott$hab)



ggplot(data = eycott, aes(x = hab, y = eycott.edge.extent)) +
  geom_bar(stat="identity", width=0.5) +
  labs(x = NULL, y = "Edge extent (m)")+
  scale_x_discrete(guide = guide_axis(angle = -60)) +
  theme_pubr()
