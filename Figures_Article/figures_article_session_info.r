# script session info
pathout <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/"
library("devtools")
library(magrittr)
session_info() %>% capture.output(file=paste0(pathout,"figures_article_session_info.txt"))