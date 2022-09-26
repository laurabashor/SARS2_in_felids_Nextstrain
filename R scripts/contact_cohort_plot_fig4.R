###variant plot mimicking outbreak.info plots ########
# see: https://github.com/outbreak-info/R-outbreak-info
# figure 4 in ms

library(tidyverse)
library(ggpubr)
library(ragg)

MUTATIONPALETTE = c('#fff7f3','#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a')
borderColour = "#555555"

dfA <- df2 %>%
  filter(dataset_ID %in% c("Passage_3", "Cat_1","Cat_5", "Cat_6")) %>%
  na.omit() %>%
  mutate(dataset_ID = fct_recode(dataset_ID,
                                 "P3" = "Passage_3",
                                 "Cat 1" = "Cat_1",
                                 "Cat 5" = "Cat_5",
                                 "Cat 6*" = "Cat_6"),
         dataset_ID = fct_relevel(dataset_ID, "P3",
                                  "Cat 1",
                                  "Cat 5",
                                  "Cat 6*"))
dfB <- df2 %>%
  filter(dataset_ID %in% c("Passage_2", "Cat_7", "Cat_8", "Cat_9","Cat_10",  
                           "Cat_11" ,  "Cat_12")) %>%
  na.omit() %>%
  mutate(dataset_ID = fct_recode(dataset_ID,
                                 "P2" = "Passage_2",
                                 "Cat 7" = "Cat_7",
                                 "Cat 8" = "Cat_8",
                                 "Cat 9" = "Cat_9",
                                 "Cat 10*" ="Cat_10",
                                 "Cat 11*" ="Cat_11",
                                 "Cat 12*" ="Cat_12"),
         dataset_ID = fct_relevel(dataset_ID, c("P2", "Cat 7","Cat 8",
                                                "Cat 9","Cat 10*","Cat 11*",
                                                "Cat 12*")))

dfC <- df2 %>%
  filter(dataset_ID %in% c("Passage_3",
                           "Cat_22", "Cat_23",  "Cat_24"  ,
                           "Cat_25" , "Cat_26")) %>%
  na.omit() %>%
  mutate(dataset_ID = fct_recode(dataset_ID,
                                 "P3" = "Passage_3",
                                 "Cat 23" = "Cat_23",  
                                 "Cat 24" = "Cat_24",
                                 "Cat 25" = "Cat_25",
                                 "Cat 22" = "Cat_22",
                                 "Cat 26*" ="Cat_26"),
         dataset_ID = fct_relevel(dataset_ID, "P3"))

cohort_plot <- function(df) {
  ggplot(df, aes(x= fct_reorder(variant, position),
                 y = dataset_ID)) +
    geom_tile(aes(color="#dedede"), color=borderColour) +
    geom_tile(aes(fill = frequency),
              color = borderColour) +
    scale_fill_gradientn(colours = MUTATIONPALETTE, limits = c(0,1),
                         # na.value = "white", 
                         labels = scales::percent) +
    scale_color_manual(values=NA) +              
    guides(color=guide_legend("No data", 
                              override.aes=list(color="white"))) +
    labs(fill="Variant\nallele\nfrequency\n", color = "No data") +
    theme_classic() +
    coord_fixed()+
    theme_minimal() +
    theme(axis.line = element_blank(),
          axis.title= element_blank(),
          axis.ticks=element_blank(),
          axis.text.x = element_text(angle = 45,
                                     vjust = 0.9,
                                     hjust=1,
                                     size=10),
          panel.grid = element_blank(),
          legend.position = "top",
          legend.text  = element_text(size = 10),
          legend.title = element_text(size=12),
          axis.text = element_text(size=12),
          legend.background=element_rect(color="black"))
}

pA <- cohort_plot(dfA)
pA
pB <- cohort_plot(dfB)
pB
pC <- cohort_plot(dfC)
pC

x <- ggplot() + theme_void()

contact_cohorts_plot <- ggarrange(pA, x, pB, x, pC, ncol=1, 
                     legend = "right", common.legend = TRUE, heights = c(1, -0.1, 2,-0.1, 2 ),
                     labels=c("a", "", "b", "", "c")) +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) 

contact_cohorts_plot

ragg::agg_tiff("fig4.tiff", width = 12, height = 9, 
               units = "in", res = 500, compression = "lzw")
contact_cohorts_plot
dev.off()

