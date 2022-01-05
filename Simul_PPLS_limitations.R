### The following script allows to replicate the simulations studies proposed by Etievant and Viallon in Section 3 of "On some 
# limitations of probabilistic models for dimension-reduction: illustration in the case of probabilistic formulations of Partial 
# Least Squares". ###

source('helper_functions.R')

## Specify the following arguments, then run
dir_data_sim1 = paste0("~/DIRECTORY_TO_SPECIFY/") # directory where to save the true parameters values, generated data
                                                  # and estimates for the first simulation study.

dir_data_sim2 = paste0("~/DIRECTORY_TO_SPECIFY/") # directory where to save the true parameters values, generated data
                                                  # and estimates for the second simulation study.

Ncores = NUMBER_TO_SPECIFY                        # number of cores to use to run the data generation and parameters 
                                                  # estimation simultaneously.


## Default values to replicate the simulation studies
Nsim    = 1000
pX      = 20
pM      = 20
r       = 3
RATIO   = 0.25
N       = c(50, 250, 500, 10^3, 5000)


## Steps
run_Simul(dir_data_sim1, dir_data_sim2, Nsim, pX, pM, r, RATIO, N, Ncores)

dir_data = dir_data_sim1
setwd(dir_data)
load(paste0("Dataframe_weights_comparison-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", RATIO, ".RData"))

## Figures
res2$n = factor(res2$n, levels = c("50", "250", "500", "1000", "5000"))
res2$n = factor(res2$n,levels(res2$n)[c(3,2,4,1,5)])
res2$Loading = factor(res2$Loading, labels = c("B", "A"))
res2$Column = factor(res2$Column, labels = c("Column 1", "Column 2", "Column 3"))
res3 = res2[res2$`Comparison with` == "true", ]
res4 = res2[res2$`Comparison with` == "PPLS", ]
res5 = res4[res4$`data generated from` == "misspecified", ]
#ggplot(res5, aes(x = n, y = `Inner Product`, color = Method)) + geom_boxplot(outlier.size = 0.1) + facet_grid(`Loading` ~ `Column`) + theme_light() + theme(legend.position="bottom", plot.title = element_text(size = 15), axis.title.y = element_text(size = 0), axis.title.x = element_text(size = 20), axis.text.y = element_text(size = 18), axis.text.x = element_text(size = 18), legend.title = element_text(size = 20), legend.text = element_text(size = 16), legend.spacing = unit(3, "lines"), strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20)) + labs(color = "Method") + scale_linetype_discrete(labels = c(expression(C[i]), expression(W[i])))  + ggtitle("Comparison with PPLS loadings") 
#ggsave("Comparison_Distribution.pdf", width = 30, height = 20, units = "cm")
#ggplot(res3, aes(x = n, y = `Median Inner Product`, color = Method, linetype = Loading, group = Type_tot, shape = Method)) + geom_point(size = 3) + geom_line(size = 0.5) +  scale_shape_manual(values = c(0,1,3)) + facet_grid(`data generated from` ~ `Column`) + theme_light() + theme(legend.position="bottom", plot.title = element_text(size = 15), axis.title.y = element_text(size = 0), axis.title.x = element_text(size = 20), axis.text.y = element_text(size = 18), axis.text.x = element_text(size = 18), legend.title = element_text(size = 20), legend.text = element_text(size = 16), legend.spacing = unit(3, "lines"), strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20)) + labs(color = "Method") + scale_linetype_discrete(labels = c(expression(C[i]), expression(W[i])))
#ggsave("Comparison_True_Weights.pdf", width = 30, height = 20, units = "cm")
#ggplot(res4, aes(x = n, y = `Median Inner Product`, color = Method, linetype = Loading, group = Type_tot, shape = Method)) + geom_point(size = 3) + geom_line(size = 0.5) +  scale_shape_manual(values = c(0,1,3)) + facet_grid(`data generated from` ~ `Column`) + theme_light() + theme(legend.position="bottom", plot.title = element_text(size = 15), axis.title.y = element_text(size = 0), axis.title.x = element_text(size = 20), axis.text.y = element_text(size = 18), axis.text.x = element_text(size = 18), legend.title = element_text(size = 20), legend.text = element_text(size = 16), legend.spacing = unit(3, "lines"), strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20)) + labs(color = "Method") + scale_linetype_discrete(labels = c(expression(C[i]), expression(W[i])))
#ggsave("Comparison_EMestimates_Weights.pdf", width = 30, height = 20, units = "cm")

# The same graphs, but without the comparison with the weights obtained via PLS-W2A on (X,M), are given below.
res8 = res5[-which(res5$Method=="PLS-W2A"),]
res6 = res3[-which(res3$Method=="PLS-W2A"),]
res7 = res4[-which(res4$Method=="PLS-W2A"),]
ggplot(res8, aes(x = n, y = `Inner Product`, color = Method)) + geom_boxplot(outlier.size = 0.1) + scale_color_manual(values=c("red","blue")) + facet_grid(`Loading` ~ `Column`) + theme_light() + theme(legend.position="bottom", plot.title = element_text(size = 15), axis.title.y = element_text(size = 0), axis.title.x = element_text(size = 20), axis.text.y = element_text(size = 18), axis.text.x = element_text(size = 18), legend.title = element_text(size = 20), legend.text = element_text(size = 16), legend.spacing = unit(3, "lines"), strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20)) + labs(color = "Method") + scale_linetype_discrete(labels = c(expression(C[i]), expression(W[i])))  + ggtitle("Comparison with PPLS loadings") 
ggsave("Comparison_Distribution_woPLSW2A.pdf", width = 30, height = 20, units = "cm")
ggplot(res6, aes(x = n, y = `Median Inner Product`, color = Method, linetype = Loading, group = Type_tot, shape = Method)) + geom_point(size = 3) + geom_line(size = 0.95) +  scale_shape_manual(values = c(0,1,3)) + scale_color_manual(values=c("red","blue", "black")) + facet_grid(`data generated from` ~ `Column`) + theme_light() + theme(legend.position="bottom", plot.title = element_text(size = 15), axis.title.y = element_text(size = 0), axis.title.x = element_text(size = 20), axis.text.y = element_text(size = 18), axis.text.x = element_text(size = 18), legend.title = element_text(size = 20), legend.text = element_text(size = 16), legend.spacing = unit(3, "lines"), strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20)) + labs(color = "Method") + scale_linetype_discrete(labels = c(expression(C[i]), expression(W[i])))
ggsave("Comparison_True_Weights_woPLSW2A.pdf", width = 30, height = 20, units = "cm")
ggplot(res7, aes(x = n, y = `Median Inner Product`, color = Method, linetype = Loading, group = Type_tot, shape = Method)) + geom_point(size = 3) + geom_line(size = 0.95) +  scale_shape_manual(values = c(0,1)) + scale_color_manual(values=c("red","blue")) + facet_grid(`data generated from` ~ `Column`) + theme_light() + theme(legend.position="bottom", plot.title = element_text(size = 15), axis.title.y = element_text(size = 0), axis.title.x = element_text(size = 20), axis.text.y = element_text(size = 18), axis.text.x = element_text(size = 18), legend.title = element_text(size = 20), legend.text = element_text(size = 16), legend.spacing = unit(3, "lines"), strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20)) + labs(color = "Method") + scale_linetype_discrete(labels = c(expression(C[i]), expression(W[i])))
ggsave("Comparison_EMestimates_Weights_woPLSW2A.pdf", width = 30, height = 20, units = "cm")
