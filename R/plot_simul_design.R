################################################################
#######
####### Plots in Section 5.1
#######
################################################################

library("tidyverse")

load("../data/simul_design.RData")
exprs <- expand.grid(vtype = c("CTE", "PTA"),
                     atype = c("constant", "uniform"))
exprs <- exprs[1:3, ]

for (k in 1:nrow(exprs)){
    plot <- results %>%
        filter(vtype == exprs$vtype[k],
               atype == exprs$atype[k]) %>%
        mutate(method = factor(method, levels = c("UNW", "IPW", "RIPW"),
                               labels = c("Unweighted", "IPW", "RIPW"))) %>%
        ggplot(aes(x = method, y = est)) +
        geom_boxplot() +
        labs(x = "Method", y = "Estimates") +
        geom_hline(yintercept = 0, color = "red") +
        theme_bw() +
        theme(panel.grid = element_blank(),
              text = element_text(size = 20))
    filename <- paste0("../figs/simul_design_", exprs$vtype[k], "_", exprs$atype[k], "_bias.pdf")
    ggsave(file = filename, plot, width = 6, height = 5)
}

results %>% group_by(vtype, atype, method) %>%
    summarize(coverage = mean(pval >= 0.05))
