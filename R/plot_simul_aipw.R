library("tidyverse")

load("../data/simul_aipw.RData")
plot <- results %>%
    filter(method != "IPW") %>%
    mutate(method = factor(method, levels = c("UNW", "AIPW_fe", "AIPW_nfe", "AIPW_mundlak", "RIPW"),
                           labels = c("Unweighted", "AIPW (w/ FE)", "AIPW (w/o FE)", "AIPW (A.I.)", "RIPW"))) %>%
    mutate(out = factor(out, levels = c(TRUE, FALSE),
                        labels = c("Correct Outcome Model", "Wrong Outcome Model")),
           treat = factor(treat, levels = c(TRUE, FALSE),
                          labels = c("Correct Treatment Model", "Wrong Treatment Model"))) %>%
    ggplot(aes(x = method, y = est)) +
    geom_boxplot() +
    labs(x = "Method", y = "Estimates") +
    geom_hline(yintercept = 0, color = "red") +
    facet_grid(out ~ treat) + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 20),
          axis.text = element_text(size = 10))
filename <- paste0("../figs/simul_aipw.pdf")
ggsave(file = filename, plot, width = 12, height = 6)
