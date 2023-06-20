library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggplot2)

dat <- read.csv("forBOX.csv", header=T, row.names = 1)
data <- as.data.frame(dat)
data %>%
    ggplot( aes(x=name, y=value, fill=name)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme_ipsum() +
    theme(
        legend.position="none",
        plot.title = element_text(size=11)
    ) +
    ggtitle("A boxplot with jitter") +
    xlab("")

ggplot(dat, aes(x = as.factor(gene), y = TPM)) +
    geom_boxplot(aes(fill = group), position = position_dodge(0.9)) +
    facet_wrap(~ gene, scales = "free") +
    scale_fill_manual(values = c("#09E359", "#E31009", "#E309DF", "#E3D809")) + 
    theme_bw()
