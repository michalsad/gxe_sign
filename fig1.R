library(ggplot2)


N <- 1000

set.seed(1)

G <- rbinom(N, size = 2, prob = .4)
E <- rnorm(N)

y <- 10 + 4*G + E + rnorm(N)

pldat <- data.frame(G = as.factor(G), E = E, Y = y)
pldatl <- data.frame(G = as.factor(G), E = E, Y = log(y))

lmod <- lm(y ~ G*E)
coefs <- coef(summary(lmod))

lmodl <- lm(log(y) ~ G*E)
coefsl <- coef(summary(lmodl))

vcols <- c("#1B9E77", "#D95F02", "#4169E1")
           
ggplot(pldat, aes(x=E, y=Y, color=G)) + 
  geom_point(position = position_jitter(w = 0.08, h = 0)) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  theme_bw(base_size = 20) +
  theme(axis.text = element_text(colour = "black"),
        plot.margin = unit(c(.2, .2, .2, .2), "cm"),
        legend.title = element_blank(),
        legend.margin = margin(c(0, 0, -10, 0)),
        legend.position = "top") + 
  scale_color_manual(values = vcols,
                     labels = c("AA", "AB", "BB")) +
  annotate("text", x=Inf, y=-Inf, hjust=1.2, vjust=-1, size = 5,
           label = sprintf("paste(italic(P) [GxE], \" = %.2f\")", coefs[4, 4]), 
           parse = TRUE) +
  xlab(expression(italic("E"))) +
  ylab(expression(italic("Y")))
ggsave("plot_gxe.pdf", width = 6.5, height = 6.5)

ggplot(pldatl, aes(x=E, y=Y, color=G)) + 
  geom_point(position = position_jitter(w = 0.08, h = 0)) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  ylab("log(Y)") +
  theme_bw(base_size = 20) +
  theme(axis.text = element_text(colour = "black"),
        plot.margin = unit(c(.2, .2, .2, .2), "cm"),
        legend.title = element_blank(),
        legend.margin = margin(c(0, 0, -10, 0)),
        legend.position = "top") + 
  scale_color_manual(values = vcols,
                     labels = c("AA", "AB", "BB")) +
  annotate("text", x=Inf, y=-Inf, hjust=1.2, vjust=-1, size=5,
           label = sprintf("paste(italic(P) [GxE], \" = %.2e\")", coefsl[4, 4]), 
           parse = TRUE) +
  xlab(expression(italic("E"))) +
  ylab(expression(paste("log(", italic("Y"), ")")))
ggsave("plot_gxe_log.pdf", width = 6.5, height = 6.5)
