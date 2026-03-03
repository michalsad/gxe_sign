library(dplyr)
library(ggplot2)

h2 <- .5
s2 <- 1
s2.gxe <- .5
n0 <- 5000
n1 <- 5000
n <- n0 + n1
p <- 200
frac.causal.gxe <- .5
frac.positive <- .5
alpha <- 2
pthr <- .05

phi <- function(x) log(x - min(x) + 1e-6)

b <- rnorm(p, sd = sqrt(s2/p))
    
p.gxe <- round(p*frac.causal.gxe)
b.gxe <- rnorm(p.gxe, sd = sqrt(s2.gxe/p.gxe))
b.tmp <- rep(0, p)
b.tmp[1:p.gxe] <- b.gxe
b.gxe <- b.tmp

type <- sign(b)*sign(b.gxe)
type.tmp <- rep(1, length(type))
type.tmp[type == 1] <- 2
type.tmp[type == -1] <- 3
type <- type.tmp

freq <- runif(p, min = .1, max = .5)
X <- sapply(freq, function(x) rbinom(n, size = 2, prob = x))
X <- scale(X)
E <- rep(c(0, 1), c(n0, n1))
E <- c(scale(E))
    
sg <- var(X %*% b)
se <- sg*((1 - h2)/h2)
    
y <- alpha*E + X %*% b + (X*E) %*% b.gxe + rnorm(n, sd = sqrt(se))
y.tilde <- phi(y)

coefs <- lapply(1:ncol(X), function(ind) coef(summary(lm(y ~ X[, ind]*E))))
coefs <- do.call(
  "rbind", 
  lapply(coefs, function(x) c(x[2, c(3, 4)], x[4, c(3, 4)])))

coefs.tilde <- lapply(
  1:ncol(X), 
  function(ind) coef(summary(lm(y.tilde ~ X[, ind]*E))))
coefs.tilde <- do.call(
  "rbind", 
  lapply(coefs.tilde, function(x) c(x[2, c(3, 4)], x[4, c(3, 4)])))

cols <- rep("gray", p)
cols[type == 1] <- "black"
cols[type == 2] <- "royalblue"
cols[type == 3] <- "firebrick2"

coefs <- as.data.frame(coefs)
colnames(coefs) <- c("estG", "pG", "estGxE", "pGxE")
coefs$color <- factor(
  cols,
  levels = c("black", "royalblue", "firebrick2"))

coefs.tilde <- as.data.frame(coefs.tilde)
colnames(coefs.tilde) <- c("estG", "pG", "estGxE", "pGxE")
coefs.tilde$color <- factor(
  cols, 
  levels = c("black", "royalblue", "firebrick2"))

q <- qt(.025, df = (n-p), lower.tail = FALSE)

pdf("plot_main_vs_inter.pdf")
ggplot(coefs, aes(x = estG, y = estGxE)) +
  geom_point(aes(colour = color), size = 2) +
  geom_smooth(
    data = coefs %>% filter(color %in% c("royalblue", "firebrick2")),
    aes(colour = color),
    method = lm,
    se = FALSE,
    fullrange = TRUE,
    show.legend = FALSE) +
  scale_color_identity() +
  geom_hline(yintercept = 0, color = "gray38") + 
  geom_vline(xintercept = 0, color = "gray38") +
  geom_hline(yintercept = q, color = "gray38", linetype = "dashed") + 
  geom_vline(xintercept = q, color = "gray38", linetype = "dashed") +
  geom_hline(yintercept = -q, color = "gray38", linetype = "dashed") + 
  geom_vline(xintercept = -q, color = "gray38", linetype = "dashed") +
  xlab(expression(italic("G"))) +
  ylab(expression(italic("G x E"))) +
  annotate("text", x = -Inf, y = -q, hjust = -0.1, vjust = -0.4, size = 5,
           label = "italic(P) < 0.05", parse = TRUE, color = "gray38") +
  theme_set(theme_bw(base_size = 20)) +
  theme(axis.text = element_text(colour='black'),
        panel.border = element_rect(color = "black", linewidth = 1),
        plot.margin = unit(c(.1, .1, .1, .1), "cm"),
        plot.title = element_text(margin = margin(b = 0), size = 20)) +
  ggtitle(expression(paste("Outcome: ", italic("Y"))))
dev.off()

pdf("plot_main_vs_inter_trans.pdf")
ggplot(coefs.tilde, aes(x = estG, y = estGxE)) +
  geom_point(aes(colour = color), size = 2) +
  geom_smooth(
    data = coefs.tilde %>% filter(color %in% c("royalblue", "firebrick2")),
    aes(colour = color),
    method = lm,
    se = FALSE,
    fullrange = TRUE) +
  scale_color_identity() +
  geom_hline(yintercept = 0, color = "gray38") + 
  geom_vline(xintercept = 0, color = "gray38") +
  geom_hline(yintercept = q, color = "gray38", linetype = "dashed") + 
  geom_vline(xintercept = q, color = "gray38", linetype = "dashed") +
  geom_hline(yintercept = -q, color = "gray38", linetype = "dashed") + 
  geom_vline(xintercept = -q, color = "gray38", linetype = "dashed") +
  xlab(expression(italic("G"))) +
  ylab(expression(italic("G x E"))) +
  annotate("text", x = -Inf, y = -q, hjust = -0.1, vjust = -0.4, size = 5,
           label = "italic(P) < 0.05", parse = TRUE, color = "gray38") +
  theme_set(theme_bw(base_size = 20)) +
  theme(axis.text = element_text(colour='black'),
        panel.border = element_rect(color = "black", linewidth = 1),
        plot.margin = unit(c(.1, .1, .1, .1), "cm"),
        plot.title = element_text(margin = margin(b = 0), size = 20)) +
  ggtitle(expression(paste("Outcome: ", "log(", italic("Y"), " - min(", 
                           italic("Y"), "))")))
dev.off()
