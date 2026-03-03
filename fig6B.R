library(ggplot2)

h2 <- .5
s2 <- 1
n0 <- 5000
n1 <- 5000
n <- n0 + n1
p <- 200
frac.causal.gxe <- .5
frac.positive <- .5
alpha <- 2
pthr <- .05

phi.x <- function(x) x
phi.log <- function(x) log(x - min(x) + 1e-6)
funcs <- list(phi.x, phi.log)

res <- list()

for (s2.gxe in seq(0, 1.5, .1)) {
  print(s2.gxe)
  res.iter <- list()
  
  for (i in 1:100) { # Reduce the number of replications to speed up the simulation
    b <- abs(rnorm(p, sd = sqrt(s2/p)))
    
    p.gxe <- round(p*frac.causal.gxe)
    b.gxe <- -abs(rnorm(p.gxe, sd = sqrt(s2.gxe/p.gxe)))
    n.positive <- round(p.gxe*frac.positive)
    
    if (n.positive > 0) {
      b.gxe[1:n.positive] <- -1*b.gxe[1:n.positive]
    }
    
    b.tmp <- rep(0, p)
    b.tmp[1:p.gxe] <- b.gxe
    b.gxe <- b.tmp
    
    freq <- runif(p, min = .1, max = .5)
    X <- sapply(freq, function(x) rbinom(n, size = 2, prob = x))
    X <- scale(X)
    E <- rep(c(0, 1), c(n0, n1))
    E <- c(scale(E))
    
    sg <- var(X %*% b)
    se <- sg*((1 - h2)/h2)
    
    y <- alpha*E + X %*% b + (X*E) %*% b.gxe + rnorm(n, sd = sqrt(se))
    ys <- matrix(sapply(funcs, function(f) f(y)), ncol = length(funcs))
    
    for (j in 1:ncol(ys)) {
      coefs <- lapply(1:p, 
                      function(ind) coef(summary(lm(ys[, j] ~ X[, ind]*E))))
      coefs <- do.call(
        "rbind", 
        lapply(coefs, function(x) c(x[2, c(3, 4)], x[4, c(3, 4)])))
      
      coefs <- cbind(coefs, sign(coefs[, 1]*coefs[, 3]))
      coefs.signif <- matrix(coefs[coefs[, 4] < pthr,], ncol = 5)
      coefs.sc <- matrix(coefs.signif[coefs.signif[, 2] < pthr,], ncol = 5)
      
      if (nrow(coefs.sc) > 1) {
        sc = max(sum(coefs.sc[, 5] == 1), sum(coefs.sc[, 5] == -1))/
          nrow(coefs.sc)
      } else {
        sc = NA
      }
      
      res.iter[[length(res.iter) + 1]] <- c(nrow(coefs.signif), sc, j)
    }
  }
  
  res.iter <- do.call("rbind", res.iter)
  
  for (i in unique(res.iter[, 3])) {
    stats <- c(apply(
      res.iter[res.iter[, 3] == i, 1:2], 
      2, 
      function(x) c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))))
    res[[length(res) + 1]] <- c(stats, i, s2.gxe)
  }
}

res <- as.data.frame(do.call("rbind", res))
colnames(res) <- c("detected", "SDdetected", "sc", "SDsc", "trans", "s2GxE")
res$trans <- as.factor(res$trans)

pdf("plot_sign_consistency_rate.pdf")
ggplot(res, aes(x = s2GxE, y = sc, color = trans)) +
  geom_point(size = 2) + 
  geom_line(linewidth = 0.5) +
  geom_errorbar(aes(ymin = sc-SDsc, ymax = sc+SDsc), width = .04,
                linetype = "solid") +
  scale_color_manual(values = c("gray45", "black")) +
  xlab(expression(sigma[GxE]^2)) +
  ylab("Sign consistency rate") +
  theme_set(theme_bw(base_size = 20)) +
  theme(axis.text = element_text(colour='black'),
        panel.border = element_rect(color = "black", linewidth = 1),
        legend.position = "none",
        plot.margin = unit(c(.1, .1, .1, .1), "cm"))
dev.off()

pdf("plot_sign_consistency_n.pdf")
ggplot(res, aes(x = s2GxE, y = detected, color = trans)) +
  geom_point(size = 2) + 
  geom_line(linewidth = 0.5) +
  geom_errorbar(aes(ymin = detected-SDdetected, ymax = detected+SDdetected), 
                width = .04,
                linetype = "solid") +
  scale_color_manual(values = c("gray45", "black")) +
  xlab(expression(sigma[GxE]^2)) + 
  ylab(expression(paste("Number of detected ", italic("G x E")))) + 
  theme_set(theme_bw(base_size = 20)) +
  theme(axis.text = element_text(colour='black'),
        panel.border = element_rect(color = "black", linewidth = 1),
        legend.position = "none",
        plot.margin = unit(c(.1, .1, .1, .1), "cm"))
dev.off()
