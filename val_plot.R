#' Create a calibration plot and calculate different measures of model 
#' performance
#'
#' @param lp linear predictor of the model to be validated
#' @param y observed outcome
#' @param lim limits of the axes of the plotting area
#' @param xlab label of the x-axis
#' @param ylab label of the y-axis
#' @param main title of the plot
#' @param g number of grouped patients added to the validation plot
#' @param statloc location of the calculated performance measures
#' @param text.cex relative size of text added to the plot
#' @param d0lab label for patients without an event
#' @param d1lab label for patients with an event
#' @param xoutcome x-coordinate at which labels for patients with and without the outcome are plotted
#' @param x_leg x-coordinate of the legend
#' @param y_leg y-coordinate of the legend
#' @param conf.int How the confidence intervals around the non-parametric smoother are calculated
#' @param B number of bootstrap samples are drawn when the confidence intervals around the smoother are calculated using the bootstrap method
#' @return A validation plot of the model showing the agreement between predicted probabilities and observed outcomes.
val_plot <- function(lp, y, lim = c(0, 1), xlab = 'Predicted probability',
                     ylab = 'Observed proportion', main = '', g = 5, 
                     statloc = c(0.1, 0.7), text.cex = 0.7, d0lab = '0',
                     d1lab = '1', xoutcome = 0.9, x_leg = 0.75, y_leg = 0.5,
                     conf.int = c('delta', 'bootstrap', 'non-parametric', 'none'), 
                     B = 2000){
  library(rms)
  
  conf.int <- match.arg(conf.int)
  plot(0, 0, type = 'n', xlim = lim, ylim = c(-0.2 * (lim[2] - lim[1]), lim[2]),
       xlab = xlab, ylab = ylab, main = main, xaxs = 'i', yaxs = 'i')
  abline(0, 1, lty = 2)  
  abline(0, 0, lty = 3)
  
  p_hat <- plogis(lp)
  
  Sm <- loess(y ~ p_hat, iter = 0)
  
  # Restrict plotting to above the 1st percentile and below the 99th percentile
  bounds <- quantile(p_hat, probs = c(0.01, 0.99))
  
  sm_x <- sort(Sm$x)
  I1   <- sm_x>=bounds[1]&sm_x<=bounds[2]
  
  sm_y <- predict(Sm , se = TRUE, newdata = sm_x)
  I2   <- I1&(sm_y$fit>=0&sm_y$fit<=1)
  
  sm_x       <- sm_x[I2]
  smooth_y <- sm_y$fit[I2]
  if(conf.int=='delta'){
    # Confidence bounds for the LOESS estimate. Use delta method to ensure
    # that confidence bounds stay between 0 and 1.
    se_logit <- sm_y$se.fit/(sm_y$fit * (1 - sm_y$fit))
    sm_lower <- plogis(qlogis(sm_y$fit[I2]) - 1.96 * se_logit[I2])
    sm_upper <- plogis(qlogis(sm_y$fit[I2]) + 1.96 * se_logit[I2])
  }
  if(conf.int=='bootstrap'){
    # Confidence intervals for the LOESS estimate using bootstrap  
    boot_loess <- matrix(nrow = B, ncol = length(sm_x))
    for(j in 1:B){
      index <- sample(1:length(y), replace = TRUE)
      p_hat_b <- p_hat[index]
      y_b     <- y[index]
      Sm_b    <- loess(y_b ~ p_hat_b, iter = 0)
      boot_loess[j, ] <- predict(Sm_b, newdata = sm_x)
    }
    sm_lower <- apply(boot_loess, 2, quantile, na.rm = T, probs = 0.025)
    sm_upper <- apply(boot_loess, 2, quantile, na.rm = T, probs = 0.975)
  }
  if(conf.int=='non-parametric'){
    # Confidence intervals based on the standard error of the LOESS function
    sm_lower <- sm_y$fit - 1.96 * sm_y$se.fit
    sm_upper <- sm_y$fit + 1.96 * sm_y$se.fit
    
    sm_lower <- sm_lower[I2]
    sm_upper <- sm_upper[I2]
  }
  polygon(x = c(sm_x, rev(sm_x)), y = c(sm_upper, rev(sm_lower)), 
          col = rgb(177, 177, 177, 177, maxColorValue = 255), border = NA)
  lines(sm_x, smooth_y, lwd = 2)
  
  # Add grouped observations
  q     <- cut2(x = p_hat, levels.mean = TRUE, g = g)
  means <- as.double(levels(q))
  prop  <- tapply(y, q, function(x) mean(x, na.rm = TRUE))
  points(means, prop, pch = 2)
  
  n_group <- as.numeric(table(q))
  se_g   <- sqrt(prop * (1 - prop)/n_group)
  for(i in 1:g){
    lines(c(means[i], means[i]), c(prop[i] - 1.96 * se_g[i], prop[i] + 1.96 * se_g[i]))
  }
  
  # Calculate performance measures
  fit1 <- lrm(y ~ lp)
  fit2 <- lrm(y ~ offset(lp))
  
  stats <- paste('Calibration\n',
                 '...in the large: ', sprintf('%.2f', fit2$coefficients), '\n',
                 '...slope : ', sprintf('%.2f', fit1$coefficients[2]), '\n',
                 'Discrimination\n',
                 '...c-statistic : ', sprintf('%.2f', fit1$stats['C']), sep = '')
  text(statloc[1], statloc[2], stats, pos = 4, cex = text.cex)
  
  # Add 'spikes' showing distribution of predicted probabilities 
  bins <- seq(0, 1, length = 99)
  x <- p_hat
  x <- x[x >= 0 & x <= 1]
  
  risk_dist_level <- .9 * (par()$usr[3])/2
  
  f0	  <- table(cut(x[y==0],bins))
  f1	  <- table(cut(x[y==1],bins))
  j0	  <- f0 > 0
  j1	  <- f1 > 0
  bins0 <- (bins[-99])[j0]
  bins1 <- (bins[-99])[j1]
  f0	  <- f0[j0]
  f1	  <- f1[j1]
  maxf  <- max(f0,f1)
  f0	  <- (0.1*f0)/maxf
  f1	  <- (0.1*f1)/maxf
  
  segments(bins1, risk_dist_level, bins1, f1 + risk_dist_level)
  segments(bins0,risk_dist_level,bins0,-f0 + risk_dist_level)
  lines(c(min(bins0,bins1)-0.01,max(bins0,bins1)+0.01),c(risk_dist_level,risk_dist_level))
  text(xoutcome, risk_dist_level + 0.5 * abs(risk_dist_level), d1lab, cex = text.cex)
  text(xoutcome, risk_dist_level - 0.5 * abs(risk_dist_level), d0lab, cex = text.cex)
  
  # Add a legend to the plot
  legend(x_leg, y_leg, c('Ideal', 'Non-parametric', 'Grouped observations'), 
         lty = c(2, 1, NA), pch = c(NA, NA, 2), lwd = c(1, 2, NA), bty = 'n', 
         cex = text.cex)
}