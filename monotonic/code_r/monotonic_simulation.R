#### ************************ Header ************************ ####
# File: monotonic_simulation.R                                   #
# Author: ian2j                                                  #
# Last Modified: 06AUG2016                                       #
# Purpose: Test an idea to force a fitted curve to be monotonic  #
# ************************************************************** #

#### Change working directory ####
setwd("C:/Users/Ian/Desktop/Master/Research/ADN/Monotonic/")

#### Set Random Seed ####
set.seed(31416)

#### Create sample data ####
N <- 1e2 # sample size
LB <- -1 # independent variable lower bound
UB <- 5 # independent variable upper bound
X <- runif(N,
           LB,
           UB) # independent variable
ALPHA <- 1 # scaling factor for dependent variable mean
BETA <- 2 # additional scaling factor for inverse logit
TAU <- .1 # dependent variable stdev
Y <- rnorm(N,
           ALPHA *
             exp(X * BETA) /
             (1 + exp(X * BETA)),
           TAU) # dependent variable

#### Setup plot options ####
# Background color of plots
BGCOLOR <- "white"
# Color of points
DOTCOLOR <- "navy"
# Color for starting curve (default regression)
STCOLOR <- "lightgreen"
# Color for final curve (our new regression)
EDCOLOR <- "purple"

# Make initial scatterplot
png("scatterplot.png",
    width = 8,
    height = 6,
    units = "in",
    res = 300)
# Set plot background color
par(bg = BGCOLOR)
# Plot points
plot(X,
     Y,
     pch = 16,
     col = DOTCOLOR)
dev.off()

#### Newton's method functions ####
# Objective function for which we seek B so that f(B) = 0
# NOTE: In our case, B consists of three values [b0, b1, b2]
#       and this function is actually the collection of partial
#       derivatives with respect to b0, b1, and b2 of our
#       loss + penalty function
f <- function(X,
              Y,
              cbeta,
              LAMBDA,
              MID) {
  d <- cbind(1,
             X,
             X ^ 2)
  common1 <- Y - d %*% cbeta
  common2 <- 2 * LAMBDA *
    (cbeta[2] / (2 * cbeta[3]) + MID) ^ (-3)
  b0 <- -2 * sum(common1)
  b1 <- -2 * sum(common1 * X)
  b1 <- b1 - common2 * (1 / (2 * cbeta[3]))
  b2 <- -2 * sum(common1 * X ^ 2)
  b2 <- b2 - common2 * (-cbeta[2] / (2 * cbeta[3] ^ 2))
  return(c(b0, 
           b1, 
           b2))
}
# Derivative of our objective function to be used in Newton's 
# method. Because our objective function is already the partial
# derivatives of the loss + penalty function with respect to 
# b0, b1, and b2, this is really the Hessian matrix of second 
# partial derivatives
fp <- function(X,
               Y,
               cbeta,
               LAMBDA,
               MID) {
  H <- matrix(0, 
              nrow = 3, 
              ncol = 3)
  common <- (cbeta[2] / (2 * cbeta[3]) + MID)
  H[1,1] <- 2 * length(X)
  H[1,2] <- 2 * sum(X)
  H[1,3] <- 2 * sum(X ^ 2)
  H[2,1] <- H[1,2]
  H[2,2] <- H[1,3] + 1.5 * 
    LAMBDA / cbeta[3] ^ 2 * common ^ (-4)
  H[2,3] <- 2 * sum(X ^ 3) + 
    LAMBDA / cbeta[3] ^ 2 * 
    (common ^ (-3) - 
       cbeta[2] * 1.5 / cbeta[3] * common ^ (-4))
  H[3,1] <- H[1,3]
  H[3,2] <- 2 * sum(X ^ 3) - LAMBDA *
    (common ^ (-3) *
       (-1 / cbeta[3] ^ 2) + 
       1.5 * cbeta[2] / cbeta[3] ^ 3 * common ^ (-4))
  H[3,3] <- 2 * sum(X ^ 4) + 
    LAMBDA * cbeta[2] * 
    (common ^ (-3) * 
       (-2 / cbeta[3] ^ 3) + 
       1.5 * cbeta[2] / cbeta[3] ^ 4 * common ^ (-4))
  return(H)
}
# Basic implementation of Newton's method for our problem
solver <- function(X, 
                   Y, 
                   LAMBDA, 
                   MID = mean(range(X)),
                   MAXITER = 100, 
                   TOL = 1e-4) {
  # Use a standard least squares regression model 
  # to select an initial set of coefficients
  d <- cbind(X, 
             X ^ 2)
  ibeta <- as.numeric(coef(lm(Y ~ d)))
  # Initialize variables for the algorithm
  ITER <- 1
  NORM <- 1e4
  while(ITER < MAXITER & NORM > TOL) {
    # While we have not crossed the maximum number 
    # of iterations AND remain above the tolerance 
    # for convergence,
    
    # Calculate the current value of the 
    # objective function
    fval <- f(X,
              Y,
              ibeta,
              LAMBDA,
              MID)
    
    # Calculate how far we are away from the 
    # optimal values of [0, 0, 0]
    NORM <- sum(abs(fval))
    
    # Calculate the Hessian matrix
    fpval <- fp(X,
                Y,
                ibeta,
                LAMBDA,
                MID)
    
    # Update the values of the coefficients
    # using Newton's method
    ibeta <- ibeta - solve(fpval, fval)
    
    # Increment the iteration count
    ITER <- ITER + 1
  }
  # Return the final solution
  return(ibeta)
}

#### Experiments ####
# Set up a sequence of penalty terms
LAMSEQ <- seq(from = 0,
              to = 5,
              len = 25)
# Get the number of lambda values to try
NL <- length(LAMSEQ) 
# Initialize matrix for output
BMAT <- matrix(0, 
               nrow = NL,
               ncol = 3)
# Loop through and calculate estimated coefficients using 
# Newton's method for each value of lambda
for(l in 1:NL) {
  BMAT[l,] <- solver(X, 
                     Y, 
                     LAMBDA = LAMSEQ[l], 
                     MID = mean(range(X)),
                     MAXITER = 20, 
                     TOL = 1e-4)
}

# Overlay all solutions on the initial scatterplot
png("all_solutions.png",
    width = 8,
    height = 6,
    units = "in",
    res = 300)
COLSEQ <- colorRampPalette(c(STCOLOR,
                             EDCOLOR))(NL)
par(bg = BGCOLOR)
plot(X, 
     Y,
     pch = 16,
     col = DOTCOLOR)
for(l in 1:NL) {
  curve(BMAT[l, 1] + 
          BMAT[l, 2] * x + 
          BMAT[l, 3] * x ^ 2,
        from = LB,
        to = UB,
        add = TRUE,
        lwd = 2,
        col = COLSEQ[l])
}
# Emphasize the initial (unpenalized) solution
curve(BMAT[1, 1] + 
        BMAT[1, 2] * x + 
        BMAT[1, 3] * x ^ 2,
      from = LB,
      to = UB,
      add = TRUE,
      lwd = 4,
      col = STCOLOR)
# Emphasize the final (most penalized) solution
curve(BMAT[NL,1] + 
        BMAT[NL,2] * x + 
        BMAT[NL,3] * x ^ 2,
      from = LB,
      to = UB,
      add = TRUE,
      lwd = 4,
      col = EDCOLOR)
# Plot the points again to put them on top of the lines
points(X, 
       Y, 
       pch = 16, 
       col = DOTCOLOR)
dev.off()

# Plot only the initial curve
png("initial_results.png",
    width = 8,
    height = 6,
    units = "in",
    res = 300)
# Set plot background color
par(bg = BGCOLOR)
# Plot the data
plot(X,
     Y,
     pch = 16,
     col = DOTCOLOR)
# Add the initial curve
curve(BMAT[1, 1] + 
        BMAT[1, 2] * x + 
        BMAT[1, 3] * x ^ 2,
      from = LB,
      to = UB,
      add = TRUE,
      lwd = 4,
      col = STCOLOR)
# Plot the points once more
points(X, 
       Y, 
       pch = 16, 
       col = DOTCOLOR)
dev.off()

# Plot only the initial and final curves
png("simple_results.png",
    width = 8,
    height = 6,
    units = "in",
    res = 300)
# Set plot background color
par(bg = BGCOLOR)
# Plot the data
plot(X,
     Y,
     pch = 16,
     col = DOTCOLOR)
# Add the initial curve
curve(BMAT[1, 1] + 
        BMAT[1, 2] * x +
        BMAT[1, 3] * x ^ 2,
      from = LB,
      to = UB,
      add = TRUE,
      lwd = 4,
      col = STCOLOR)
# Add the final curve
curve(BMAT[NL, 1] + 
        BMAT[NL, 2] * x + 
        BMAT[NL, 3] * x ^ 2,
      from = LB,
      to = UB,
      add = TRUE,
      lwd = 4,
      col = EDCOLOR)
# Plot the points once more
points(X, 
       Y, 
       pch = 16, 
       col = DOTCOLOR)
# Add a legend
legend("bottomright",
       legend = c(expression(lambda * " = 0"),
                  expression(lambda * " = 5")),
       lwd = c(4, 4),
       lty = c(1, 1),
       col = c(STCOLOR, 
               EDCOLOR))
dev.off()
