## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# install.packages("FMM")

## ----eval=FALSE---------------------------------------------------------------
# install.packages("devtools")
# devtools::install_github("FMMGroupVa/FMM")

## -----------------------------------------------------------------------------
library("FMM")

## -----------------------------------------------------------------------------
generateFMM(M=2,A=3,alpha=1.5,beta=2.3,omega=0.1,
            plot=TRUE,outvalues=FALSE)

## -----------------------------------------------------------------------------
fmm2.data <-generateFMM(M=0,A=c(2,2),alpha=c(1.5,3.4),beta=c(0.2,2.3),omega=c(0.1,0.2),
                         plot=FALSE, outvalues=TRUE)
str(fmm2.data)

## -----------------------------------------------------------------------------
set.seed(15)
fmm2.data <-generateFMM(M=0,A=c(2,2),alpha=c(1.5,3.4),beta=c(0.2,2.3),omega=c(0.1,0.2),
                        plot=TRUE, outvalues=TRUE,
                        sigmaNoise=0.3)

## -----------------------------------------------------------------------------
fit.fmm2 <- fitFMM(vData = fmm2.data$y, timePoints = fmm2.data$t, nback = 2)
summary(fit.fmm2)

## -----------------------------------------------------------------------------
getFMMPeaks(fit.fmm2, timePointsIn2pi = TRUE)

## -----------------------------------------------------------------------------
fit1 <- fitFMM(vData = fmm2.data$y, timePoints = fmm2.data$t, nback = 2, 
               lengthAlphaGrid = 48, lengthOmegaGrid = 24, showTime = TRUE)
fit2 <- fitFMM(vData = fmm2.data$y, timePoints = fmm2.data$t, nback = 2, 
               lengthAlphaGrid = 14, lengthOmegaGrid = 7, showTime = TRUE)
getR2(fit1)
getR2(fit2)


## -----------------------------------------------------------------------------
fit3 <- fitFMM(vData = fmm2.data$y, timePoints = fmm2.data$t, nback = 2, 
               maxiter = 5, stopFunction = R2(difMax = 0.01),
               showTime = TRUE, showProgress = TRUE)

## -----------------------------------------------------------------------------
set.seed(1115)
rfmm.data <-generateFMM(M = 3, A = c(4,3,1.5,1), 
                        alpha = c(3.8,1.2,4.5,2),
                        beta = c(rep(3,2),rep(1,2)),
                        omega = c(rep(0.1,2),rep(0.05,2)),
                        plot = TRUE, outvalues = TRUE,
                        sigmaNoise = 0.3)

## -----------------------------------------------------------------------------
fit.rfmm <- fitFMM(vData = rfmm.data$y, timePoints = rfmm.data$t, nback = 4,
                   betaOmegaRestrictions = c(1, 1, 2, 2),
                   lengthAlphaGrid = 24, lengthOmegaGrid = 15, numReps = 5)
summary(fit.rfmm)

## ----message=FALSE, fig.height=5, fig.width=10--------------------------------
titleText <- "Two FMM waves"
par(mfrow=c(1,2))
# default plot
plotFMM(fit.fmm2, textExtra = titleText)
# component plot
plotFMM(fit.fmm2, components = TRUE, textExtra = titleText, legendInComponentsPlot = TRUE)

## ----message=FALSE, fig.height=5, fig.width=10--------------------------------
library("RColorBrewer")
library("ggplot2")
library("gridExtra")

titleText <- "Four restricted FMM waves"
# default plot
defaultrFMM2 <- plotFMM(fit.rfmm, use_ggplot2 = TRUE, textExtra = titleText)
defaultrFMM2 <- defaultrFMM2 + theme(plot.margin=unit(c(1,0.25,1.3,1), "cm")) +
  ylim(-5, 6)
# component plot
comprFMM2 <- plotFMM(fit.rfmm, components=TRUE, use_ggplot2 = TRUE, textExtra = titleText)
comprFMM2 <- comprFMM2 + theme(plot.margin=unit(c(1,0.25,0,1), "cm")) +
  ylim(-5, 6) + scale_color_manual(values = brewer.pal("Set1",n = 8)[3:6])

grid.arrange(defaultrFMM2, comprFMM2, nrow = 1)


