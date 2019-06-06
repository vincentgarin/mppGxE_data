################################################################################
################################################################################
#
# The usefulness of multi-parent multi-environment QTL analysis: an illustration
# in different NAM populations
#
# Vincent Garin, Marcos Malosetti, Fred van Eeuwijk
#
################################################################################
################################################################################

# Library

# you can find in the ~mppGxE_data/software a tar.gz version of the packages
# mppR (1.1.11) and mppGxE (1.0.0) used in this research


library(mppR)
library(mppGxE)
library(asreml)

setwd('E:/mppGxE_data') # set your working directory here

# create folder to store the results

res_fold <- file.path(getwd(),'results')
res_EUNAM <- file.path(res_fold, 'EUNAM')
res_USNAM <- file.path(res_fold, 'USNAM')

dir.create(path = res_fold)
dir.create(path = res_EUNAM)
dir.create(path = res_USNAM)

######################
# 1. EU-NAM analysis #
######################

load(file = './data/EUNAM_mppData.RData')
load(file = './data/EUNAM_plot_data.RData')

trait <- "DMY"
E1 <- "CIAM"
E2 <- "TUM"

Q_eff <- c("par", "anc", "biall")

res_fold <- "./results/EUNAM"

t_E1 <- paste(trait, E1, sep = "_")
t_E2 <- paste(trait, E2, sep = "_")
t_E1E2 <- paste(trait, paste(c(E1, E2), collapse = "_"), sep = "_")


# 1.1 Complete data QTL detection
#################################

# create some folder to save the results

res_fold_i <- file.path(res_fold, paste0('QTL_an_', t_E1E2))
res_fold_i_M1 <- file.path(res_fold_i, paste0("M", 1))
res_fold_i_M2 <- file.path(res_fold_i, paste0("M", 2))
res_fold_i_M3 <- file.path(res_fold_i, paste0("M", 3))
res_fold_i_M4 <- file.path(res_fold_i, paste0("M", 4))

dir.create(path = res_fold_i)
dir.create(path = res_fold_i_M1)
dir.create(path = res_fold_i_M2)
dir.create(path = res_fold_i_M3)
dir.create(path = res_fold_i_M4)

##### M1

for(i in 1:3){

  QTL <- mpp_proc(pop.name = "EUNAM", trait = t_E1E2,
                  trait.name =  t_E1E2, mppData = mppData,
                  Q.eff = Q_eff[i], VCOV = 'cr.err', win.cof = 50,
                  plot.gen.eff = (Q_eff[i] != 'biall'),
                  thre.cof = 4, thre.QTL = 4, CI = FALSE,
                  verbose = FALSE, output.loc = res_fold_i_M1)

}

#### M2 E1

for(i in 1:3){

  QTL <- mpp_proc(pop.name = "EUNAM", trait = t_E1,
                  trait.name =  t_E1, mppData = mppData,
                  Q.eff = Q_eff[i], VCOV = "cr.err", win.cof = 50,
                  plot.gen.eff = (Q_eff[i] != 'biall'),
                  thre.cof = 4, thre.QTL = 4, CI = FALSE,
                  verbose = FALSE, output.loc = res_fold_i_M2)

}

#### M2 E2

for(i in 1:3){

QTL <- mpp_proc(pop.name = "EUNAM", trait = t_E2,
                trait.name =  t_E2, mppData = mppData,
                Q.eff = Q_eff[i], VCOV = "cr.err", win.cof = 50,
                plot.gen.eff = (Q_eff[i] != 'biall'),
                thre.cof = 4, thre.QTL = 4, CI = FALSE,
                verbose = FALSE, output.loc = res_fold_i_M2)

}

#### M3

for(i in 1:3){

QTL <- mppGE_proc(pop.name = 'EUNAM', trait.name = "DMY",
                  mppData = mppData, trait = c(t_E1, t_E2),
                  EnvNames = c(E1, E2), Q.eff = Q_eff[i],
                  VCOV = "CS_CSRT", plot.gen.eff = TRUE,
                  output.loc = res_fold_i_M3, verbose = FALSE)

}

#### M4

for(i in 1:3){

QTL <- mppGE_oneS_proc(pop.name = 'EUNAM', trait.name = "DMY",
                      plot_data = plot_data, mppData = mppData,
                      trait = trait, EnvNames = c(E1, E2),
                      exp_des_form = 'env:Rep + env:Rep:Block',
                      Q.eff = Q_eff[i],
                      VCOV = "CS_CSRT", plot.gen.eff = TRUE,
                      verbose = FALSE, output.loc = res_fold_i_M4,
                      workspace = 20e7)
}

# 1.2 Cross-validation
######################

CV_rep <- 3
k_fold <- 3
VCOV <- 'cr.err'
VCOV2 <- 'CS_CSRT'

# create some folder to save the results

res_fold_i <- file.path(res_fold, paste0('QTL_CV_', t_E1E2))
res_fold_i_M1 <- file.path(res_fold_i, paste0("M", 1))
res_fold_i_M2 <- file.path(res_fold_i, paste0("M", 2))
res_fold_i_M3 <- file.path(res_fold_i, paste0("M", 3))
res_fold_i_M4 <- file.path(res_fold_i, paste0("M", 4))

dir.create(path = res_fold_i)
dir.create(path = res_fold_i_M1)
dir.create(path = res_fold_i_M2)
dir.create(path = res_fold_i_M3)
dir.create(path = res_fold_i_M4)

##### M1

for(j in 1:3){

  CV <- mpp_CV(pop.name = 'EUNAM', trait.name = t_E1E2,
               mppData = mppData, trait = t_E1E2, cv.ref = c(t_E1, t_E2),
               Rep = CV_rep, k = k_fold, Q.eff = Q_eff[j],
               VCOV = VCOV, thre.cof = 4, her = 1, win.cof = 50,
               window = 20, thre.QTL = 4, win.QTL = 20, alpha.bk = 0.01,
               output.loc = res_fold_i_M1, verbose = FALSE)

}

#### M2 E1

for(j in 1:3){

  CV <- mpp_CV(pop.name = 'EUNAM', trait.name = t_E1,
               mppData = mppData, trait = t_E1, cv.ref = NULL, Rep = CV_rep,
               k = k_fold, Q.eff = Q_eff[j], VCOV = VCOV,
               thre.cof = 4, her = 1, win.cof = 50, window = 20,
               thre.QTL = 4, win.QTL = 20, alpha.bk = 0.01,
               output.loc = res_fold_i_M2, verbose = FALSE)

}

#### M2 E2

for(j in 1:3){

  CV <- mpp_CV(pop.name = 'EUNAM', trait.name = t_E2,
               mppData = mppData, trait = t_E2, cv.ref = NULL, Rep = CV_rep,
               k = k_fold, Q.eff = Q_eff[j], VCOV = VCOV,
               thre.cof = 4, her = 1, win.cof = 50, window = 20,
               thre.QTL = 4, win.QTL = 20, alpha.bk = 0.01,
               output.loc = res_fold_i_M2, verbose = FALSE)

}

#### M3 (Very long - several days)

for(j in 1:3){

  CV <- mppGE_CV(pop.name = 'EUNAM', trait.name = trait,
                 mppData = mppData, trait = c(t_E1, t_E2), cv.ref = NULL,
                 Rep = CV_rep, k = k_fold, EnvNames = c(E1, E2), Q.eff = Q_eff[j],
                 VCOV = VCOV2, output.loc = res_fold_i_M3, verbose = FALSE,
                 workspace = 20e7)

}

#### M4 (Very long - several days)

for(j in 1:3){

  CV <- mppGE_oneS_CV(pop.name = 'EUNAM', trait.name = trait,
                      plot_data = plot_data, mppData = mppData,
                      trait = trait, cv.ref = c(t_E1, t_E2), Rep = CV_rep,
                      k = k_fold, EnvNames = c(E1, E2), Q.eff = Q_eff[j],
                      VCOV = VCOV2, exp_des_form = 'env:Rep + env:Rep:Block',
                      verbose = TRUE, output.loc = res_fold_i_M4,
                      workspace = 50e7)
}

######################
# 2. US-NAM analysis #
######################

load(file = './data/USNAM_mppData.RData')
load(file = './data/USNAM_plot_data.RData')

trait <- "DTA"
E1 <- "NC07"
E2 <- "NY07"

Q_eff <- "par"

res_fold <- "./results/USNAM"

t_E1 <- paste(trait, E1, sep = "_")
t_E2 <- paste(trait, E2, sep = "_")
t_E1E2 <- paste(trait, paste(c(E1, E2), collapse = "_"), sep = "_")

# 2.1 Complete data QTL detection
#################################

res_fold_i <- file.path(res_fold, paste0('QTL_an_', t_E1E2))
res_fold_i_M1 <- file.path(res_fold_i, paste0("M", 1))
res_fold_i_M2 <- file.path(res_fold_i, paste0("M", 2))
res_fold_i_M3 <- file.path(res_fold_i, paste0("M", 3))
res_fold_i_M4 <- file.path(res_fold_i, paste0("M", 4))

dir.create(path = res_fold_i)
dir.create(path = res_fold_i_M1)
dir.create(path = res_fold_i_M2)
dir.create(path = res_fold_i_M3)
dir.create(path = res_fold_i_M4)

##### M1

QTL <- mpp_proc(pop.name = "USNAM", trait = t_E1E2,
                trait.name =  t_E1E2, mppData = mppData,
                Q.eff = Q_eff, VCOV = "cr.err", win.cof = 50,
                plot.gen.eff = TRUE, thre.cof = 4, thre.QTL = 4,
                CI = FALSE, verbose = FALSE, output.loc = res_fold_i_M1)

#### M2 E1

QTL <- mpp_proc(pop.name = "USNAM", trait = t_E1,
                trait.name =  t_E1, mppData = mppData,
                Q.eff = Q_eff, VCOV = "cr.err", win.cof = 50,
                plot.gen.eff = TRUE, thre.cof = 4, thre.QTL = 4,
                CI = FALSE, verbose = FALSE, output.loc = res_fold_i_M2)

#### M2 E2

QTL <- mpp_proc(pop.name = "USNAM", trait = t_E2,
                trait.name =  t_E2, mppData = mppData,
                Q.eff = Q_eff, VCOV = "cr.err", win.cof = 50,
                plot.gen.eff = TRUE, thre.cof = 4, thre.QTL = 4,
                CI = FALSE, verbose = FALSE, output.loc = res_fold_i_M2)

#### M3

QTL <- mppGE_proc(pop.name = 'USNAM', trait.name = "DTA",
                 mppData = mppData, trait = c(t_E1, t_E2),
                 EnvNames = c(E1, E2), Q.eff = Q_eff,
                 VCOV = "CS_CSRT", plot.gen.eff = TRUE,
                 output.loc = res_fold_i_M3, verbose = FALSE,
                 workspace = 50e+7)

#### M4

e_form <- 'env:set + env:set:block + env:row + env:col'

QTL <- mppGE_oneS_proc(pop.name = 'USNAM', trait.name = "DTA",
                     plot_data = plot_data, mppData = mppData,
                     trait = trait, EnvNames = c(E1, E2),
                     Q.eff = Q_eff, VCOV = "CS_CSRT", exp_des_form = e_form,
                     plot.gen.eff = TRUE, verbose = FALSE,
                     output.loc = res_fold_i_M4, workspace = 50e+7)

# 2.2 Cross-validation
######################

CV_rep <- 3
k_fold <- 3
VCOV <- 'cr.err'
VCOV2 <- 'CS_CSRT'

res_fold_i <- file.path(res_fold, paste0('QTL_CV_', t_E1E2))
res_fold_i_M1 <- file.path(res_fold_i, paste0("M", 1))
res_fold_i_M2 <- file.path(res_fold_i, paste0("M", 2))
res_fold_i_M3 <- file.path(res_fold_i, paste0("M", 3))
res_fold_i_M4 <- file.path(res_fold_i, paste0("M", 4))

dir.create(path = res_fold_i)
dir.create(path = res_fold_i_M1)
dir.create(path = res_fold_i_M2)
dir.create(path = res_fold_i_M3)
dir.create(path = res_fold_i_M4)

##### M1

CV <- mpp_CV(pop.name = 'USNAM', trait.name = t_E1E2,
             mppData = mppData, trait = t_E1E2, cv.ref = c(t_E1, t_E2),
             Rep = CV_rep, k = k_fold, Q.eff = 'par', VCOV = VCOV,
             thre.cof = 4, her = 1, win.cof = 50, window = 20,
             thre.QTL = 4, win.QTL = 20, alpha.bk = 0.01,
             output.loc = res_fold_i_M1, verbose = FALSE)

#### M2 E1

CV <- mpp_CV(pop.name = 'USNAM', trait.name = t_E1,
             mppData = mppData, trait = t_E1, cv.ref = NULL, Rep = CV_rep,
             k = k_fold, Q.eff = 'par', VCOV = VCOV,
             thre.cof = 4, her = 1, win.cof = 50, window = 20,
             thre.QTL = 4, win.QTL = 20, alpha.bk = 0.01,
             output.loc = res_fold_i_M2, verbose = FALSE)

#### M2 E2

CV <- mpp_CV(pop.name = 'USNAM', trait.name = t_E2,
             mppData = mppData, trait = t_E2, cv.ref = NULL, Rep = CV_rep,
             k = k_fold, Q.eff = 'par', VCOV = VCOV,
             thre.cof = 4, her = 1, win.cof = 50, window = 20,
             thre.QTL = 4, win.QTL = 20, alpha.bk = 0.01,
             output.loc = res_fold_i_M2, verbose = FALSE)

#### M3 (Very long - several days)

CV <- mppGE_CV(pop.name = 'USNAM', trait.name = "DTA",
               mppData = mppData, trait = c(t_E1, t_E2), cv.ref = NULL,
               Rep = CV_rep, k = k_fold, EnvNames = c(E1, E2), Q.eff = 'par',
               VCOV = VCOV2, output.loc = res_fold_i_M3, verbose = FALSE,
               workspace = 20e7)

#### M4 (Extremely long - several days)

e_form <- 'env:set + env:set:block + env:row + env:col'

CV <- mppGE_oneS_CV(pop.name = 'USNAM', trait.name = "DTA",
                    plot_data = plot_data, mppData = mppData,
                    trait = trait, cv.ref = c(t_E1, t_E2), Rep = CV_rep,
                    k = k_fold, EnvNames = c(E1, E2), Q.eff ='par',
                    VCOV = VCOV2, exp_des_form = e_form, verbose = FALSE,
                    output.loc = res_fold_i_M4, workspace = 50e7)
