## ----brglm2--------------------------------------------------------------
library(brglm2)
data("endometrial", package = "brglm2")
endometrialML <- glm(HG ~ NV + PI + EH, data = endometrial,family = binomial("probit"))
endometrialBR_mean <- update(endometrialML, method = "brglmFit", type = "AS_mean")

## ----AD, message = FALSE, warning = FALSE, results = FALSE---------------
library(BRTemplate)
X <- model.matrix(endometrialML)
y <- as.vector(endometrial$HG)
parain <- rep(0, ncol(X))
file.name <- "Probit.cpp"
file.remote <- system.file("TMB", file.name, package = "BRTemplate")
file.copy(file.remote, getwd())
obj <- get_AD(file.name, X, y, parain)

## ----probit gen----------------------------------------------------------
print(genProbit)

## ----est, results = FALSE------------------------------------------------
library(nleqslv)
library(parallel)
nc <- detectCores()
br.emp <- nleqslv(rep(0,ncol(X)), gsolv, ADobj = obj, dll="Probit", datagen=genProbit,
                ncores = nc, R=0, seed=NULL, y=y, observed=TRUE, X=X, 
                method="Newton", jac=jsolv)

## ----res1----------------------------------------------------------------
mat <- cbind(endometrialML$coef, endometrialBR_mean$coef, br.MC$x, br.emp$x)
colnames(mat) <- c("ML","BR (exact)", "BR (Monte Carlo)",  "BR (empirical)")
print(mat)

## ----q FS----------------------------------------------------------------
br.emp2 <- nleqslv(rep(0,ncol(X)), gsolv, ADobj = obj, dll="Probit", datagen=genProbit,
                ncores = nc, R=0, seed=NULL, y=y, observed=FALSE, X=X, 
                method="Newton", jac=jsolv)
mat <- cbind(mat, br.emp2$x)
colnames(mat) <-  c("ML","BR (exact)", "BR (Monte Carlo)",  "BR (empirical, NR)",  "BR (empirical, FS)")
print(mat)

