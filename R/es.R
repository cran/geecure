es <- function(Time, Status, X, Z, id, model, corstr, stdz, esmax, eps) {
    Kn <- length(id)
    K <- length(unique(id))
    n <- as.vector(table(id))
    Z1 <- Z
    X1 <- X
    if (stdz) {
        for (i in 2:ncol(Z1)) {
            Z[, i] <- (Z[, i] - mean(Z[, i]))/sd(Z[, i])
        }
        if (model == "para") {
            for (i in 2:ncol(X1)) {
                X[, i] <- (X[, i] - mean(X[, i]))/sd(X[, i])
            }
        }
        if (model == "semi") {
            for (i in 1:ncol(X1)) {
                X[, i] <- (X[, i] - mean(X[, i]))/sd(X[, i])
            }
        }
    }
    ppmt2 <- c(rep(0, dim(Z)[2]), rep(0, dim(X)[2]))
    w <- Status
    Lambda <- initial_Lambda(Time, Status, X, Z, id, model, corstr = "independence")$Lambda
    KK1 <- 1
    repeat {
        gamma1 <- eval(parse(text = paste("geese", "(", "w~Z[,-1]", ",id=", "id", ",family = binomial", ",corstr='", 
            corstr, "'", ")", sep = "")))$beta
        # gamma1<- geega(w,Z,gamma=rep(0,dim(Z)[2]),id,corstr)$gamma
        SK2 <- 1
        beta2 <- rep(0, dim(X)[2])
        Lambda <- Lambda
        repeat {
            beta1 <- geebt(Status, Lambda, X, beta = rep(0, dim(X)[2]), w, id, corstr)$beta
            gSS3 <- basesurv(Time, Status, X = X1, beta = beta1, Lambda, w, model)$bcumhaz
            if ((any(abs(beta1 - beta2) > 1e-06)) & (SK2 <= 20)) {
                beta2 <- beta1
                Lambda <- gSS3
                SK2 <- SK2 + 1
            } else {
                ppmt1 <- c(gamma1, beta1)
                survival <- basesurv(Time, Status, X = X1, beta = beta1, Lambda, w, model)$uncuresurv
                w <- Status + ((1 - Status) * exp(gamma1 %*% t(Z1)) * survival)/(1 + exp(gamma1 %*% t(Z1)) * survival)
                w <- as.vector(w)
                break
            }
        }
        if (any(abs(ppmt1 - ppmt2) > eps) && (KK1 < esmax)) {
            ppmt2 <- ppmt1
            KK1 <- KK1 + 1
        } else break
    }
    kappa <- basesurv(Time, Status, X = X1, beta = beta1, Lambda, w, model)$kappa
    if (stdz == "TRUE") {
        gamma1 <- c(gamma1[1] - sum((gamma1[-1] * apply(Z1[, -1, drop = FALSE], 2, mean)/apply(Z1[, -1, drop = FALSE], 2, sd))), gamma1[-1]/apply(Z1[, 
            -1, drop = FALSE], 2, sd))
        beta1 <- beta1/apply(X1, 2, sd)
    }
    if (corstr == "independence") {
        alpha <- 0
        phi <- 1
    }
    if (corstr == "exchangeable") {
        alpha <- eval(parse(text = paste("geese", "(", "w~Z[,-1]", ",id=", "id", ",family = binomial", ",corstr='", 
            corstr, "'", ")", sep = "")))$alpha
        phi <- eval(parse(text = paste("geese", "(", "w~Z[,-1]", ",id=", "id", ",family = binomial", ",corstr='", 
            corstr, "'", ")", sep = "")))$gamma
    }
    rho <- geebt(Status, Lambda, X, beta = beta1, w, id, corstr)$rho
    pphi <- geebt(Status, Lambda, X, beta = beta1, w, id, corstr)$pphi
    uncureprob <- as.vector(exp(gamma1 %*% t(Z1))/(1 + exp(gamma1 %*% t(Z1))))
    convergence <- sum((ppmt1 - ppmt2)^2)
    es <- list(gamma = gamma1, beta = beta1, kappa = kappa, Lambda = Lambda, gcor = alpha, gphi = phi, bcor = rho, 
        bphi = pphi, w = w, Survival = survival, Uncureprob = uncureprob, tau = convergence)
}