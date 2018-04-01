emes <- function(Time, Status, X, Z, id, corstr, stdz, esmax, eps) {
    Kn <- length(id)
    K <- length(unique(id))
    n <- as.vector(table(id))
    Z1 <- Z
    X1 <- X
    t11 <- sort(Time)
    c11 <- Status[order(Time)]
    kk <- length(table(t11[c11 == 1]))
    gg1 <- Status
    gg2 <- log(Time)
    gg1[gg1 < 1e-06] <- 1e-06
    gg3 <- log(gg1) + log(Time)
    if (stdz) {
        for (i in 2:ncol(Z1)) {
            Z[, i] <- (Z[, i] - mean(Z[, i]))/sd(Z[, i])
        }
        for (i in 1:ncol(X1)) {
            X[, i] <- (X[, i] - mean(X[, i]))/sd(X[, i])
        }
    }
   
    ww1 <- eval(parse(text = paste("geese", "(", "gg1 ~ Z[,-1]", ", id = ", "id", ", family = binomial", ", corstr = '", 
            corstr, "'", ")", sep = "")))
    ww2 <- eval(parse(text = paste("geese", "(", "Status ~ X + offset(gg3) - 1", ", id = ", "id", ", family = poisson", ", corstr = '", 
            corstr, "'", ")", sep = "")))

    pmt1 <- c(ww1$beta, ww2$beta)
   
    KK1 <- 1

    repeat {
       gSSS1 <- rep(0, kk)
       repeat{  gSS1 <- basesurv(Time, Status, X, beta = pmt1[(1+dim(Z)[2]):(dim(X)[2]+dim(Z)[2])], w = gg1, model = "semi")$bsurv
                gSS2 <- basesurv(Time, Status, X, beta = pmt1[(1+dim(Z)[2]):(dim(X)[2]+dim(Z)[2])], w = gg1, model = "semi")$uncuresurv
                gSS3 <- basesurv(Time, Status, X, beta = pmt1[(1+dim(Z)[2]):(dim(X)[2]+dim(Z)[2])], w = gg1, model = "semi")$bcumhaz
                gg2 <- log(gSS3)
                gg3 <- as.vector(log(gg1)) + gg2
                
                ww2 <- eval(parse(text = paste("geese", "(", "Status ~ X + offset(gg3) - 1", ", id = ", "id", ", family = poisson", ", corstr = '", 
                corstr, "'", ")", sep = "")))
                if (any(abs(ww2$beta - pmt1[(1+dim(Z)[2]):(dim(X)[2]+dim(Z)[2])]) > 1e-06) || any(abs(gSS1 - gSSS1) > 1e-06)) {
                pmt1[(1+dim(Z)[2]):(dim(X)[2]+dim(Z)[2])] <- ww2$beta
                gSSS1 <- gSS1
                } else {
                gg1 <- as.vector(Status + ((1 - Status) * exp(pmt1[1:dim(Z)[2]] %*% t(Z)) * gSS2)/(1 + exp(pmt1[1:dim(Z)[2]] %*% t(Z)) * gSS2))
                g11 <- gg1[order(Time)]
                gg1[gg1 < 1e-06] <- 1e-06
                gg3 <- as.vector(log(gg1)) + gg2
                break
               }
            }
        ww1 <- eval(parse(text = paste("geese", "(", "gg1 ~ Z[,-1]", ", id = ", "id", ", family = binomial", ", corstr = '", 
            corstr, "'", ")", sep = "")))
        ww2 <- eval(parse(text = paste("geese", "(", "Status ~ X + offset(gg3) - 1", ", id = ", "id", ", family = poisson", ", corstr = '", 
            corstr, "'", ")", sep = "")))
        
        pmt2 <- c(ww1$beta, ww2$beta)
        if (any(abs(pmt2 - pmt1) > eps)) {
            pmt1 <- pmt2
            KK1 <- KK1 + 1
        } else break
    }
        gamma1 <- pmt1[1:dim(Z)[2]]
        beta1 <- pmt1[(1+dim(Z)[2]):(dim(X)[2]+dim(Z)[2])]

       if (stdz == "TRUE") {
        gamma1 <- c(gamma1[1] - sum((gamma1[-1] * apply(Z1[, -1, drop = FALSE], 2, mean)/apply(Z1[, -1, drop = FALSE], 2, sd))), gamma1[-1]/apply(Z1[, 
            -1, drop = FALSE], 2, sd))
        beta1 <- beta1/apply(X1, 2, sd)
    }
    if (corstr == "independence") {
        alpha <- 0
        phi <- 1
        rho <- 0
        pphi <- 1
    }
    if (corstr == "exchangeable") {
        alpha <- eval(parse(text = paste("geese", "(", "gg1 ~ Z[,-1]", ", id = ", "id", ", family = binomial", ", corstr= '", 
            corstr, "'", ")", sep = "")))$alpha
        phi <- eval(parse(text = paste("geese", "(", "gg1 ~ Z[,-1]", ", id = ", "id", ", family = binomial", ", corstr= '", 
            corstr, "'", ")", sep = "")))$gamma
        rho <- eval(parse(text = paste("geese", "(", "Status ~ X + offset(gg3) - 1", ",id=", "id", ",family = poisson", ",corstr='", 
            corstr, "'", ")", sep = "")))$alpha
        pphi <- eval(parse(text = paste("geese", "(", "Status ~ X + offset(gg3) - 1", ",id=", "id", ",family = poisson", ",corstr='", 
            corstr, "'", ")", sep = "")))$gamma
    }
 
    uncureprob <- as.vector(exp(gamma1 %*% t(Z1))/(1 + exp(gamma1 %*% t(Z1))))
    convergence <- sum((pmt1 - pmt2)^2)
    emes <- list(gamma = gamma1, beta = beta1, bsurv = gSS1, gcor = alpha, gphi = phi, bcor = rho, 
        bphi = pphi, w = gg1, Uncureprob = uncureprob, tau = convergence)
}