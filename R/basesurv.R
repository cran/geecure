basesurv <- function(Time, Status, X, beta, Lambda, w, model) {
    if (model == "para") {
        kappa <- -sum(Status)/sum(w * log(Time) * (Status - Lambda * exp(beta %*% t(X))))
        bcumhaz <- Time^kappa
        uncuresurv <- exp(-bcumhaz * exp(beta %*% t(X)))
    }
    if (model == "semi") {
        t2 <- Time
        t11 <- sort(Time)
        c11 <- Status[order(Time)]
        x111 <- X[order(Time), ]
        g11 <- w[order(Time)]
        tt1 <- unique(t11[c11 == 1])
        kk <- length(table(t11[c11 == 1]))
        dd <- as.matrix(table(t11[c11 == 1]))
        gSS <- rep(0, kk)
        gSS1 <- rep(1, kk)
        Kn <- length(Time)
        gSS[1] <- dd[1]/(sum(g11[min((1:Kn)[t11 == tt1[1]]):Kn] * exp(beta %*% t(x111[min((1:Kn)[t11 == tt1[1]]):Kn, ]))))
        for (i in 1:(kk - 1)) {
            gSS[i + 1] <- gSS[i] + dd[i + 1]/(sum(g11[min((1:Kn)[t11 == tt1[i + 1]]):Kn] * exp(beta %*% t(x111[min((1:Kn)[t11 == tt1[i + 1]]):Kn, ]))))
        }
        gSS1 <- exp(-gSS)        
        gSS2 <- rep(0, Kn)
        gSS3 <- rep(0, Kn)
        for (i in 1:(Kn)) {
            kk1 <- 1
            if (t2[i] < tt1[1]) {
                gSS2[i] <- 1
                gSS3[i] <- 1e-08
            } else {
                if (t2[i] >= tt1[kk]) {
                  gSS2[i] <- 0
                  gSS3[i] <- gSS[kk]
                } else {
                  repeat {
                    if (t2[i] >= tt1[kk1]) 
                      kk1 <- kk1 + 1 else break
                  }
                  { gSS2[i] <- (gSS1[kk1 - 1])^(exp(beta %*% X[i, ]))
                    gSS3[i] <- gSS[kk1 - 1]
                  }
                }
            }
        }
        bcumhaz <- gSS3
        uncuresurv <- gSS2
        kappa <- 1
    }
    list(uncuresurv = uncuresurv, bcumhaz = bcumhaz, kappa = kappa)
}