initial_Lambda <- function(Time, Status, X, Z, id, model, corstr) {
    w <- Status
    t2 <- Time
    K <- length(unique(id))
    n <- as.vector(table(id))
    Kn <- sum(n)
    cens <- Status
    t11 <- sort(Time)
    c11 <- Status[order(Time)]
    x111 <- X[order(Time), ]
    g11 <- w[order(Time)]
    tt1 <- unique(t11[c11 == 1])
    kk <- length(table(t11[c11 == 1]))
    dd <- as.matrix(table(t11[c11 == 1]))
    gg1 <- Status
    gg2 <- log(Time)
    gg1[gg1 < 1e-06] <- 1e-06
    gg3 <- log(gg1) + log(Time)
    pmt1c <- eval(parse(text = paste("geese", "(", "w~Z[,-1]", ",id=", "id", ",family = binomial", ",corstr='", corstr, "'", ")", sep = "")))$beta
    pmt1s <- eval(parse(text = paste("geese", "(", "w~X-1+", "offset(", "gg3", ")", ",id=", "id", ",family = poisson", ",corstr='", corstr, "'", ")", sep = "")))$beta
    ppmt2 <- c(pmt1c, pmt1s)
    KK <- 1
    repeat {
        gSSS1 <- rep(0, kk)
        KK1 <- 1
        repeat {
            gSS <- rep(0, kk)
            gSS1 <- rep(1, kk)            
            gSS[1] <- dd[1]/(sum(g11[min((1:(Kn))[t11 == tt1[1]]):(Kn)] * exp(x111[min((1:(Kn))[t11 == tt1[1]]):(Kn), ] %*% pmt1s)))
            for (i in 1:(kk - 1)) {
                gSS[i + 1] <- gSS[i] + dd[i + 1]/(sum(g11[min((1:(Kn))[t11 == tt1[i + 1]]):(Kn)] * exp(x111[min((1:(Kn))[t11 == tt1[i + 1]]):(Kn), ] %*% pmt1s)))
            }
            gSS1 <- exp(-gSS)
            gSS2 <- rep(0, Kn)
            gSS3 <- rep(0, Kn)
            for (i in 1:Kn) {
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
                    { gSS2[i] <- (gSS1[kk1 - 1])^(exp(X[i, ] %*% pmt1s))
                      gSS3[i] <- gSS[kk1 - 1]
                    }
                  }
                }
            }
            gg2 <- log(gSS3)
            gg3 <- log(gg1) + gg2
            ww2 <- eval(parse(text = paste("geese", "(", "w~X-1+", "offset(", "gg3", ")", ",id=", "id", ",family = poisson", ",corstr='", corstr, "'", ")", sep = "")))
            if (KK1 < 100 && (any(abs(ww2$beta - pmt1s) > 1e-06) || any(abs(gSS1 - gSSS1) > 1e-06))) {
                pmt1s <- c(ww2$beta)
                gSSS1 <- gSS1
                KK1 <- KK1 + 1
            } else {
                gg1 <- Status + ((1 - Status) * exp(Z %*% pmt1c) * gSS2)/(1 + exp(Z %*% pmt1c) * gSS2)
                g11 <- gg1[order(t2)]
                gg1[gg1 < 1e-06] <- 1e-06
                gg3 <- log(gg1) + gg2
                break
            }
        }
        pmt2c <- eval(parse(text = paste("geese", "(", "Status~Z[,-1]", ",id=", "id", ",family = binomial", ",corstr='", corstr, "'", ")", sep = "")))$beta
        pmt2s <- eval(parse(text = paste("geese", "(", "Status~X-1+", "offset(", "gg3", ")", ",id=", "id", ",family = poisson", ",corstr='", corstr, "'", ")", sep = "")))$beta
        if (any(abs(pmt2c - pmt1c) > 1e-06) || max((pmt2s - pmt1s)^2) > 1e-08) {
            pmt1c <- pmt2c
            pmt1s <- pmt2s
            KK <- KK + 1
        } else break
    }
    Lambda <- gSS3
    initial_Lambda <- list(Lambda = Lambda)
}