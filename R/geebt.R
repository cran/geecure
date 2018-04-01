geebt <- function(Status, Lambda, X, beta, w, id, corstr) {
    K <- length(unique(id))
    n <- as.vector(table(id))
    gbeta <- rep(0, length(beta))
    newY1 <- Status/Lambda
    W1 <- diag(w * Lambda)
    mu <- exp(beta %*% t(X))
    SK1 <- 1
    repeat {
        if (corstr == "independence") {
            pphi <- 1
            rho <- 0
        }
        if (corstr == "exchangeable") {
            res <- as.vector((newY1 - mu)/sqrt(mu))
            rres <- 0
            pphi <- sum(res^2)/(sum(n) - dim(X)[2])
            resm <- matrix(0, ncol = K, nrow = max(n))
            for (i in 1:K) {
                resm[1:n[i], i] <- res[id == i]
            }
            res <- resm
            res <- t(res)
            for (i in 1:K) {
                if (n[i] == 1) 
                  rres <- rres + res[i, 1] else {
                  for (j in 1:(n[i] - 1)) rres <- rres + res[i, j] * sum(res[i, (j + 1):n[i]])
                }
            }
            rho <- (pphi^(-1)) * rres/(sum(n * (n - 1))/2 - dim(X[, -1, drop = FALSE])[2])
        }
        SK <- 1
        repeat {
            D1 <- diag(c(mu[id == 1], 0))[1:n[1], 1:n[1]] %*% diag(rep(1, n[1])) %*% (X[id == 1, ])
            for (i in 2:K) {
                D1 <- rbind(D1, diag(c(mu[id == i], 0))[1:n[i], 1:n[i]] %*% diag(rep(1, n[i])) %*% (X[id == i, ]))
            }
            S1 <- newY1 - mu
            R1 <- matrix(rho, n[1], n[1])
            diag(R1) <- 1
            V1 <- sqrt(diag(c(mu[id == 1], 0))[1:n[1], 1:n[1]]) %*% R1 %*% sqrt(diag(c(mu[id == 1], 0))[1:n[1], 1:n[1]]) * pphi
            for (i in 2:K) {
                R1 <- matrix(rho, n[i], n[i])
                diag(R1) <- 1
                V1 <- bdiag(V1, sqrt(diag(c(mu[id == i], 0))[1:n[i], 1:n[i]]) %*% R1 %*% sqrt(diag(c(mu[id == i], 0))[1:n[i], 1:n[i]]) * pphi)
            }
            V1 <- as.matrix(V1)
            ZU1 <- beta %*% t(D1) + S1
            newbeta <- t(ginv(t(D1) %*% ginv(V1) %*% W1 %*% D1) %*% t(D1) %*% ginv(V1) %*% W1 %*% t(ZU1))
            if (any(abs(newbeta - beta) > 1e-06) && (SK <= 500)) {
                beta <- newbeta
                mu <- exp(beta %*% t(X))
                SK <- SK + 1
            } else break
        }
        if (any(abs(beta - gbeta) > 1e-06) && (SK1 < 30)) {
            gbeta <- beta
            mu <- exp(gbeta %*% t(X))
            SK1 <- SK1 + 1
        } else break
    }
    list(beta = gbeta, rho = rho, pphi = pphi)
}