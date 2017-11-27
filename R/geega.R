geega <- function(w, Z, gamma, id, corstr) {
    gamma1 <- rep(0, length(gamma))
    id <- id
    K <- length(unique(id))
    n <- as.vector(table(id))
    alpha <- 0
    phi <- 1
    P <- exp(gamma %*% t(Z))/(1 + exp(gamma %*% t(Z)))
    PP <- P * (1 - P)
    repeat {
        repeat {
            D <- diag(c(PP[id == 1], 0))[1:n[1], 1:n[1]] %*% diag(rep(1, n[1])) %*% Z[id == 1, ]
            for (i in 2:K) {
                D <- rbind(D, diag(c(PP[id == i], 0))[1:n[i], 1:n[i]] %*% diag(rep(1, n[i])) %*% Z[id == i, ])
            }
            S <- w - P
            R <- matrix(alpha, n[1], n[1])
            diag(R) <- 1
            V <- sqrt(diag(c(PP[id == 1], 0))[1:n[1], 1:n[1]]) %*% R %*% sqrt(diag(c(PP[id == 1], 0))[1:n[1], 1:n[1]]) * phi
            for (i in 2:K) {
                R <- matrix(alpha, n[i], n[i])
                diag(R) <- 1
                V <- bdiag(V, sqrt(diag(c(PP[id == i], 0))[1:n[i], 1:n[i]]) %*% R %*% sqrt(diag(c(PP[id == i], 0))[1:n[i], 1:n[i]]) * phi)
            }
            V <- as.matrix(V)
            ZU <- gamma %*% t(D) + S
            newgamma <- t(ginv(t(D) %*% ginv(V) %*% D) %*% t(D) %*% ginv(V) %*% t(ZU))
            if (any(abs(newgamma - gamma) > 1e-06)) {
                gamma <- newgamma
                P <- exp(gamma %*% t(Z))/(1 + exp(gamma %*% t(Z)))
                PP <- P * (1 - P)
            } else break
        }
        if (any(abs(gamma - gamma1) > 1e-06)) {
            gamma1 <- gamma
            if (corstr == "exchangeable") {
                P <- exp(gamma1 %*% t(Z))/(1 + exp(gamma1 %*% t(Z)))
                PP <- P * (1 - P)
                rr <- as.vector((w - P)/sqrt(PP))
                rr1 <- 0
                phi <- sum(rr^2)/(sum(n) - dim(Z)[2])
                rrm <- matrix(0, ncol = K, nrow = max(n))
                for (i in 1:K) {
                  rrm[(1:n[i]), i] <- rr[id == i]
                }
                rr <- rrm
                rr <- t(rr)
                for (i in 1:K) {
                  if (n[i] == 1) {
                    rr1 <- rr1 + rr[i, 1]
                  } else {
                    for (j in 1:(n[i] - 1)) rr1 <- rr1 + rr[i, j] * sum(rr[i, (j + 1):n[i]])
                  }
                }
                alpha <- (phi^(-1)) * rr1/(sum(n * (n - 1))/2 - dim(Z)[2])
            }
            if (corstr == "independence") {
                phi <- 1
                alpha <- 0
            }
        } else break
    }
    list(gamma = gamma1, alpha = alpha, phi = phi)
}