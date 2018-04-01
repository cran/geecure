geecure2 <- function(formula, cureform, data, id, corstr = c("independence", "exchangeable"), Var = TRUE, stdz = FALSE, boots = FALSE, nboot = 100, esmax = 100, eps = 1e-06) {
    call <- match.call()
    data <- data
    id <- id
    uid <- sort(unique(id))
    newid <- rep(0, length(id))
    for (i in 1:length(id)) {
        j <- 1
        repeat {
            if (id[i] != uid[j]) 
                j <- j + 1 else {
                newid[i] <- j
                break
            }
        }
    }
    data$id <- newid    
    data1 <- data[data$id == 1, ]
    for (i in 2:length(uid)) {
        data1 <- rbind(data1, data[data$id == i, ])
    }    
    data <- data1
    id <- data$id    
    Kn <- length(id)
    K <- length(unique(id))
    n <- as.vector(table(id))
    mf <- model.frame(formula, data)
    mf1 <- model.frame(cureform, data)
    Z <- model.matrix(attr(mf1, "terms"), mf1)  
    X <- model.matrix(attr(mf, "terms"), mf)[, -1]
    X <- as.matrix(X)
    colnames(X) <- colnames(model.matrix(attr(mf, "terms"), mf))[-1]         
    gamma_name <- colnames(Z)  
    gamma_length <- ncol(Z)  
    beta_name <- colnames(X) 
    beta_length <- ncol(X)
    Y <- model.extract(mf, "response")
    if (!inherits(Y, "Surv")) 
        stop("Response must be a survival object")
    Time <- Y[, 1]
    Status <- Y[, 2]    
    stdz <- stdz
    esmax <- esmax
    eps <- eps
    esfit <- emes(Time, Status, X, Z, id, corstr, stdz, esmax, eps)
    gamma <- esfit$gamma
    names(gamma) <- gamma_name
    beta <- esfit$beta
    gcor <- esfit$gcor
    gphi <- esfit$gphi
    bcor <- esfit$bcor
    bphi <- esfit$bphi
    w <- esfit$w
    bsurv <- esfit$bsurv
    if (Var) {
        varfit <- varest2(Time, Status, X, Z, id, gamma, beta, bsurv, w)
        var_gamma <- varfit$varga
        var_beta <- varfit$varbe
        sd_gamma <- varfit$sdga
        sd_beta <- varfit$sdbe
    }    
    if (boots) {
        Bootsample <- nboot
        stdz <- stdz
        corstr <- corstr       
        BMc <- matrix(0, Bootsample, ncol(Z))
        BMs <- matrix(0, Bootsample, ncol(X))
        BMnus <- matrix(0, Bootsample, 4)       
        for (rrn in 1:Bootsample) {
        repeat{  bootid<- sample((1:K), replace = TRUE)
                 bootdata <- data[id == bootid[1], ]
                 bootdata$id <- rep(1, sum(id == bootid[1]))
                 for (ll in 2:K) {
                 bootdata1 <- data[id == bootid[ll], ]
                 bootdata1$id <- rep(ll, sum(id == bootid[ll]))
                 bootdata <- rbind(bootdata, bootdata1)
                 }            
                 id_boot <- bootdata$id
                 Kn_boot <- length(id_boot)
                 K_boot <- length(unique(id_boot))
                 n_boot <- as.vector(table(id_boot))
                 mf <- model.frame(formula, bootdata)
                 mf1 <- model.frame(cureform, bootdata)
                 Z <- model.matrix(attr(mf1, "terms"), mf1)
                 X <- model.matrix(attr(mf, "terms"), mf)[, -1]
                 X <- as.matrix(X)
                 colnames(X) <- colnames(model.matrix(attr(mf, "terms"), mf))[-1]
                 Y <- model.extract(mf, "response")
                 if (!inherits(Y, "Surv")) 
                 stop("Response must be a survival object")
                 Time <- Y[, 1]
                 Status <- Y[, 2]  
                 tryboot <- try(emes(Time, Status, X, Z, id = id_boot, corstr, stdz, esmax = 100, eps = 1e-06), silent = TRUE)  
                 if(is(tryboot,"try-error") == FALSE)
                 break
              }
            esfitboot <- tryboot
            BMc[rrn, ] <- esfitboot$gamma
            BMs[rrn, ] <- esfitboot$beta
            BMnus[rrn, ] <- c(esfitboot$gcor, esfitboot$bcor, esfitboot$gphi, esfitboot$bphi)
        }
        var_gamma_boots <- apply(BMc, 2, var)
        sd_gamma_boots <- sqrt(var_gamma_boots)
        var_beta_boots <- apply(BMs, 2, var)
        sd_beta_boots <- sqrt(var_beta_boots)
        var_gcor_boots <- var(BMnus[, 1])
        sd_gcor_boots <- sqrt(var_gcor_boots)
        var_bcor_boots <- var(BMnus[, 2])
        sd_bcor_boots <- sqrt(var_bcor_boots)
    }    
    fit <- list()
    class(fit) <- c("geecure2")    
    fit$gamma <- gamma
    if (Var) {
        fit$gamma_var <- var_gamma
        fit$gamma_sd <- sd_gamma
        fit$gamma_zvalue <- gamma/sd_gamma
        fit$gamma_pvalue <- (1 - pnorm(abs(fit$gamma_zvalue))) * 2
    }
    fit$beta <- beta
    if (Var) {
        fit$beta_var <- var_beta
        fit$beta_sd <- sd_beta
        fit$beta_zvalue <- beta/sd_beta
        fit$beta_pvalue <- (1 - pnorm(abs(fit$beta_zvalue))) * 2
    }    
    fit$alpha <- gcor
    fit$rho <- bcor
    fit$bphi <- bphi
    fit$gphi <- gphi    
    fit$num_of_clusters <- K
    fit$max_cluster_size <- max(n)    
    if (boots) {
        fit$boots_gamma_sd <- sd_gamma_boots
        fit$boots_beta_sd <- sd_beta_boots
        fit$boots_gcor_sd <- sd_gcor_boots
        fit$boots_bcor_sd <- sd_bcor_boots
        fit$boots_gamma_zvalue <- gamma/sd_gamma_boots
        fit$boots_beta_zvalue <- beta/sd_beta_boots
        fit$boots_gcor_zvalue <- gcor/sd_gcor_boots
        fit$boots_bcor_zvalue <- bcor/sd_bcor_boots
        fit$boots_gamma_pvalue <- (1 - pnorm(abs(fit$boots_gamma_zvalue))) * 2
        fit$boots_beta_pvalue <- (1 - pnorm(abs(fit$boots_beta_zvalue))) * 2
        fit$boots_gcor_pvalue <- (1 - pnorm(abs(fit$boots_gcor_zvalue))) * 2
        fit$boots_bcor_pvalue <- (1 - pnorm(abs(fit$boots_bcor_zvalue))) * 2
    }    
    fit$call <- call   
    fit$gamma_name <- gamma_name
    fit$beta_name <- beta_name    
    fit$Time <- Time
    fit$Var <- Var
    fit$boots <- boots
    class(fit) = "geecure2"
    fit
}

print.geecure2 <- function(x, ...) {
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
    }
    cat("\nCure Probability Model:\n")
    if (x$Var) {
        if (x$boots) {
            gm <- array(x$gamma, c(length(x$gamma), 4))
            rownames(gm) <- x$gamma_name
            colnames(gm) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
            gm[, 2] <- x$boots_gamma_sd
            gm[, 3] <- x$boots_gamma_zvalue
            gm[, 4] <- x$boots_gamma_pvalue
        }
		else {
            gm <- array(x$gamma, c(length(x$gamma), 4))
            rownames(gm) <- x$gamma_name
            colnames(gm) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
            gm[, 2] <- x$gamma_sd
            gm[, 3] <- x$gamma_zvalue
            gm[, 4] <- x$gamma_pvalue
        }
    }
	else {
        if (x$boots) {
            gm <- array(x$gamma, c(length(x$gamma), 4))
            rownames(gm) <- x$gamma_name
            colnames(gm) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
            gm[, 2] <- x$boots_gamma_sd
            gm[, 3] <- x$boots_gamma_zvalue
            gm[, 4] <- x$boots_gamma_pvalue
        }
		else {
            gm <- array(x$gamma, c(length(x$gamma), 1))
            rownames(gm) <- x$gamma_name
            colnames(gm) <- "Estimate"
        }
    }
    print(gm)
    cat("\n")
    cat("\nFailure Time Distribution Model:\n")
    if (x$Var) {
            if (x$boots) {
                bt <- array(x$beta, c(length(x$beta), 4))
                rownames(bt) <- x$beta_name
                colnames(bt) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
                bt[, 2] <- x$boots_beta_sd
                bt[, 3] <- x$boots_beta_zvalue
                bt[, 4] <- x$boots_beta_pvalue
            }
			else {
                bt <- array(x$beta, c(length(x$beta), 4))
                rownames(bt) <- x$beta_name
                colnames(bt) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
                bt[, 2] <- x$beta_sd
                bt[, 3] <- x$beta_zvalue
                bt[, 4] <- x$beta_pvalue
            }
    }
	else {
            if (x$boots) {
                bt <- array(x$beta, c(length(x$beta), 4))
                rownames(bt) <- x$beta_name
                colnames(bt) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
                bt[, 2] <- x$boots_beta_sd
                bt[, 3] <- x$boots_beta_zvalue
                bt[, 4] <- x$boots_beta_pvalue
            }
			else {
                bt <- array(x$beta, c(length(x$beta), 1))
                rownames(bt) <- x$beta_name
                colnames(bt) <- "Estimate"
            }        
    }
    print(bt)
    cat("\n")
    cat("\nEstimated Correlation Parameters:\n")
    if (x$boots) {
        dep <- array(0, c(2, 4))
        rownames(dep) <- c("rho_1", "rho_2")
        colnames(dep) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
        dep[, 1] <- c(x$alpha, x$rho)
        dep[, 2] <- c(x$boots_gcor_sd, x$boots_bcor_sd)
        dep[, 3] <- c(x$boots_gcor_zvalue, x$boots_bcor_zvalue)
        dep[, 4] <- c(x$boots_gcor_pvalue, x$boots_bcor_pvalue)
    }
	else {
        dep <- array(0, c(2, 1))
        rownames(dep) <- c("rho_1", "rho_2")
        colnames(dep) <- "Estimate"
        dep[, 1] <- c(x$alpha, x$rho)
    }
    print(dep)
    cat("\n")    
    cat("Number of clusters:", x$num_of_clusters)
    cat("       Maximum cluster size:", x$max_cluster_size)
    cat("\n")
    invisible(x)
}
