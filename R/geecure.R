geecure <- function(formula, cureform, data, id, model = c("para", "semi"), corstr = c("independence", "exchangeable"), Var = TRUE, stdz = FALSE, boots = FALSE, nboot = 100, esmax = 100, eps = 1e-06) {
    call <- match.call()
    model <- match.arg(model)
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
    if (model == "para") {
        X <- model.matrix(attr(mf, "terms"), mf)
    }
    if (model == "semi") {
        X <- model.matrix(attr(mf, "terms"), mf)[, -1]
    }    
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
    esfit <- es(Time, Status, X, Z, id, model, corstr, stdz, esmax, eps)
    gamma <- esfit$gamma
    names(gamma) <- gamma_name
    beta <- esfit$beta
    if (model == "para") {
        kappa <- esfit$kappa
    }
    gcor <- esfit$gcor
    gphi <- esfit$gphi
    bcor <- esfit$bcor
    bphi <- esfit$bphi
    w <- esfit$w
    Lambda <- esfit$Lambda
    if (Var) {
        varfit <- varest(Time, Status, X, Z, id, gamma, beta, kappa, gphi, gcor, bphi, bcor, Lambda, w, model)
        var_gamma <- varfit$varga
        var_beta <- varfit$varbe
        var_kappa <- varfit$varka
        sd_gamma <- varfit$sdga
        sd_beta <- varfit$sdbe
        sd_kappa <- varfit$sdka
    }    
    if (boots) {
        Bootsample <- nboot
        stdz <- stdz
        corstr <- corstr
        model <- model        
        BMc <- matrix(0, Bootsample, ncol(Z))
        BMs <- matrix(0, Bootsample, ncol(X))
        BMnus <- matrix(0, Bootsample, 4)
        if (model == "para") {
            BMnus <- matrix(0, Bootsample, 5)
        }        
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
                 if (model == "para") {
                 X <- model.matrix(attr(mf, "terms"), mf)
                 }
                 if (model == "semi") {
                 X <- model.matrix(attr(mf, "terms"), mf)[, -1]
                 }            
                 Y <- model.extract(mf, "response")
                 if (!inherits(Y, "Surv")) 
                 stop("Response must be a survival object")
                 Time <- Y[, 1]
                 Status <- Y[, 2]  
                 tryboot <- try(es(Time, Status, X, Z, id = id_boot, model, corstr, stdz, esmax = 100, eps = 1e-06), silent = TRUE)  
                 if(is(tryboot,"try-error") == FALSE)
                 break
              }
            esfitboot <- tryboot
            BMc[rrn, ] <- esfitboot$gamma
            BMs[rrn, ] <- esfitboot$beta
            if (model == "semi") {
                BMnus[rrn, ] <- c(esfitboot$gcor, esfitboot$bcor, esfitboot$gphi, esfitboot$bphi)
            }
            if (model == "para") {
                BMnus[rrn, ] <- c(esfitboot$gcor, esfitboot$bcor, esfitboot$gphi, esfitboot$bphi, esfitboot$kappa)
            }
        }
        var_gamma_boots <- apply(BMc, 2, var)
        sd_gamma_boots <- sqrt(var_gamma_boots)
        var_beta_boots <- apply(BMs, 2, var)
        sd_beta_boots <- sqrt(var_beta_boots)
        var_gcor_boots <- var(BMnus[, 1])
        sd_gcor_boots <- sqrt(var_gcor_boots)
        var_bcor_boots <- var(BMnus[, 2])
        sd_bcor_boots <- sqrt(var_bcor_boots)
        if (model == "para") {
            var_kappa_boots <- var(BMnus[, 5])
            sd_kappa_boots <- sqrt(var_kappa_boots)
        }
    }    
    fit <- list()
    class(fit) <- c("geecure")    
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
    if (model == "para") {
        fit$kappa <- kappa
        if (Var) {
            fit$kappa_var <- var_kappa
            fit$kappa_sd <- sd_kappa
            fit$kappa_zvalue <- kappa/sd_kappa
            fit$kappa_pvalue <- (1 - pnorm(abs(fit$kappa_zvalue))) * 2
        }
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
        if (model == "para") {
            fit$boots_kappa_sd <- sd_kappa_boots
            fit$boots_kappa_zvalue <- kappa/sd_kappa_boots
            fit$boots_kappa_pvalue <- (1 - pnorm(abs(fit$boots_kappa_zvalue))) * 2
        }
    }    
    fit$call <- call   
    fit$gamma_name <- gamma_name
    fit$beta_name <- beta_name    
    fit$Time <- Time
    fit$model <- model
    fit$Var <- Var
    fit$boots <- boots
	class(fit) = "geecure"
    fit
#    printgeecure(fit, model, Var, boots)
}

print.geecure <- function(x, ...) {
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
        if (x$model == "semi") {
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
                bt <- array(0, c(length(x$beta) + 1, 4))
                rownames(bt) <- c(x$beta_name, "kappa")
                colnames(bt) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
                bt[, 1] <- c(x$beta, x$kappa)
                bt[, 2] <- c(x$boots_beta_sd, x$boots_kappa_sd)
                bt[, 3] <- c(x$boots_beta_zvalue, x$boots_kappa_zvalue)
                bt[, 4] <- c(x$boots_beta_pvalue, x$boots_kappa_pvalue)
            }
			else {
                bt <- array(0, c(length(x$beta) + 1, 4))
                rownames(bt) <- c(x$beta_name, "kappa")
                colnames(bt) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
                bt[, 1] <- c(x$beta, x$kappa)
                bt[, 2] <- c(x$beta_sd, x$kappa_sd)
                bt[, 3] <- c(x$beta_zvalue, x$kappa_zvalue)
                bt[, 4] <- c(x$beta_pvalue, x$kappa_pvalue)
            }
        }
    }
	else {
        if (x$model == "semi") {
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
		else {
            if (x$boots) {
                bt <- array(0, c(length(x$beta) + 1, 4))
                rownames(bt) <- c(x$beta_name, "kappa")
                colnames(bt) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
                bt[, 1] <- c(x$beta, x$kappa)
                bt[, 2] <- c(x$boots_beta_sd, x$boots_kappa_sd)
                bt[, 3] <- c(x$boots_beta_zvalue, x$boots_kappa_zvalue)
                bt[, 4] <- c(x$boots_beta_pvalue, x$boots_kappa_pvalue)
            }
			else {
                bt <- array(0, c(length(x$beta) + 1, 1))
                bt[, 1] <- c(x$beta, x$kappa)
                rownames(bt) <- c(x$beta_name, "kappa")
                colnames(bt) <- "Estimate"
            }
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
