varest2 <- function(Time, Status, X, Z, id, gamma, beta, bsurv, w) {
        K <- length(unique(id))
        n <- as.vector(table(id))
        Kn <- length(id)
        zzzz <- Z
        xxxx <- X
        gg1 <- w
        ga <- gamma
        be <- beta
                 
        t2 <- Time
        c1 <- Status
        t11 <- sort(Time)
        c11 <- Status[order(Time)]
        tt1 <- unique(t11[c11 == 1])
        kk <- length(table(t11[c11 == 1]))

        gSS <- -log(bsurv)
        gS <- c(gSS[1], gSS[2:kk] - gSS[1:(kk - 1)])   

        ss <- 0
    
        for (k in 1:K) {

        xxx1 <- xxxx[id == k, , drop = FALSE]
        zzz1 <- zzzz[id == k, , drop = FALSE]
        gg11 <- gg1[id == k]
        c111 <- c1[id == k]
        t21 <- t2[id == k]        
        
        a11 <- rep(0, dim(zzzz)[2])
        for (s in 1:dim(zzzz)[2]) {
            a11[s] <- sum((gg11 - (1 - (1 + exp(zzz1 %*% ga))^(-1))) * zzz1[, s])
        }
       
        b11=rep(0,dim(xxxx)[2])
      
        for(s1 in 1:dim(xxxx)[2])
       {  
        for (s in 1:kk) {
            
            AA <- sum(gS[s] * exp(xxx1[(c111 == 1) & (t21 == tt1[s]), , drop = FALSE] %*% be) * xxx1[(c111 == 1) & (t21 == tt1[s]), s1, drop = FALSE]/(1 - exp(-gS[s] * exp(xxx1[(c111 == 1) & (t21 == tt1[s]), , drop = FALSE] %*% be)))) - gS[s] * 
                sum(gg11[t21 >= tt1[s]] * exp(xxx1[t21 >= tt1[s], , drop = FALSE] %*% be) * xxx1[t21 >= tt1[s], s1, drop = FALSE])
            b11[s1] <- b11[s1] + AA
        }
       }
                
        d11 <- rep(0, kk)
        
        for (s in 1:kk) {
            d11[s] <- sum(exp(xxx1[(c111 == 1) & (t21 == tt1[s]), , drop = FALSE] %*% be)/(1 - exp(-gS[s] * exp(xxx1[(c111 == 1) & (t21 == tt1[s]), , drop = FALSE] %*% be)))) - sum(gg11[t21 >= tt1[s]] * exp(xxx1[t21 >= tt1[s], , drop = FALSE] %*% be))
        }
        
        ss <- ss + t(t(c(a11, b11, d11))) %*% t(c(a11, b11, d11))        
    } 

        AAA <- 0
    
    for (k in 1:Kn) {
        Aa <- ((1 - (1 + exp(zzzz[k, , drop = FALSE] %*% ga))^(-1)) * ((1 + exp(zzzz[k, , drop = FALSE] %*% ga))^(-1)))[1, 1]
        AAA <- AAA + Aa * (t(zzzz[k, , drop = FALSE]) %*% zzzz[k, , drop = FALSE])
    }
    
     BBB=matrix(0,dim(xxxx)[2],dim(xxxx)[2])
    
    for (k in 1:dim(xxxx)[2])
    {  for(s1 in 1:dim(xxxx)[2])
      { for (s in 1:(kk)) {
           Bm <- gS[s] * exp(xxxx[(c1 == 1) & (t2 == tt1[s]), , drop = FALSE] %*% be)
           BBB[k, s1] <- BBB[k, s1] + (sum((Bm * (exp(-Bm) + Bm * exp(-Bm) - 1)/(1 - exp(-Bm))^2) * (xxxx[(c1 == 1) & (t2 == tt1[s]), k, drop = FALSE] * xxxx[(c1 == 1) & (t2 == tt1[s]), s1, drop = FALSE])) + gS[s] * sum(gg1[t2 >= tt1[s]] * 
           exp(xxxx[t2 >= tt1[s], , drop = FALSE] %*% be) * (xxxx[t2 >= tt1[s], k, drop = FALSE] * xxxx[t2 >= tt1[s], s1, drop = FALSE])))
        } 
      }
    }

    BBC <- matrix(0, dim(xxxx)[2], kk)
  
    for(s1 in 1:dim(xxxx)[2])
    {  
      for (s in 1:kk) {
        BCm <- gS[s] * exp(xxxx[(c1 == 1) & (t2 == tt1[s]), , drop = FALSE] %*% be)
        BBC[s1, s] <- sum(exp(xxxx[(c1 == 1) & (t2 == tt1[s]), , drop = FALSE] %*% be) * (exp(-BCm) + BCm * exp(-BCm) - 1)/((1 - exp(-BCm))^2) * xxxx[(c1 == 1) & (t2 == tt1[s]), s1, drop = FALSE]) + sum(gg1[t2 >= tt1[s]] * 
        exp(xxxx[t2 >= tt1[s], , drop = FALSE] %*% be) * xxxx[t2 >= tt1[s], s1, drop = FALSE])
      }
    }
    
    CCC <- rep(0, kk)
    
    for (s in 1:kk) {
        CCm <- gS[s] * exp(xxxx[(c1 == 1) & (t2 == tt1[s]), , drop = FALSE] %*% be)
        CCC[s] <- sum(exp(2 * (xxxx[(c1 == 1) & (t2 == tt1[s]), , drop = FALSE] %*% be) - CCm)/(1 - exp(-CCm))^2)
    }

    AAA1 <- 0
    
    for (k in 1:Kn) {
        AAA1 <- AAA1 + (gg1[k] * (1 - gg1[k]) * (t(zzzz[k, , drop = FALSE]) %*% zzzz[k, , drop = FALSE]))        
    }
    
    tt1 <- c(tt1, Inf)
    
    BBB1 <- matrix(0,dim(xxxx)[2],dim(xxxx)[2])

    for(k in 1:dim(xxxx)[2])          
     {
       for(s1 in 1:dim(xxxx)[2])               
         {
            for (s in 1:kk) {
            BBB1[k, s1] <- BBB1[k, s1] + sum(gg1[(t2 == tt1[s]) | (t2 >= tt1[s] & t2 < tt1[s + 1] & c1 == 0)] * (1 - gg1[(t2 == tt1[s]) | (t2 >= tt1[s] & t2 < tt1[s + 1] & c1 == 0)]) * 
            (exp(xxxx[(t2 == tt1[s]) | (t2 >= tt1[s] & t2 < tt1[s + 1] & c1 == 0), , drop = FALSE] %*% be) * sum(gS[1:s]))^2 * (xxxx[(t2 == tt1[s]) | (t2 >= tt1[s] & t2 < tt1[s + 1] & c1 == 0), k, drop = FALSE] * 
            xxxx[(t2 == tt1[s]) | (t2 >= tt1[s] & t2 < tt1[s + 1] & c1 == 0), s1, drop = FALSE]))
           }
         }
     }
    
    BBC1 <- matrix(0,dim(xxxx)[2], kk)
    
    for(s1 in 1:dim(xxxx)[2])  
      {
        for (s in 1:kk) {
          BBC11 <- 0
          for (j in s:kk) {
            BBC11 <- BBC11 + sum(gg1[(t2 == tt1[j]) | (t2 >= tt1[j] & t2 < tt1[j + 1] & c1 == 0)] * (1 - gg1[(t2 == tt1[j]) | (t2 >= tt1[j] & t2 < tt1[j + 1] & c1 == 0)]) * 
            (exp(xxxx[(t2 == tt1[j]) | (t2 >= tt1[j] & t2 < tt1[j + 1] & c1 == 0), , drop = FALSE] %*% be))^2 * xxxx[(t2 == tt1[j]) | (t2 >= tt1[j] & t2 < tt1[j + 1] & c1 == 0), s1, drop = FALSE] * sum(gS[1:j]))
          }
          BBC1[s1, s] <- BBC11
        }
      }

    BBA1 <- matrix(0,dim(xxxx)[2],dim(zzzz)[2])
    
    for(s1 in 1:dim(xxxx)[2])  
     {
       for (j in 1:dim(zzzz)[2]) {
          for (s in 1:kk) {
            BBA1[s1, j] <- BBA1[j] + sum(gg1[(t2 == tt1[s]) | (t2 >= tt1[s] & t2 < tt1[s + 1] & c1 == 0)] * (1 - gg1[(t2 == tt1[s]) | (t2 >= tt1[s] & t2 < tt1[s + 1] & c1 == 0)]) * 
            xxxx[(t2 == tt1[s]) | (t2 >= tt1[s] & t2 < tt1[s + 1] & c1 == 0), s1, drop = FALSE] * zzzz[(t2 == tt1[s]) | (t2 >= tt1[s] & t2 < tt1[s + 1] & c1 == 0), j, drop = FALSE] * 
            exp(xxxx[(t2 == tt1[s]) | (t2 >= tt1[s] & t2 < tt1[s + 1] & c1 == 0), , drop = FALSE] %*% be) * sum(gS[1:s]))
            }
        }
     }
    
    CCC1 <- matrix(0, kk, kk)
    
    for (i in 1:kk) {
        for (j in 1:kk) {
            CCC1[i, j] <- sum(gg1[t2 >= tt1[max(i, j)]] * (1 - gg1[t2 >= tt1[max(i, j)]]) * exp(2 * (xxxx[t2 >= tt1[max(i, j)], , drop = FALSE] %*% be)))
        }
    }
    
    CCA1 <- matrix(0, kk, dim(zzzz)[2])
    
    for(s in 1:kk)
    {
       for (j in 1:dim(zzzz)[2]) 
       {
           CCA1[s, j] <- sum(gg1[t2 >= tt1[s]] * (1 - gg1[t2 >= tt1[s]]) * exp(xxxx[t2 >= tt1[s], , drop = FALSE] %*% be) * zzzz[t2 >= tt1[s], j, drop = FALSE])
       }
    }
 
    M11 <- AAA - AAA1
    M12 <- BBA1
    M13 <- CCA1
    
    M21 <- t(M12)
    M22 <- BBB - BBB1
    M23 <- t(BBC - BBC1)
    
    M31 <- t(M13)
    M32 <- t(M23)
    M33 <- diag(CCC) - CCC1
    
    M <- cbind(rbind(M11, M12, M13), rbind(M21, M22, M23), rbind(M31, M32, M33))
    M1 <- solve(M) %*% ss %*% solve(M)
        
    vcm <- M1
    var_gamma <- diag(vcm)[1:dim(zzzz)[2]]
    var_beta <- diag(vcm)[(dim(zzzz)[2] + 1):(dim(zzzz)[2] + dim(xxxx)[2])]
    
    sd_gamma <- sqrt(var_gamma)
    sd_beta <- sqrt(var_beta)
    
    list(varga = var_gamma, varbe = var_beta, sdga = sd_gamma, sdbe = sd_beta)
}