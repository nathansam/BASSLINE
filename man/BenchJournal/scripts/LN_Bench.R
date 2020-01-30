library(BASSLINE)
LN <- MCMC_LN(1000, 20, 10, cancer[,1], cancer[,2], cancer[,3:11])
LML <- LML_LN(thin = 20, Time = cancer[,1], Cens = cancer[,2],
               X = cancer[,3:11], chain = LN)
LN.DIC <- DIC_LN(Time = cancer[,1], Cens = cancer[,2], X = cancer[,3:11],
                   chain = LN)
LN.CD <- CaseDeletion_LN(Time = cancer[,1], Cens = cancer[,2],
                           cancer[,3:11], chain = LN)
