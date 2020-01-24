library(BASSLINE)
LST <- MCMC_LST(1000,20,10,cancer[,1], cancer[,2], cancer[,3:11])
LML <- LML_LST(thin = 20, Time = cancer[,1], Cens = cancer[,2],
               X = cancer[,3:11], chain = LST)
LST.DIC <- DIC_LST(Time = cancer[,1], Cens = cancer[,2], X = cancer[,3:11],
                   chain = LST)
LST.CD <- CaseDeletion_LST(Time = cancer[,1], Cens = cancer[,2],
                           cancer[,3:11], chain = LST)
