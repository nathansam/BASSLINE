library(BASSLINE)
LEP <- MCMC_LEP(100,20,1,cancer[,1], cancer[,2], cancer[,3:11])
LML <- LML_LEP(thin = 2, Time = cancer[,1], Cens = cancer[,2],
               X = cancer[,3:11], chain = LEP)
LEP.DIC <- DIC_LEP(Time = cancer[,1], Cens = cancer[,2], X = cancer[,3:11],
                   chain = LEP)

