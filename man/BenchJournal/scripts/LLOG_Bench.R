library(BASSLINE)
LLOG <- MCMC_LLOG(N = 1000, thin = 20, burn = 40, Time = cancer[,1],
                  Cens = cancer[,2], X = cancer[,3:11])

LML <- LML_LLOG(thin = 20, Time = cancer[,1], Cens = cancer[,2],
                X = cancer[,3:11], chain = LLOG)

LLOG.DIC <- DIC_LLOG(Time = cancer[,1], Cens = cancer[,2], X = cancer[,3:11],
                     chain = LLOG)
