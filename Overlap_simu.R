#estimate A -> B effect; 
library(mvtnorm);
library(getopt);

trait <- 0;
#m1 <- matrix(NA, nrow=5, ncol=11); 
#m2 <- matrix(NA, nrow=5, ncol=11);
mylist <- data.frame(m1=c(), m2=c());
rho <- 0.5
sigma <- matrix(c(1, rho, rho, 1), ncol=2)

for (tt in seq(20, 100, by=20)){
  trait <- trait + 1
  count <- 0;
  for (j in seq(0, 0.1, by=0.01)){
    count <- count + 1
    p1 <- 0; p2 <- 0
    
    for (i in 1:5) {
      beta.AB <- j;
      N.size <- 1000;
      beta.GA <- rnorm(tt, 0, 0.2);
      # beta.GB <- beta.GA * beta.AB;
      N <- c(N.size, N.size, N.size, N.size);
      G <- matrix(nrow=sum(N),ncol=tt)
      YA <- 0;YB <- 0;
      for(ii in 1:4) {
        G[((ii-1)*N.size + 1):(ii*N.size),] <- matrix(rbinom(N.size*tt, 2, 0.25), ncol=tt);
        E <- rmvnorm(n=N.size, mean=c(1,2), sigma=sigma);
        YA[((ii-1)*N.size + 1):(ii*N.size)] <- G[((ii-1)*N.size+1):(ii*N.size),]%*%beta.GA + E[,1];
        YB[((ii-1)*N.size + 1):(ii*N.size)] <- beta.AB*YA[((ii-1)*N.size + 1):(ii*N.size)] + E[,2];         
      }
      
      G1 <- G[1:(2*N.size), ];
      Y1 <- YA[1:(2*N.size)];
      
      G2 <- G[(N.size + 1):(3*N.size), ];
      Y2 <- YB[(N.size + 1):(3*N.size)];
      
      G3 <- G[(2*N.size + 1):(4*N.size), ];
      Y3 <- YB[(2*N.size + 1):(4*N.size)];
      
      beta1 <- 0; beta2 <- 0; beta3 <- 0;
      for(ii in 1:tt) {
        beta1[ii] <- summary(lm(Y1 ~ G1[,ii]))$coefficients[2,1]
        beta2[ii] <- summary(lm(Y2 ~ G2[,ii]))$coefficients[2,1]
        beta3[ii] <- summary(lm(Y3 ~ G3[,ii]))$coefficients[2,1]
      }
      
      p1[i] <- summary(lm(beta2 ~ beta1))$coefficients[2,4];
      p2[i] <- summary(lm(beta3 ~ beta1))$coefficients[2,4];
    }
    #m1[trait, count] <- mean(p1 < 0.05);
    #m2[trait, count] <- mean(p2 < 0.05);
    res <- list();
    res$m1 <- p1
    res$m2 <- p2
    mylist <- rbind(mylist, res);
  }
}
