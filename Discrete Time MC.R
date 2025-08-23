P=t(matrix(c(c(0.5,0.4,0.1),c(0.3,0.4,0.3),c(0.2,0.3,0.5)),nrow=3))
# Check sum across = 1
apply(P,1,sum) 
x0=c(0.1,0.2,0.7)
sum(x0)
#After one step
x0%*%P
#The two-steo prob trans matrix
P%*%P
x0%*%P%*%P
# First time only (installs package)
install.packages("expm")
library("expm")
# Verify the second power is P%*%P
P%^%100
c(1,0,0)%*%(P%^%100)
c(0.2,0.5,0.3) %*%(P%^%100)

#STATIONARY DISTRIBUTION
qr.solve
A        <- matrix(c(-0.3, 0.2, 0.1, 1, 0.4, -0.4, 0, 1, 0, 1, -1, 1 ), ncol=3,nrow=4)
b        <- c(0,0,0, 1)
pi       <- qr.solve(A,b)
names(pi) <- c('state.1', 'state.2', 'state.3')
pi
sum(pi)
stationary <- function(transition) {
  stopifnot(is.matrix(transition) &&
              nrow(transition) == ncol(transition) &&
              all(transition >= 0 & transition <= 1))
  
  n <- nrow(transition)
  
  # Stationary distribution solves: pi * P = pi
  # Equivalent: t(P) - I, with constraint sum(pi)=1
  A <- rbind(t(transition) - diag(n), rep(1, n))
  b <- c(rep(0, n), 1)
  
  res <- qr.solve(A, b)
  names(res) <- paste0("state.", 1:n)
  return(res)
}

# Test
stationary(matrix(c(0.7, 0.2, 0.1,
                    0.4, 0.6, 0,
                    0,   1,   0),
                  nrow=3, byrow=TRUE))
#To check the formula 
pi <- stationary(P)
pi %*% P
pi

#Computing Stationary Distributions of a Discrete Markov Chain
P=t(matrix(c(c(0.5,0.4,0.1),c(0.3,0.4,0.3),c(0.2,0.3,0.5)),nrow=3))
# Check sum across = 1
apply(P,1,sum)

x0=c(0.1,0.2,0.7)
# Check sums to 1
sum(x0)

#Brute-force solution
pi_bru <- (P %^% 100)[1,]
pi_bru

pi_bru - pi_bru%*%P

#Solving ia eigendecomposition
library(MASS)
r=eigen(P) #gives both the eigenvectors and eigenvalues
rvec=r$vectors #Ax=Λx
lvec=ginv(r$vectors) #inverse of right eigen vec gives left eigenvec
lam <- r$values
# Two ways of checking the spectral decomposition
## Standard definition
rvec %*%diag(lam)%*%ginv(rvec) #P=RΛR^-1
## With left eigenvectors (trivial chang)
rvec%*%diag(lam)%*%lvec
lam

#Another Procedure of finding Stationary Distribution
pi_eig<-lvec[1,]/sum(lvec[1,])
pi_eig
sum(pi_eig)
pi_eig %*% P
r<-eigen(t(P))
V<-r$vectors
lam<-r$values
V%*%diag(lam)%*%ginv(V)
pi_eig2 <- V[,1]/sum(V[,1])

#Another approach to solve the simultaneous linear equations
K<-3
A_basic <- t(diag(rep(1,K))-P)
b_basic <- rep(0,K)

# Now add the constraint 
A_constr <- rbind(A_basic,rep(1,K))
b_constr <- c(b_basic,1)

pi_lineq <- t(solve(t(A_constr)%*%A_constr,t(A_constr)%*%b_constr))
pi_lineq%*%P
pi_lineq
sum(pi_lineq)
















