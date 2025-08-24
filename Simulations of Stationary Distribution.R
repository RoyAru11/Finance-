library(ggplot2)
# simulate discrete Markov chains according to transition matrix P
run.mc.sim <- function( P, num.iters = 50 ) {
  
  # number of possible states
  num.states <- nrow(P)
  
  # stores the states X_t through time
  states     <- numeric(num.iters) #array of size 50
  
  # initialize variable for first state 
  states[1]    <- 1
  
  for(t in 2:num.iters) {
    
    # probability vector to simulate next state X_{t+1}
    p  <- P[states[t-1], ]
    
    ## draw from multinomial and determine state
    states[t] <-  which(rmultinom(1, 1, p) == 1)
  }
  return(states)
}
P <- t(matrix(c( 0.7, 0.2, 0.1, 
                 0.4, 0.6,   0, 
                 0,   1,   0  ), nrow=3, ncol=3))

num.chains     <- 5
nm.iterations  <- 50

# each column stores the sequence of states for a single chains
chain.states  <- matrix(NA, ncol=num.chains, nrow=num.iterations)

# simulate chains
for(c in seq_len(num.chains)){
  chain.states[,c] <- run.mc.sim(P)
}

matplot(chain.states, type='l', lty=1, col=1:5, ylim=c(0,4), xlim=c(0,50), ylab='state', xlab='time')
abline(h=1, lty=3)
abline(h=3, lty=3)

P <- t(matrix(c( 1/3, 2/3,   0,   0,  0,   0,   0,   0,
                 1/3, 1/3, 1/3,   0,  0,   0,   0,   0,
                 0, 1/3, 1/3, 1/3,  0,   0,   0,   0,
                 0,   0, 1/3, 1/3, 1/3,  0,   0,   0,
                 0,   0,   0, 1/3, 1/3, 1/3,  0,   0,
                 0,   0,   0,   0, 1/3, 1/3, 1/3,  0,
                 0,   0,   0,   0,   0, 1/3, 1/3, 1/3,
                 0,   0,   0,   0,   0,   0, 2/3, 1/3), nrow=8, ncol=8))
num.chains     <- 5
num.iterations <- 50
chain.states <- matrix(NA, ncol=num.chains, nrow=num.iterations)
for(c in seq_len(num.chains)){
  chain.states[,c] <- run.mc.sim(P)
}
matplot(chain.states, type='l', lty=1, col=1:5, ylim=c(0,10), xlim=c(0,50), ylab='state', xlab='time')
abline(h=1, lty=3)
abline(h=8, lty=3)


# Simulating Discrete Markov Chains: Limiting Disributions
run.mc.sim <- function( P,   # probability transition matrix
                        num.iters=50, 
                        num.chains=150 )
{
  
  # number of possible states
  num.states <- nrow(P)
  
  # states X_t for all chains
  states     <- matrix(NA, ncol=num.chains, nrow=num.iters)
  
  # probability vectors pi^n through time
  all_probs  <- matrix(NA, nrow=num.iters, ncol=num.states)
  
  # forces chains to start in state 1
  pi_0      <- c(1, rep(0, num.states-1))
  
  # initialize variables for first state 
  P_n           <- P
  all_probs[1,] <- pi_0
  states[1,]    <- 1
  
  for(t in 2:num.iters) {
    
    # pi^n for this iteration
    pi_n           <- pi_0 %*% P_n
    all_probs[t,]  <- pi_n
    
    for(chain_num in seq_len(num.chains)) {
      # probability vector to simulating next state 
      p                     <- P[ states[t-1,chain_num], ]
      states[t,chain_num]   <- which(rmultinom(1, 1, p) == 1)
    }
    
    # update probability transition matrix
    P_n           <- P_n %*% P
  }
  return(list(all.probs=all_probs, states=states))
}

# setup transition matrix 
P <- t(matrix(c( 0.7, 0.2, 0.1, 
                 0.4, 0.6,   0, 
                 0,   1,   0  ), nrow=3, ncol=3))

sim1 <- run.mc.sim(P)
states <- sim1[[2]]
matplot(states[,1:5], type='l', lty=1, col=1:5, ylim=c(0,4), ylab='state', xlab='time')
abline(h=1, lty=3)
abline(h=3, lty=3)

all.probs <- sim1[[1]]
matplot(all.probs, type='l', col=1:3, lty=1, ylab='probability', xlab='time')
legend('topright', c('state.1', 'state.2', 'state.3'), lty=1, col=1:3)

A <- matrix(c(-0.3, 0.2, 0.1, 1, 0.4, -0.4, 0, 1, 0, 1, -1, 1 ), 
                    ncol=3,nrow=4)
b <- c(0,0,0, 1)
pi <- drop(solve(t(A) %*% A, t(A) %*% b)) #t means transpose

# show comparison
results1           <- t(data.frame(pi_n = all.probs[50,], pi = pi))
colnames(results1) <- c('state.1', 'state.2', 'state.3')
results1

# Finally, we can also plot the proportion of chains that are in each state through time.
# These should roughly equal the probability vectors above, with some noise due to random chance:
state.probs <- t(apply(apply(sim1[[2]], 1, function(x) table(factor(x, levels=1:3))),
                       2, function(x) x/sum(x)))
matplot(state.probs[1:50,], col=1:3, lty=1, type='l', ylab='empirical probability', xlab='time')
legend('topright', c('state.1', 'state.2', 'state.3'), lty=1, col=1:3)


# Simulation 2: 8x8 example
P <- t(matrix(c( 1/3, 2/3,   0,   0,  0,   0,   0,   0,
                 1/3, 1/3, 1/3,   0,  0,   0,   0,   0,
                 0, 1/3, 1/3, 1/3,  0,   0,   0,   0,
                 0,   0, 1/3, 1/3, 1/3,  0,   0,   0,
                 0,   0,   0, 1/3, 1/3, 1/3,  0,   0,
                 0,   0,   0,   0, 1/3, 1/3, 1/3,  0,
                 0,   0,   0,   0,   0, 1/3, 1/3, 1/3,
                 0,   0,   0,   0,   0,   0, 2/3, 1/3), nrow=8, ncol=8))

sim2a <- run.mc.sim(P)
states <- sim2a[[2]]
matplot(states[,1:5], type='l', lty=1, col=1:5, ylim=c(0,9), ylab='state', xlab='time')
abline(h=1, lty=3)
abline(h=8, lty=3)

all.probs <- sim2a[[1]]
matplot(all.probs, type='l', col=1:8, lty=1, ylab='probability',
        xlab='time', ylim=c(0, 0.5))
legend('topright', paste('state.', 1:8, sep=''), lty=1, col=1:8)

# Now we alter the transition matrix above to encourage the chain to stay in states 4 and 5:
P <- t(matrix(c( 1/3,   2/3,    0,    0,    0,    0,    0,   0,
                 1/3,   1/3,  1/3,    0,    0,    0,    0,   0,
                 0,  .5/6, .5/6,  5/6,    0,    0,    0,   0,
                 0,     0, .5/6,  5/6, .5/6,    0,    0,   0,
                 0,     0,    0, .5/6,  5/6, .5/6,    0,   0,
                 0,     0,    0,    0,  5/6, .5/6, .5/6,   0,
                 0,     0,    0,    0,    0,  1/3,  1/3, 1/3,
                 0,     0,    0,    0,    0,    0,  2/3, 1/3 ), nrow=8, ncol=8))
sim2b <- run.mc.sim(P)
all.probs <- sim2b[[1]]
matplot(all.probs, type='l', col=1:8, lty=1, ylab='probability',
        xlab='time', ylim=c(0,1.2))
legend('topright', paste('state.', 1:8, sep=''), lty=1, col=1:8)

state.probs <- t(apply(apply(sim2b[[2]], 1, function(x) table(factor(x, levels=1:8))), 2, function(x) x/sum(x)))
matplot(state.probs[1:50,], col=1:8, lty=1, type='l', ylab='empirical probability', xlab='time', ylim=c(0,1.2))
legend('topright', paste('state.', 1:8, sep=''), lty=1, col=1:8)










