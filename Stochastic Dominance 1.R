install.packages("RSD")
library(RSD)
library(ggplot2)
outcome1 <- c(5,3,5,7)
outcome2 <- c(2,4,6,8)
pr <- rep(1/4,4) #replicate
#Create SD objecte
sd.obj <- createStochasticDominance(outcome1,outcome2,pr,pr)
#Compare distributions based on FSD rule
fsd.test(sd.obj)
#Visualize the FSD values
fsd.plot(sd.obj)

# Another Example
outcome1 <- c(1,3,5,7,9)
outcome2 <- c(2,4,6,8)
pr1 <- rep(1/5,5)
pr2 <- rep(1/4,4)
# Create SD object
sd.obj <- createStochasticDominance(outcome1,outcome2,pr1,pr2)
# Compare distributions based on SSD rule
ssd.test(sd.obj) # 2
fsd.test(sd.obj)
# Visualize the SSD values
ssd.plot(sd.obj)

# Compare distributions based on ASSD rule (LL version)
assd.test(sd.obj, type = 'll')
# Compare distributions based on ASSD rule (THS version)
assd.test(sd.obj, type = 'ths')


