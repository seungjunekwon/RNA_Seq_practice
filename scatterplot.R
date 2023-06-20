library(ggplot2)

non <- read.csv("non.csv", header=T, row.names=1)
up <- read.csv("up.csv", header=T, row.names=1)
down <- read.csv("down.csv", header=T, row.names=1)

p <- ggplot(non, aes(x=aC, y=bA)) + geom_point(size=1, color="gray")

g <- ggplot(up, aes(x=aC, y=bA)) + geom_point(size=2, color="red")

plot <- ggplot(width=5, height=5) +
         geom_point(data=non, aes(x=aC, y=bA), size=1, color="grey77", alpha=0.3) + geom_point(data=up, aes(x=aC, y=bA), size=2, color="red3", alpha=0.7) + geom_point(data=down, aes(x=aC, y=bA), size=2, color="blue3", alpha=0.7)
plot + theme_classic()
