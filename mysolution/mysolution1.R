library(igraph)


g <- sample_gnp(100, 0.05)


summary(g)

# IGRAPH 8d67f24 U--- 100 231 -- Erdos-Renyi (gnp) graph
# + attr: name (g/c), type (g/c), loops (g/l), p (g/n)

# graf nie jest ważony


# wierzchołki
V(g)

# krawędzie
E(g)

E(g)$weight <- runif(length(E(g)), 0.01, 1)

summary(g)
# > summary(g)
# IGRAPH 81e826b U-W- 100 267 -- Erdos-Renyi (gnp) graph
# + attr: name (g/c), type (g/c), loops (g/l), p (g/n), weight (e/n)

# jest W czyli graf jest ważony

# stopień każdego węzła
degree(g)

# histogram
hist(degree(g))

cl <- components(g)
cl
# są dwa klastry, jeden ma 99 node'ów a drugi 1 node

pr <- page_rank(g)$vector

plot(g, vertex.size = pr * 200)
