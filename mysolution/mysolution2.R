
set.seed(42)
g <- sample_pa(1000)

layout <- layout.fruchterman.reingold(g)

plot(g,
  layout = layout, vertex.size = 2,
  vertex.label = NA, edge.arrow.size = .2
)

bt <- betweenness(g)

centr <- which.max(bt)

centr

# ma indeks nr 2

diameter(g)

# średnica 10

# Barabási-Albert (BA) i Erdős-Rényi (ER) różnia się tym, że BA jest mniej losowy i tworzy zgrupowania node'ów, co lepiej odzwierciedla prawidziwe social networks. ER jest całkowicie losowy i nie widać tam żadnych schematów.