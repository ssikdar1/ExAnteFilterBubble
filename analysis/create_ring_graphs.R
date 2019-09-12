library(viridis)
library(tidygraph)
library(ggraph)
library(igraph)

N = 10
ring <- create_ring(N)

hop_distance <- function(i,j,N){min(abs(i-j), abs(j-i), abs(j-i-N), abs(i-j-N))}

viridis_pallete <- viridis::viridis(N)
viridis_pallete_rev <- rev(viridis_pallete)

white_pallete <- c()
for (i in 1:10){white_pallete <- c(white_pallete, "white")}

png("circle1.png")
plot(ring, 
     vertex.size=10, vertex.label=NA, layout=layout_in_circle, vertex.color = white_pallete)
dev.off()

white_pallete[N] <- viridis_pallete[N]

png("circle2.png")
plot(ring, 
     vertex.size=10, vertex.label=NA, layout=layout_in_circle, vertex.color = white_pallete)
dev.off()

new_pal <- c()
for (i in 1:10){new_pal <- c(new_pal, viridis_pallete_rev[hop_distance(1,i,10)])}
new_pal[1] <- viridis_pallete[N]
png("circle3.png")
plot(ring, vertex.size=10, vertex.label=NA, layout=layout_in_circle, vertex.color = new_pal)
dev.off()



