install.packages("igraphdata") # Red de interacción proteína-proteína (PPI)
library(igraphdata)
library(igraph)

#Cargar el dataset
data (yeast) 

proteinas <- igraph::upgrade_graph(yeast) # Actualizacion

red<-graph_from_adjacency_matrix()
plot(proteinas,vertex.size=15,vertex.size=5,
     edge.arrow.size=0.25,layout=layout_nicely,vertex.size.label=0.25)
#################################################################################
#                               CALCULOS
#################################################################################

# Tipo de distribución de conectividades que tiene
hist(degree(proteinas), breaks=10, col="violet", main="Red de PPI", )

#################################################################################
# Ajuste en log-log para ver que tipo de distribución podría ser.
# install.packages("poweRlaw")
library(poweRlaw)

deg_prot <- degree(proteinas) # Grados de la red
logarit <- displ$new(deg_prot) # Crear modelo

est_xmin <- estimate_xmin(logarit)# Estimar xmin
logarit$setXmin(est_xmin)

logarit$setPars(estimate_pars(logarit)) # Estimar y fijar los parámetros de la ley de potencia

plot(logarit) #Grafico
lines(logarit, col = "red")

#################################################################################

# Las diez proteínas más conectadas
sort(degree(proteinas), decreasing = TRUE)[1:10]

#################################################################################

# Diámetro y promedio de las distancias
diameter(proteinas) # Diametro

mean_distance(proteinas) # Promedio de la distancia

#################################################################################
#                               FUNCIONES
#################################################################################
proteinas1 <- igraph::upgrade_graph(yeast) # Actualizacion
# Q a partir de eliminar al azar un noodo de la red, genere el promedio de las distancias después de eliminar n=1,2 - 100 nodos al azar

distancia <- matrix(nrow = 100, ncol = 2) # Para guardar los datos (matriz vacia)
colnames(distancia) <- c("Nodo_elim", "Distancia") 

nodo <- sample(1:100, replace = FALSE) # Lista de nodos

nodos_i <- vcount(proteinas1) # Numero de nodos inicial

for (x in 1:100) {
  proteinas1 <- delete_vertices(proteinas1, nodo[x]) # Iniciar con la eliminacion de nodos
  distancia[x,] <- rbind(c(nodo[x], mean_distance(proteinas1)))# Registro de la distancia y nodo eliminado
  }

distancia
print(paste("Num de nodos al inicio:", nodos_i, "; Nodos al final", vcount(proteinas1))) # Comprobacion.

#################################################################################
proteinas2 <- igraph::upgrade_graph(yeast) # Actualizacion
# Q elimine las proteínas más conectadas y calcule el promedio de las distancias cada vez que se remueve un nodo

hubs <- which(degree(proteinas2) >= sort(degree(proteinas2), decreasing = TRUE)[100]) # Lista de nodos centrales
# Aquellos que su degree sea igual a alguno de los 100 valores mas altos

distancia2 <- matrix(nrow = length(hubs), ncol = 2) # Para guardar los datos (matriz vacia)
colnames(distancia2) <- c("Nodo_elim", "Distancia") 

for (x in 1:length(hubs)) {
  proteinas2 <- delete_vertices(proteinas2, hubs[x]) # Iniciar con la eliminacion de nodos
  distancia2[x,] <- rbind(c(hubs[x], mean_distance(proteinas2)))# Registro de la distancia y nodo eliminado
}

distancia2

 print(paste("Num de nodos al inicio:", nodos_i, "; Nodos al final:", vcount(proteinas2))) # Comprobacion.

#################################################################################
#                                 CALCULOS
#################################################################################
# Calcula el proemdio del coeficiente de clusterización. 
transitivity(proteinas)

# ¿Hay proteínas que tengan un coeficiente de clusterización de 1? Eso qué significa.
transitivity(proteinas, type = "local") # CC local para todos los nodos
cc1 <- which(transitivity(proteinas, type = "local") == 1) # Nodos con un CC = 1
cc1
