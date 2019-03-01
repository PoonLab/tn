library(igraph)

#Expecting the output from a tn93 run formatted to a csv file.
#Expecting patient information in the format ID_Date
args = commandArgs(trailingOnly = T)
input <- read.csv(args[1], stringsAsFactors = F)

#Creates a graph based on the inputted data frame. The tn93 Distances become edge4 attributes
g <- graph_from_data_frame(input, directed=FALSE, vertices=NULL)

#Adds the ID's and Sample collection years as different vertex attributes for each vertex
####- TO-DO: Assign Dates by Month -####
temp <- sapply(V(g)$name, function(x) strsplit(x, '_')[[1]])


y <- as.Date(temp[2,])
yAlt <- as.Date(temp[2,], format = "%d-%b-%y")
y[is.na(y)] <- yAlt[!is.na(yAlt)]

V(g)$name <- temp[1,]
V(g)$year <- as.integer(y)

#Obtain the range of years and the maximum input year
years <- as.integer(levels(factor(V(g)$year)))
newY <- (max(years) - 365)
newV <- V(g)[year>=newY]
oldV <- V(g)[year<newY]

ageD <- {}

for (i in cutoffs){
  print(i)
  subG <- subgraph.edges(g, E(g)[(newV%--%oldV) & (Distance <= i)], delete.vertices = F)
  
  temp <- sapply(oldV, function(x) {
    case <- V(subG)[x]
    weight <- length(E(subG)[inc(case)])
    c(weight, case$year)
  })
  
  df <- data.frame(Age=unname(temp[2,]) , Freq =unname(temp[1,])) 
  mod <- glm(Freq ~ Age, data = df, family = "poisson")
  
  ageD[as.character(i)] <- mod
}
