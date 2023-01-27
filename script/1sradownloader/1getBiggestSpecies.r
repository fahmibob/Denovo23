library(dplyr)
library(readr)

species<-read.csv("./../data/species.txt")
sraData<-read.csv("./../data/rna_sra.csv")

sraSpecies<-sraData[sraData$Organism.Name %in% species$Organism.Name,]
sranotSpecies<-species[!species$Organism.Name %in% sraData$Organism.Name,]

sraSpeciesMax<-sraSpecies%>%
  group_by(Organism.Name)%>%
  top_n(1,Total.Size..Mb)
  
write.csv(sraSpeciesMax,"output.csv")
