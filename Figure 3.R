###The blue module genes from SLE and turquoise module genes from AF were emerged and the he turquoise genes from SLE and turquoise module genes from AF were emerged respectively for PPI analysis
rm(list = ls()) 

load("AF_module.Rdata")
load("SLE_module.Rdata")

#SLE_blue和AF_turquoise合并
SLE_blue_and_AF_turquoise<-c(SLE_blue,AF_turquoise)
SLE_blue_and_AF_turquoise<-SLE_blue_and_AF_turquoise[!duplicated(SLE_blue_and_AF_turquoise)]
write.table(SLE_blue_and_AF_turquoise,
            file="SLE_blue_and_AF_turquoise.txt",
            row.names = F,
            col.names = F,
            quote = F)
#SLE_Turquoise和AF_turquoise
SLE_turquoise_and_AF_turquoise<-c(SLE_turquoise,AF_turquoise)
SLE_turquoise_and_AF_turquoise<-SLE_turquoise_and_AF_turquoise[!duplicated(SLE_turquoise_and_AF_turquoise)]
write.table(SLE_turquoise_and_AF_turquoise,
            file="SLE_turquoise_and_AF_turquoise.txt",
            row.names = F,
            col.names = F,
            quote = F)

###51 hub genes were identified.
rm(list = ls()) 

blue_betweenness<-blue_betweenness$Name
blue_closeness<-blue_closeness$Name
blue_degree<-blue_degree$Name
blue_stress<-blue_stress$Name

k1 = blue_betweenness%in% blue_closeness;table(k1)
blue_betweenness_closeness<- blue_betweenness[k1]

k2 = blue_betweenness_closeness%in% blue_degree;table(k2)
blue_betweenness_closeness_degree<- blue_betweenness_closeness[k2]

k3 = blue_betweenness_closeness_degree%in% blue_stress;table(k3)
blue_betweenness_closeness_degree_stress<- blue_betweenness_closeness_degree[k3]

write.table(blue_betweenness,
            file="blue_betweenness.txt",
            row.names = F,
            col.names = F,
            quote = F)
write.table(blue_closeness,
            file="blue_closeness.txt",
            row.names = F,
            col.names = F,
            quote = F)
write.table(blue_degree,
            file="blue_degree.txt",
            row.names = F,
            col.names = F,
            quote = F)
write.table(blue_stress,
            file="blue_stress.txt",
            row.names = F,
            col.names = F,
            quote = F)


turquoise_betweenness<-turquoise_betweenness$Name
turquoise_closeness<-turquoise_closeness$Name
turquoise_degree<-turquoise_degree$Name
turquoise_stress<-turquoise_stress$Name

k1 = turquoise_betweenness%in% turquoise_closeness;table(k1)
turquoise_betweenness_closeness<- turquoise_betweenness[k1]

k2 = turquoise_betweenness_closeness%in% turquoise_degree;table(k2)
turquoise_betweenness_closeness_degree<- turquoise_betweenness_closeness[k2]

k3 = turquoise_betweenness_closeness_degree%in% turquoise_stress;table(k3)
turquoise_betweenness_closeness_degree_stress<- turquoise_betweenness_closeness_degree[k3]

write.table(turquoise_betweenness,
            file="turquoise_betweenness.txt",
            row.names = F,
            col.names = F,
            quote = F)
write.table(turquoise_closeness,
            file="turquoise_closeness.txt",
            row.names = F,
            col.names = F,
            quote = F)
write.table(turquoise_degree,
            file="turquoise_degree.txt",
            row.names = F,
            col.names = F,
            quote = F)
write.table(turquoise_stress,
            file="turquoise_stress.txt",
            row.names = F,
            col.names = F,
            quote = F)

PPI<-c(blue_betweenness_closeness_degree_stress,
       turquoise_betweenness_closeness_degree_stress)
PPI<-PPI[!duplicated(PPI)]
write.table(PPI,
            file="PPI.txt",
            row.names = F,
            col.names = F,
            quote = F)
