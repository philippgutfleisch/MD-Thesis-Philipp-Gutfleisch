#consensus cluster

df <- as.matrix(matrix)
df <- t(df)
df <- log(df+1)
df <- scale(df)
df <- t(df)

library(ConsensusClusterPlus)
results <-  ConsensusClusterPlus(df,
                               maxK=10,
                               reps=1000,
                               pItem=0.8,
                               pFeature=1,
                               clusterAlg="km",
                               distance="euclidean",
                               seed=123,
                               plot= "png")


#output
#get cases of each cluster               
cluster <- as.data.frame(results[[3]][["consensusClass"]])
cluster <- cluster%>%
        rownames_to_column("Patient_ID")
cluster <- merge(column_annotations, cluster)

write.table(cluster,file="consensuscluster.txt",sep="\t",quote=F,col.names=F)



icl <-  calcICL(results,plot="png")


                                                       



