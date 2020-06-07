setwd("C:/google drive/Coding workspace/works/14.COVID19-p2_April/COVID2_Shared/motif analysis")


#adj_list1=make_kmer_network()
print('select the kmer dataset')
fname=file.choose()
kmer_result=read_csv(fname)
kmer=kmer_result



make_kmer_network=function(country,month,kmer)
{
  #it will create a network from a distance matrix to a binary matrix to be considered as 
  #adjacency matrix of the kmer_network. The cutoff value is chosen the median of the input values
  library(igraph)
  library(readr)
  library(gdata)
  library(philentropy)
  #library(distances)
  
  mon_seq=which(kmer$Months==month)
  kmer=kmer[mon_seq,]
  
  contry_seq=which(kmer$Country==country)
  if (length(contry_seq)>100)
  { contry_seq=sample(contry_seq,size = 75)}
  kmer=kmer[contry_seq,]
  kmerdata=kmer[,1:84]
  kmermeta=kmer[,85:91]
  
  kmer_dist=distance(kmerdata, method = "euclidean",use.row.names = TRUE)
  cut_off=median(kmer_dist,diag=FALSE)
  #cut_off=0.01
  hist(kmer_dist)
  
  t1=which(kmer_dist<cut_off)    #make those 1
  t2=which(kmer_dist>=cut_off)   #make those 0
  
  kmer.binary=kmer_dist #matrix(data = NA,nrow = dim(kmerdata)[1],ncol = dim(kmerdata)[1])
  kmer.binary[t1]=0
  kmer.binary[t2]=1
  warning('similarity value less than cutoff is converted to 1, NOTE this')
  rm(t1,t2)
  
  adj_list=gather(as.data.frame(kmer.binary))
  adj_list=adj_list[which(adj_list$value==1),]
  adj_list=as.matrix(adj_list)
  adj_list[,1]=substring(adj_list[,1],2,5)  #delete 'V' from the column 'key'
  adj_list=adj_list[,-2]
  write.csv(adj_list,file = paste0(country,"_",month,"_network.csv"))
  #adj_list=separate(adj_list,key,c("V1,V2"),sep = "")
  
  return(adj_list)
  
}














  #plot(density(dt.mat))
  #hist(dt.mat)
# 
#   g1=graph_from_adjacency_matrix(kmer.binary,mode="undirected",diag=FALSE)
#   plot(g1)
#   
#   ######make two col adj list from dt.binary for use in FANMOD
#   row_vec=NULL
#   col_vec=NULL
#   for (i in 1:dim(kmer.binary)[1])
#   {
#     t1=which(kmer.binary[i,]==1)
#     if(length(t1>0))
#     {
#       col_vec=c(col_vec,t1)
#       row_vec=c(row_vec,i)
#     }
#   }
#   kmer.binary.list=expand.grid(row_vec,col_vec)
#   #write.table(expand.grid(row_vec,col_vec),file = "graph.txt",row.names = FALSE,col.names = FALSE,sep = " ")
#   write.table(kmer.binary.list,file = "graph.txt",row.names = FALSE,col.names = FALSE,sep = " ")
#   
#   ######make two col adj list from dt.binary for use in FANMOD
#   dt.binary.list1=NULL
#   for (i in 1:dim(dt.binary)[1])
#   {
#     t1=which(dt.binary[i,]==1)
#     dt.binary.list1=rbind(dt.binary.list1,expand.grid(i,t1))
#   }
#   
#   write.table(dt.binary.list1,file = "China_graph.txt",row.names = FALSE,col.names = FALSE,sep = " ")
  
  ######################### date: 6th May 2020
  
  #g1.degree.histogram <- as.data.frame(table(degree(g1)))
  


#################
print('note distance 0 means same object, distance is complement of similarity here')
x=matrix(c(1,2,3,4,5,6,1,2,3),ncol=3,byrow = TRUE)
distance(x, method = "euclidean")