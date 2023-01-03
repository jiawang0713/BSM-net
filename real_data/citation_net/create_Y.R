#--------------------------------------------------
#--- Variable Selection for High-dimensional Nodal Attributes 
#--- in Social Networks with Degree Heterogeneity
#--------------------------------------------------
#--------- Real data analysis (Facebook friendship network)
#----------   Construct the adjacency matrix Y based on citation data set https://www.aminer.cn/citation
#--------------------------------------------------

#--------------------------------------------------
#--------- Step 1: convert raw text data 'acm.txt' into 'citation_full.RData' (file contained in the folder)
#--------------------------------------------------
setwd("~/citation_network_V9/")
con = file("~/acm.txt", "r") # 'acm.txt' could be downloaded from https://www.aminer.cn/citation ACM-Citation-network V9
num=2385022 
title<-rep("", num)
author<-rep("", num)
year<-rep(0, num)
publication<-rep("", num)
citation<-rep("", num)
year<-rep(0, num)
index<-rep(0, num)
abstract<-rep("", num)

count=1
while (count<(num+1)) {
  line = readLines(con, n = 1)
  if(line ==""){
    if(citation[count]!="")
    {
      citation[count]=substring(citation[count], 2)
    }
    count=count+1
    
  } else if(substring(line, 1,2)=="#*") 
  {
    title[count]=substring(line, 3)
  } else if(substring(line, 1,2)=="#@") {
    author[count]=substring(line, 3)
  } else if(substring(line, 1,2)=="#t") {
    year[count]=as.numeric(substring(line, 3))
  } else if(substring(line, 1,2)=="#c") {
    publication[count]=substring(line, 3)
  } else if(substring(line, 1,2)=="#%") {
    citation[count]=paste(citation[count],substring(line, 3), sep=" ")
  } else if(substring(line, 1,6)=="#index") {
    index[count]=as.numeric(substring(line, 7))
  } else if(substring(line, 1,2)=="#!") {
    abstract[count]=substring(line, 3)
  } 
}
close(con)
save(title, author, year, publication, citation, index, abstract,file = "~/citation_full.RData")

#--------------------------------------------------
#--------- Step 2: Get node degree for each paper and save it in a vector named 'count' (file contained in the folder)
#--------------------------------------------------

indx_sub = which(abstract!="")
temp = which(citation!="")
indx_sub = intersect(indx_sub,temp)
current_indx = indx_sub
current_n = length(current_indx)
current_edge_list = citation[current_indx]

# get degree for each node, save in a vector named 'count'
count = rep(0,current_n)
for (j in 1:current_n){
  temp =  as.numeric(unlist(strsplit(current_edge_list[j], split=" "))) + 1 
  count[j] = count[j]+length(temp)  
  ind = match(temp, current_indx)
  ind = ind[!is.na(ind)]
  if (length(ind)>0){count[ind] = count[ind] +1 }
}
write.csv(count,file = '~/count.csv')

#--------------------------------------------------
#--------- Step 3: Get a dense sub-network with samller size 
#--------------------------------------------------

count = read.csv("~/count.csv")
count = count$x
ind = which(count >200)
selected_indx = current_indx[ind]
selected_edge_list = current_edge_list[ind]
Y = matrix(0, ncol= length(selected_indx), nrow= length(selected_indx))
for (i in 1:length(selected_edge_list)){
  temp<- as.numeric(unlist(strsplit(as.vector(selected_edge_list[i]), split=" "))) + 1 
  temp = intersect(temp, selected_indx)
  if (length(temp)>0){
    ind = match(temp, selected_indx)
    Y[i, ind] = rep(1, length(ind))
    Y[ind, i ] = rep(1, length(ind))
  }
}
n = length(selected_indx)
write.csv(Y,file = '~/citation_Y_1207_n=476.csv')
write.csv(selected_indx,file = '~/selected_indx_1207_n=476.csv')
