#--------------------------------------------------
#--- Variable Selection for High-dimensional Nodal Attributes 
#--- in Social Networks with Degree Heterogeneity
#--------------------------------------------------
#--------- Real data analysis (Facebook friendship network)
#--------- plot the results
#--------- last edited in July 22, 2021
#--------------------------------------------------
library(readxl)
library(GGally)
library(network)
library(sna)
# read the results
d = 21 
n = 347
selected_sp <- read.table("~/Desktop/ACI-ICS/Research_NET/facebook_net/gamma_sp_abs_d=21_n=347_M=5000.csv", quote="\"", comment.char="")
beta_delta_sp <- read.table("~/Desktop/ACI-ICS/Research_NET/facebook_net/beta_delta_sp_abs_d=21_n=347_M=5000.csv", quote="\"", comment.char="")
beta_delta_sp = as.matrix(beta_delta_sp)
rank = which(apply(selected_sp, 2, mean)>0.5)
selected_sp_vb = rank[which(rank<d+1)] 
selected_nodes = rank[which(rank>d)] - d 
beta_sp = colMeans(beta_delta_sp)[as.numeric(selected_sp_vb)]

colMeans(selected_sp)[as.numeric(selected_sp_vb)]
for (i in selected_sp_vb){
  print(quantile(beta_delta_sp[,i], probs = c(0.025, 0.975)))
}

# plot Facebook friendship network 
Y =  read.csv("Documents/Research/BVS_NET/fb_net/Y.csv", row.names=1)
net = network(Y, directed = FALSE)
net %v% "group"  = col
set.seed(827)
ggnet2(net, color = '#ff9b00', 
       size = "degree", size.cut = 3 , max_size = 2.5, edge.alpha = 1/4, alpha=0.5)+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12))

# plot Facebook friendship network with active nodes marked in orange
net = network(Y, directed = FALSE)
col=rep("inactive",n)
col[selected_nodes]=rep("active",length(selected_nodes))
net %v% "group"  = col
set.seed(827)
ggnet2(net, color = 'group' , palette = c("inactive" = "#a8bed1", "active" = "#ff9b00"), 
       size = "degree", size.cut = 3 , max_size = 2.5, edge.alpha = 1/4, alpha=0.5)+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12))

# plot histograms showing degree distributions
df = data.frame(degree = c(degree(net)[selected_nodes], degree(net)[-selected_nodes])/2,
                group = c(rep('active', length(selected_nodes)), rep('inactive', n-length(selected_nodes)))   )
ggplot(df, aes(x=degree, color=group, fill=group)) + 
  geom_histogram(aes(y=..density..), alpha=0.5, bins = 30,
                 position="identity")+
  geom_density(alpha=.2)+
  scale_color_manual(values=c( "#ff9b00","#a8bed1"))+
  scale_fill_manual(values=c( "#ff9b00","#a8bed1"))+
  theme_bw()