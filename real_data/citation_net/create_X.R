#--------------------------------------------------
#--- Variable Selection for High-dimensional Nodal Attributes 
#--- in Social Networks with Degree Heterogeneity
#--------------------------------------------------
#--------- Real data analysis (Facebook friendship network)
#----------   LDA on paper abstract to construct X
#--------------------------------------------------

library(GGally)
library(network)
library(sna)
setwd("~/citation_network/")
load("citation_full.RData") # created in 'create_Y.R'
set.seed(2020)

selected_indx = read.csv(file = 'selected_indx_1207_n=476.csv')
selected_indx = as.vector(selected_indx$x)
selected_indx_sub = match(selected_indx, indx_sub)

# select more documents to build LAD
sample_ind = sample( (1:length(indx_sub))[-selected_indx_sub] , 5000-n)
sample_ind = c(selected_indx_sub, sample_ind)

# preprocess 
library(tm)
corpus<-abstract[indx_sub][sample_ind]
docs <- Corpus(VectorSource(unlist(corpus)))
docs = DocumentTermMatrix(docs, control = list(stemming = TRUE, stopwords = stopwords("english") , minWordLength = 3, removeNumbers = TRUE, removePunctuation = TRUE))

library("slam")
summary(col_sums(docs))
# term frequency-inverse document frequency
term_tfidf = tapply(docs$v/row_sums(docs)[docs$i], docs$j, mean) *log2(nDocs(docs)/col_sums(docs > 0))
docs <- docs[, term_tfidf >= 0.1]
ind_deleted = which(row_sums(docs) ==0)
docs<- docs[row_sums(docs) > 0,]
summary(col_sums(docs))
dim(docs)

library("topicmodels")
k = 500 # the number of topics
SEED = 0713
LDA.Gibbs.fit =  LDA(docs, k = k, method = "Gibbs", control = list(seed = SEED, burnin = 2000, thin = 100, iter = 1000))

Terms <- terms(LDA.Gibbs.fit , 5) # The 5 most frequent terms for each topics for interpretation. 

write.csv(Terms,file = '~/Terms_1207_n=476.csv')
write.csv(posterior(LDA.Gibbs.fit)$topic[1:n,],file = '~/citation_X_1207_n=476.csv')
