#--------------------------------------------------
#--- Variable Selection for High-dimensional Nodal Attributes 
#--- in Social Networks with Degree Heterogeneity
#--------------------------------------------------
#--------- Real data analysis (International trade network)
#--------- data pre-processing: create X and Y
#--------------------------------------------------

#--------------------------------------------------
#--------- Construct node-level features X based on WDIs
#--------------------------------------------------

country_codes <- read.csv("~/trade_data/country_codes_V202001.csv")
# remove countries without WDI
library(readxl)
WDI <- read_excel("~/trade_data/WDI.xlsx", sheet = "Data")
WDI_variable_name= unique(WDI$`Series Code` )
WDI$`Country Code`=as.factor(WDI$`Country Code`)
WDI_country_code = as.character(unique(WDI$`Country Code`))
WDI_vec = as.numeric(sub("NA", "", WDI$`2017 [YR2017]`))
WDI_mat = matrix(WDI_vec, ncol = 1429, nrow=264, byrow = TRUE)
WDI_df = as.data.frame(WDI_mat, row.names = WDI_country_code )
names(WDI_df) = WDI_variable_name
country_code_overlap = intersect(WDI_country_code, country_codes$iso_3digit_alpha)
ind_net = match(country_code_overlap, country_codes$iso_3digit_alpha)
country_codes = country_codes[ind_net,]
# remove countries not in the network
n = dim(country_codes)[1]
WDI_df_after = data.frame()
for (i in 1:n){
  code = as.character(country_codes$iso_3digit_alpha[i])
  WDI_df_after[i, 1] = code 
  WDI_df_after[i, 2:1430] = WDI_df[which(WDI_country_code== code),] 
}
names(WDI_df_after)= c('country_code', WDI_variable_name)
WDI_df = WDI_df_after

# remove NAs in the data frame
# left with 142 countries and 328 variables
while (sum(apply(WDI_df, 1, function (x) sum(is.na(x))))>0){
  ratio = dim(WDI_df)[2]/dim(WDI_df)[1]
  ind_vb = which.max(apply(WDI_df, 2, function (x) sum(is.na(x))))
  count_vb = max(apply(WDI_df, 2, function (x) sum(is.na(x))))
  ind_cy = which.max(apply(WDI_df, 1, function (x) sum(is.na(x))))
  count_cy = max(apply(WDI_df, 1, function (x) sum(is.na(x))))
  if ( count_vb*ratio*1.5 > count_cy){
    # remove variable
    WDI_df = WDI_df[, -ind_vb ]
  }else{
    # remove country
    WDI_df = WDI_df[ -ind_cy ,  ]
  }
}
# save the final list of countries and X
ind  = match(WDI_df$country_code, country_codes$iso_3digit_alpha)
country_codes = country_codes[ind,]
write.csv(country_codes, file = '~/trade_data/final_country_codes_n=142.csv')
write.csv(WDI_df[,-1], file = '~/trade_data/trade_X_n=142.csv')

#--------------------------------------------------
#--------- Construct the adjacency matrix Y based on trade values
#--------------------------------------------------

BACI_Y2017 <- read.csv("~/trade_data/BACI_HS12_V202001/BACI_HS12_Y2017_V202001.csv") 
# downloaded from http://www.cepii.fr/CEPII/en/bdd_modele/download.asp?id=37 -> Trade flows data -> HS12 (2012-2019, 799 Mo)
country_codes_IE <- read.csv("~/trade_data/country_codes_V202001.csv")
code_IE=country_codes_IE$country_code[match(country_codes$iso_3digit_alpha,country_codes_IE$iso_3digit_alpha)]
n=142
# aggregate the total value of cash for edge i->j in Year 2017
Y_2017=matrix(0, ncol =n, nrow=n)
v = BACI_Y2017[1,5]
end_loop = dim(BACI_Y2017)[1]-1
for (ind in 1:end_loop){
  index = match(BACI_Y2017[ind,c(2,3)], code_IE)
  # if edge i->j belongs to the list of countries
  if (sum(is.na(index))==0){
    if (sum(BACI_Y2017[ind,c(2,3)]==BACI_Y2017[ind+1,c(2,3)])==2){
      v = v + BACI_Y2017[ind+1,5]
    }else{ 
      Y_2017[index[1], index[2] ]=v
      v = BACI_Y2017[ind+1,5]
    }
    if (ind == end_loop){
      index = match(BACI_Y2017[ind+1,c(2,3)], code_IE)
      if (sum(is.na(index))==0) {    
          Y_2017[index[1], index[2] ]=v
          }
    }
  }
}
write.csv(Y_2017 ,file = '~/trade_data/Y_2017.csv')
# constrct the adjacency matrix Y based on trade value Y_2017
Y = matrix(0, ncol = n, nrow = n)
for (i in 2:n){
  for (j in 1:i){
    if ((Y_2017[i,j]/rowSums(Y_2017)[i]>=0.005)&(Y_2017[j,i]/rowSums(Y_2017)[j]>=0.005)){
      Y[i,j]=1
      Y[j,i]=1
    }
  }
}
sum(Y)
write.csv(Y,file = '~/trade_data/trade_X_n=142.csv')
