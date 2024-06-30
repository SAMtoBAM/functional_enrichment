

library(gtools)
library(ggplot2)
library(ggpubr)
library(dplyr)

##read in data
##data contains a factor for starship or core, the COG and the count of COG for each group

###First run ESE00019
COG=read.csv(file="~/cluster2/projects/Penicillium/functional_enrichment/ESE00019.combined_mRNA.COG.tsv", header=F, sep ='\t')

##rewrite another column to siplify core and starship categories
COG$dataset<- ifelse(grepl("starship", COG$V3), "starship", "core")
##have a look at a stacked column
ggplot()+geom_col(data=COG , aes(x=dataset, y=V2, fill=V1), position= position_stack())+theme_pubr(x.text.angle = 35)
##calculate the proportion per group (starship, core)
COG_prop= COG %>%  
  group_by(dataset, V1) %>%
  summarise(count = sum(V2)) %>% 
  mutate(proportion = count / sum(count))
##check out a stacked column of the proportions now
ggplot()+
  geom_col(data=COG_prop , aes(x=dataset, y=proportion, fill=V1), position= position_stack())+
  theme_pubr(x.text.angle = 35)


###begin generating values to be used in the hypergeomtric tests
xx = COG_prop %>%  
  group_by(V1) %>% summarise(m = sum(count))
##merge with the proportion dataset
COG_prop=merge(COG_prop, xx, by='V1')

##add the total of all other genes (not the same COG) (n)
xx=COG_prop %>%  
  reframe(n = sum(count)-m)
COG_prop$n=xx

##add the total of all starship genes (k)
##can now subset for just those lines with the starship counts; and then sum all
COG_prop2=COG_prop %>%
  filter(stringr::str_detect(dataset, 'starship'))
COG_prop2$k=sum(COG_prop2$count)


##now loop over each line and calculate the pvalue
##perform the hypergeometric COG for each COG using phyper(lower.tail = FALSE), and to COG both sides we COG q-1

##first create an empty dataframe
COG_prop3 = data.frame()
enrich_pval=""
derich_pval=""
##now loop through the lines and create an array of the resulting pvalues 
for(i in 1:nrow(COG_prop2)) {
  row <- COG_prop2[i,]
  # do stuff with row
  enrich_pval[[i]]=phyper(as.numeric(row$count)-1, as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = FALSE)
  derich_pval[[i]]=phyper(as.numeric(row$count), as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = TRUE)
}
##add the pvalues to the dataframe
COG_prop3 <- cbind(COG_prop2,enrich_pval)
COG_prop3 <- cbind(COG_prop3,derich_pval)
##now adjust the pvalues and then add that to another column
enrich_pval.adjust=p.adjust(enrich_pval, method = "BH")
derich_pval.adjust=p.adjust(derich_pval, method = "BH")

COG_prop3 <- cbind(COG_prop3,enrich_pval.adjust)
COG_prop3 <- cbind(COG_prop3,derich_pval.adjust)

##now add the adjusted pvalues, without differentiating between enriched and deriched
##we will not differentiate as the difference is clear and we are just showing statistically significant differences

COG_prop3$group1="core"
COG_prop3$group2="starship"

##add stars instead for the significance
COG_prop3$enrich_pval.adjust2 =stars.pval(COG_prop3$enrich_pval.adjust)
COG_prop3$derich_pval.adjust2 =stars.pval(COG_prop3$derich_pval.adjust)

ESE19unfilt=ggplot(data=COG_prop , aes(x=V1, y=proportion, fill=dataset))+
  geom_col(position= position_dodge())+
  theme_pubr()+
  xlab("COG category")+
  ylab("proportion")+theme(legend.title = element_blank())+
  labs(subtitle=expression(~italic(P.caseifulvum)~' ESE00019'))+
  scale_fill_manual(values = c("grey","red"))+
  stat_pvalue_manual(COG_prop3, label = "enrich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  stat_pvalue_manual(COG_prop3, label = "derich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  scale_x_discrete(labels=c("A"=expression(bold(A)), "Absent"=expression(bold(Absent)), "C"=expression(bold(C)),"E"=expression(bold(E)),"F"=expression(bold(F)),"G"=expression(bold(G)),"H"=expression(bold(H)),"J"=expression(bold(J)),"L"=expression(bold(L)),"O"=expression(bold(O)),"Q"=expression(bold(Q)),"T"=expression(bold(T)),"U"=expression(bold(U)),"V"=expression(bold(V)), parse=TRUE))



##calculate the proportion of each COG annotation with the core and starship genes seperately
##also re;ove any genes with COGs S (unknown), and L/K (transposable element associated). The latter is because the annotated regions were not softmasked before annotation due to size and databases
##also removed all the 'Absent' (genes with no COG annotation)
COG_prop4= subset(COG, V1 != "Absent" & V1 != "S" & V1 != "L" & V1 != "K") %>%  
  group_by(dataset, V1) %>%
  summarise(count = sum(V2)) %>% 
  mutate(proportion = count / sum(count))

##plot a stacked column with all the COG proportions for core vs starships
ggplot()+
  geom_col(data=COG_prop4 , aes(x=dataset, y=proportion, fill=V1), position= position_stack())+
  theme_pubr(x.text.angle = 35)

##can plot side by side proportions for each COG seperately
ggplot()+
  geom_col(data=COG_prop4 , aes(x=V1, y=proportion, fill=dataset), position= position_dodge())+
  theme_pubr()+
  xlab("COG category")+
  ylab("proportion")+theme(legend.title = element_blank())+
  scale_fill_manual(values = c("grey","red"))

##now need to calculate the hypergeometric COG p-value for all COG groups and then correct for all the COGs
##the four features needed for the hypergeometric COG are q, m, n, k = number genes in starships with the COG, the total number of genes (in and out of starships) with COG, the total count of all other genes (not with the same COG), the total of all starship genes 
##we will actually test for enrichment and deenrichment too

##get the sum of all genes with each COG (m)
xx = COG_prop4 %>%  
  group_by(V1) %>% summarise(m = sum(count))
##merge with the proportion dataset
COG_prop4=merge(COG_prop4, xx, by='V1')

##add the total of all other genes (not the same COG) (n)
xx=COG_prop4 %>%  
  reframe(n = sum(count)-m)
COG_prop4$n=xx

##add the total of all starship genes (k)
##can now subset for just those lines with the starship counts; and then sum all
COG_prop5=COG_prop4 %>%
  filter(stringr::str_detect(dataset, 'starship'))
COG_prop5$k=sum(COG_prop5$count)


##now loop over each line and calculate the pvalue
##perform the hypergeometric COG for each COG using phyper(lower.tail = FALSE), and to COG both sides we COG q-1

##first create an empty dataframe
COG_prop6 = data.frame()
enrich_pval=""
derich_pval=""
##now loop through the lines and create an array of the resulting pvalues 
for(i in 1:nrow(COG_prop5)) {
  row <- COG_prop5[i,]
  # do stuff with row
  enrich_pval[[i]]=phyper(as.numeric(row$count)-1, as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = FALSE)
  derich_pval[[i]]=phyper(as.numeric(row$count), as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = TRUE)
}
##add the pvalues to the dataframe
COG_prop6 <- cbind(COG_prop5,enrich_pval)
COG_prop6 <- cbind(COG_prop6,derich_pval)
##now adjust the pvalues and then add that to another column
enrich_pval.adjust=p.adjust(enrich_pval, method = "BH")
derich_pval.adjust=p.adjust(derich_pval, method = "BH")

COG_prop6 <- cbind(COG_prop6,enrich_pval.adjust)
COG_prop6 <- cbind(COG_prop6,derich_pval.adjust)


##now add the adjusted pvalues, without differentiating between enriched and deriched
##we will not differentiate as the difference is clear and we are just showing statistically significant differences

COG_prop6$group1="core"
COG_prop6$group2="starship"

##add stars instead for the significance
COG_prop6$enrich_pval.adjust2 =stars.pval(COG_prop6$enrich_pval.adjust)
COG_prop6$derich_pval.adjust2 =stars.pval(COG_prop6$derich_pval.adjust)

ESE19=ggplot(data=COG_prop4 , aes(x=V1, y=proportion, fill=dataset))+
  geom_col(position= position_dodge())+
  theme_pubr()+
  xlab("COG category")+
  ylab("proportion")+theme(legend.title = element_blank())+
  labs(subtitle=expression(~italic(P.caseifulvum)~' ESE00019'))+
  scale_fill_manual(values = c("grey","red"))+
  stat_pvalue_manual(COG_prop6, label = "enrich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  stat_pvalue_manual(COG_prop6, label = "derich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  scale_x_discrete(labels=c("A"=expression(bold(A)), "C"=expression(bold(C)),"J"=expression(bold(J)),"M"=expression(bold(M)),"P"=expression(bold(P)),"U"=expression(bold(U)), parse=TRUE))


######
###now LCP05531
COG=read.csv(file="~/cluster2/projects/Penicillium/functional_enrichment/LCP05531.combined_mRNA.COG.tsv", header=F, sep ='\t')

##rewrite another column to siplify core and starship categories
COG$dataset<- ifelse(grepl("starship", COG$V3), "starship", "core")
##have a look at a stacked column
ggplot()+geom_col(data=COG , aes(x=dataset, y=V2, fill=V1), position= position_stack())+theme_pubr(x.text.angle = 35)
##calculate the proportion per group (starship, core)
COG_prop= COG %>%  
  group_by(dataset, V1) %>%
  summarise(count = sum(V2)) %>% 
  mutate(proportion = count / sum(count))
##check out a stacked column of the proportions now
ggplot()+
  geom_col(data=COG_prop , aes(x=dataset, y=proportion, fill=V1), position= position_stack())+
  theme_pubr(x.text.angle = 35)

###begin generating values to be used in the hypergeomtric tests
xx = COG_prop %>%  
  group_by(V1) %>% summarise(m = sum(count))
##merge with the proportion dataset
COG_prop=merge(COG_prop, xx, by='V1')

##add the total of all other genes (not the same COG) (n)
xx=COG_prop %>%  
  reframe(n = sum(count)-m)
COG_prop$n=xx

##add the total of all starship genes (k)
##can now subset for just those lines with the starship counts; and then sum all
COG_prop2=COG_prop %>%
  filter(stringr::str_detect(dataset, 'starship'))
COG_prop2$k=sum(COG_prop2$count)


##now loop over each line and calculate the pvalue
##perform the hypergeometric COG for each COG using phyper(lower.tail = FALSE), and to COG both sides we COG q-1

##first create an empty dataframe
COG_prop3 = data.frame()
enrich_pval=""
derich_pval=""
##now loop through the lines and create an array of the resulting pvalues 
for(i in 1:nrow(COG_prop2)) {
  row <- COG_prop2[i,]
  # do stuff with row
  enrich_pval[[i]]=phyper(as.numeric(row$count)-1, as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = FALSE)
  derich_pval[[i]]=phyper(as.numeric(row$count), as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = TRUE)
}
##add the pvalues to the dataframe
COG_prop3 <- cbind(COG_prop2,enrich_pval)
COG_prop3 <- cbind(COG_prop3,derich_pval)
##now adjust the pvalues and then add that to another column
enrich_pval.adjust=p.adjust(enrich_pval, method = "BH")
derich_pval.adjust=p.adjust(derich_pval, method = "BH")

COG_prop3 <- cbind(COG_prop3,enrich_pval.adjust)
COG_prop3 <- cbind(COG_prop3,derich_pval.adjust)

##now add the adjusted pvalues, without differentiating between enriched and deriched
##we will not differentiate as the difference is clear and we are just showing statistically significant differences

COG_prop3$group1="core"
COG_prop3$group2="starship"

##add stars instead for the significance
COG_prop3$enrich_pval.adjust2 =stars.pval(COG_prop3$enrich_pval.adjust)
COG_prop3$derich_pval.adjust2 =stars.pval(COG_prop3$derich_pval.adjust)

LCP05531unfilt=ggplot(data=COG_prop , aes(x=V1, y=proportion, fill=dataset))+
  geom_col(position= position_dodge())+
  theme_pubr()+
  xlab("COG category")+
  ylab("proportion")+theme(legend.title = element_blank())+
  labs(subtitle=expression(~italic(P.biforme)~' LCP05531'))+
  scale_fill_manual(values = c("grey","red"))+
  stat_pvalue_manual(COG_prop3, label = "enrich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  stat_pvalue_manual(COG_prop3, label = "derich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  scale_x_discrete(labels=c("A"=expression(bold(A)), "Absent"=expression(bold(Absent)), "C"=expression(bold(C)),"D"=expression(bold(D)),"G"=expression(bold(G)),"H"=expression(bold(H)),"J"=expression(bold(J)),"L"=expression(bold(L)),"O"=expression(bold(O)),"Q"=expression(bold(Q)),"T"=expression(bold(T)),"U"=expression(bold(U)),"Z"=expression(bold(Z)), parse=TRUE))





##calculate the proportion of each COG annotation with the core and starship genes seperately
##also re;ove any genes with COGs S (unknown), and L/K (transposable element associated). The latter is because the annotated regions were not softmasked before annotation due to size and databases
##also removed all the 'Absent' (genes with no COG annotation)
COG_prop4= subset(COG, V1 != "Absent" & V1 != "S" & V1 != "L" & V1 != "K") %>%  
  group_by(dataset, V1) %>%
  summarise(count = sum(V2)) %>% 
  mutate(proportion = count / sum(count))

##plot a stacked column with all the COG proportions for core vs starships
ggplot()+
  geom_col(data=COG_prop4 , aes(x=dataset, y=proportion, fill=V1), position= position_stack())+
  theme_pubr(x.text.angle = 35)

##can plot side by side proportions for each COG seperately
ggplot()+
  geom_col(data=COG_prop4 , aes(x=V1, y=proportion, fill=dataset), position= position_dodge())+
  theme_pubr()+
  xlab("COG category")+
  ylab("proportion")+theme(legend.title = element_blank())+
  scale_fill_manual(values = c("grey","red"))

##now need to calculate the hypergeometric COG p-value for all COG groups and then correct for all the COGs
##the four features needed for the hypergeometric COG are q, m, n, k = number genes in starships with the COG, the total number of genes (in and out of starships) with COG, the total count of all other genes (not with the same COG), the total of all starship genes 
##we will actually test for enrichment and deenrichment too

##get the sum of all genes with each COG (m)
xx = COG_prop4 %>%  
  group_by(V1) %>% summarise(m = sum(count))
##merge with the proportion dataset
COG_prop4=merge(COG_prop4, xx, by='V1')

##add the total of all other genes (not the same COG) (n)
xx=COG_prop4 %>%  
  reframe(n = sum(count)-m)
COG_prop4$n=xx

##add the total of all starship genes (k)
##can now subset for just those lines with the starship counts; and then sum all
COG_prop5=COG_prop4 %>%
  filter(stringr::str_detect(dataset, 'starship'))
COG_prop5$k=sum(COG_prop5$count)


##now loop over each line and calculate the pvalue
##perform the hypergeometric COG for each COG using phyper(lower.tail = FALSE), and to COG both sides we COG q-1

##first create an empty dataframe
COG_prop6 = data.frame()
enrich_pval=""
derich_pval=""
##now loop through the lines and create an array of the resulting pvalues 
for(i in 1:nrow(COG_prop5)) {
  row <- COG_prop5[i,]
  # do stuff with row
  enrich_pval[[i]]=phyper(as.numeric(row$count)-1, as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = FALSE)
  derich_pval[[i]]=phyper(as.numeric(row$count), as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = TRUE)
}
##add the pvalues to the dataframe
COG_prop6 <- cbind(COG_prop5,enrich_pval)
COG_prop6 <- cbind(COG_prop6,derich_pval)
##now adjust the pvalues and then add that to another column
enrich_pval.adjust=p.adjust(enrich_pval, method = "BH")
derich_pval.adjust=p.adjust(derich_pval, method = "BH")

COG_prop6 <- cbind(COG_prop6,enrich_pval.adjust)
COG_prop6 <- cbind(COG_prop6,derich_pval.adjust)


##now add the adjusted pvalues, without differentiating between enriched and deriched
##we will not differentiate as the difference is clear and we are just showing statistically significant differences

COG_prop6$group1="core"
COG_prop6$group2="starship"

##add stars instead for the significance
COG_prop6$enrich_pval.adjust2 =stars.pval(COG_prop6$enrich_pval.adjust)
COG_prop6$derich_pval.adjust2 =stars.pval(COG_prop6$derich_pval.adjust)

LCP05531=ggplot(data=COG_prop4 , aes(x=V1, y=proportion, fill=dataset))+
  geom_col(position= position_dodge())+
  theme_pubr()+
  xlab("COG category")+
  ylab("proportion")+theme(legend.title = element_blank())+
  labs(subtitle=expression(~italic(P.biforme)~' LCP05531'))+
  scale_fill_manual(values = c("grey","red"))+
  stat_pvalue_manual(COG_prop6, label = "enrich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  stat_pvalue_manual(COG_prop6, label = "derich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  scale_x_discrete(labels=c("A"=expression(bold(A)), "C"=expression(bold(C)),"J"=expression(bold(J)),"E"=expression(bold(E)),"P"=expression(bold(P)),"Z"=expression(bold(Z)), parse=TRUE))


######
###now LCP06093
COG=read.csv(file="~/cluster2/projects/Penicillium/functional_enrichment/LCP06093.combined_mRNA.COG.tsv", header=F, sep ='\t')

##rewrite another column to siplify core and starship categories
COG$dataset<- ifelse(grepl("starship", COG$V3), "starship", "core")
##have a look at a stacked column
ggplot()+geom_col(data=COG , aes(x=dataset, y=V2, fill=V1), position= position_stack())+theme_pubr(x.text.angle = 35)
##calculate the proportion per group (starship, core)
COG_prop= COG %>%  
  group_by(dataset, V1) %>%
  summarise(count = sum(V2)) %>% 
  mutate(proportion = count / sum(count))
##check out a stacked column of the proportions now
ggplot()+
  geom_col(data=COG_prop , aes(x=dataset, y=proportion, fill=V1), position= position_stack())+
  theme_pubr(x.text.angle = 35)

###begin generating values to be used in the hypergeomtric tests
xx = COG_prop %>%  
  group_by(V1) %>% summarise(m = sum(count))
##merge with the proportion dataset
COG_prop=merge(COG_prop, xx, by='V1')

##add the total of all other genes (not the same COG) (n)
xx=COG_prop %>%  
  reframe(n = sum(count)-m)
COG_prop$n=xx

##add the total of all starship genes (k)
##can now subset for just those lines with the starship counts; and then sum all
COG_prop2=COG_prop %>%
  filter(stringr::str_detect(dataset, 'starship'))
COG_prop2$k=sum(COG_prop2$count)


##now loop over each line and calculate the pvalue
##perform the hypergeometric COG for each COG using phyper(lower.tail = FALSE), and to COG both sides we COG q-1

##first create an empty dataframe
COG_prop3 = data.frame()
enrich_pval=""
derich_pval=""
##now loop through the lines and create an array of the resulting pvalues 
for(i in 1:nrow(COG_prop2)) {
  row <- COG_prop2[i,]
  # do stuff with row
  enrich_pval[[i]]=phyper(as.numeric(row$count)-1, as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = FALSE)
  derich_pval[[i]]=phyper(as.numeric(row$count), as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = TRUE)
}
##add the pvalues to the dataframe
COG_prop3 <- cbind(COG_prop2,enrich_pval)
COG_prop3 <- cbind(COG_prop3,derich_pval)
##now adjust the pvalues and then add that to another column
enrich_pval.adjust=p.adjust(enrich_pval, method = "BH")
derich_pval.adjust=p.adjust(derich_pval, method = "BH")

COG_prop3 <- cbind(COG_prop3,enrich_pval.adjust)
COG_prop3 <- cbind(COG_prop3,derich_pval.adjust)

##now add the adjusted pvalues, without differentiating between enriched and deriched
##we will not differentiate as the difference is clear and we are just showing statistically significant differences

COG_prop3$group1="core"
COG_prop3$group2="starship"

##add stars instead for the significance
COG_prop3$enrich_pval.adjust2 =stars.pval(COG_prop3$enrich_pval.adjust)
COG_prop3$derich_pval.adjust2 =stars.pval(COG_prop3$derich_pval.adjust)

LCP06093unfilt=ggplot(data=COG_prop , aes(x=V1, y=proportion, fill=dataset))+
  geom_col(position= position_dodge())+
  theme_pubr()+
  xlab("COG category")+
  ylab("proportion")+theme(legend.title = element_blank())+
  labs(subtitle=expression(~italic(P.camemberti)~' LCP06093'))+
  scale_fill_manual(values = c("grey","red"))+
  stat_pvalue_manual(COG_prop3, label = "enrich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  stat_pvalue_manual(COG_prop3, label = "derich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  scale_x_discrete(labels=c("A"=expression(bold(A)), "Absent"=expression(bold(Absent)), "C"=expression(bold(C)),"E"=expression(bold(E)),"G"=expression(bold(G)),"H"=expression(bold(H)),"J"=expression(bold(J)),"L"=expression(bold(L)),"O"=expression(bold(O)),"Q"=expression(bold(Q)),"T"=expression(bold(T)),"U"=expression(bold(U)), parse=TRUE))




##calculate the proportion of each COG annotation with the core and starship genes seperately
##also re;ove any genes with COGs S (unknown), and L/K (transposable element associated). The latter is because the annotated regions were not softmasked before annotation due to size and databases
##also removed all the 'Absent' (genes with no COG annotation)
COG_prop4= subset(COG, V1 != "Absent" & V1 != "S" & V1 != "L" & V1 != "K") %>%  
  group_by(dataset, V1) %>%
  summarise(count = sum(V2)) %>% 
  mutate(proportion = count / sum(count))

##plot a stacked column with all the COG proportions for core vs starships
ggplot()+
  geom_col(data=COG_prop4 , aes(x=dataset, y=proportion, fill=V1), position= position_stack())+
  theme_pubr(x.text.angle = 35)

##can plot side by side proportions for each COG seperately
ggplot()+
  geom_col(data=COG_prop4 , aes(x=V1, y=proportion, fill=dataset), position= position_dodge())+
  theme_pubr()+
  xlab("COG category")+
  ylab("proportion")+theme(legend.title = element_blank())+
  scale_fill_manual(values = c("grey","red"))

##now need to calculate the hypergeometric COG p-value for all COG groups and then correct for all the COGs
##the four features needed for the hypergeometric COG are q, m, n, k = number genes in starships with the COG, the total number of genes (in and out of starships) with COG, the total count of all other genes (not with the same COG), the total of all starship genes 
##we will actually test for enrichment and deenrichment too

##get the sum of all genes with each COG (m)
xx = COG_prop4 %>%  
  group_by(V1) %>% summarise(m = sum(count))
##merge with the proportion dataset
COG_prop4=merge(COG_prop4, xx, by='V1')

##add the total of all other genes (not the same COG) (n)
xx=COG_prop4 %>%  
  reframe(n = sum(count)-m)
COG_prop4$n=xx

##add the total of all starship genes (k)
##can now subset for just those lines with the starship counts; and then sum all
COG_prop5=COG_prop4 %>%
  filter(stringr::str_detect(dataset, 'starship'))
COG_prop5$k=sum(COG_prop5$count)


##now loop over each line and calculate the pvalue
##perform the hypergeometric COG for each COG using phyper(lower.tail = FALSE), and to COG both sides we COG q-1

##first create an empty dataframe
COG_prop6 = data.frame()
enrich_pval=""
derich_pval=""
##now loop through the lines and create an array of the resulting pvalues 
for(i in 1:nrow(COG_prop5)) {
  row <- COG_prop5[i,]
  # do stuff with row
  enrich_pval[[i]]=phyper(as.numeric(row$count)-1, as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = FALSE)
  derich_pval[[i]]=phyper(as.numeric(row$count), as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = TRUE)
}
##add the pvalues to the dataframe
COG_prop6 <- cbind(COG_prop5,enrich_pval)
COG_prop6 <- cbind(COG_prop6,derich_pval)
##now adjust the pvalues and then add that to another column
enrich_pval.adjust=p.adjust(enrich_pval, method = "BH")
derich_pval.adjust=p.adjust(derich_pval, method = "BH")

COG_prop6 <- cbind(COG_prop6,enrich_pval.adjust)
COG_prop6 <- cbind(COG_prop6,derich_pval.adjust)


##now add the adjusted pvalues, without differentiating between enriched and deriched
##we will not differentiate as the difference is clear and we are just showing statistically significant differences

COG_prop6$group1="core"
COG_prop6$group2="starship"

##add stars instead for the significance
COG_prop6$enrich_pval.adjust2 =stars.pval(COG_prop6$enrich_pval.adjust)
COG_prop6$derich_pval.adjust2 =stars.pval(COG_prop6$derich_pval.adjust)

LCP06093=ggplot(data=COG_prop4 , aes(x=V1, y=proportion, fill=dataset))+
  geom_col(position= position_dodge())+
  theme_pubr()+
  xlab("COG category")+
  ylab("proportion")+theme(legend.title = element_blank())+
  labs(subtitle=expression(~italic(P.camemberti)~' LCP06093'))+
  scale_fill_manual(values = c("grey","red"))+
  stat_pvalue_manual(COG_prop6, label = "enrich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  stat_pvalue_manual(COG_prop6, label = "derich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  scale_x_discrete(labels=c("A"=expression(bold(A)), "C"=expression(bold(C)),"H"=expression(bold(H)),"I"=expression(bold(I)),"J"=expression(bold(J)),"M"=expression(bold(M)),"P"=expression(bold(P)),"U"=expression(bold(U)), parse=TRUE))




#########
###now lastly ESE00090
COG=read.csv(file="~/cluster2/projects/Penicillium/functional_enrichment/ESE00090.combined_mRNA.COG.tsv", header=F, sep ='\t')

##rewrite another column to siplify core and starship categories
COG$dataset<- ifelse(grepl("starship", COG$V3), "starship", "core")
##have a look at a stacked column
ggplot()+geom_col(data=COG , aes(x=dataset, y=V2, fill=V1), position= position_stack())+theme_pubr(x.text.angle = 35)
##calculate the proportion per group (starship, core)
COG_prop= COG %>%  
  group_by(dataset, V1) %>%
  summarise(count = sum(V2)) %>% 
  mutate(proportion = count / sum(count))
##check out a stacked column of the proportions now
ggplot()+
  geom_col(data=COG_prop , aes(x=dataset, y=proportion, fill=V1), position= position_stack())+
  theme_pubr(x.text.angle = 35)

###begin generating values to be used in the hypergeomtric tests
xx = COG_prop %>%  
  group_by(V1) %>% summarise(m = sum(count))
##merge with the proportion dataset
COG_prop=merge(COG_prop, xx, by='V1')

##add the total of all other genes (not the same COG) (n)
xx=COG_prop %>%  
  reframe(n = sum(count)-m)
COG_prop$n=xx

##add the total of all starship genes (k)
##can now subset for just those lines with the starship counts; and then sum all
COG_prop2=COG_prop %>%
  filter(stringr::str_detect(dataset, 'starship'))
COG_prop2$k=sum(COG_prop2$count)


##now loop over each line and calculate the pvalue
##perform the hypergeometric COG for each COG using phyper(lower.tail = FALSE), and to COG both sides we COG q-1

##first create an empty dataframe
COG_prop3 = data.frame()
enrich_pval=""
derich_pval=""
##now loop through the lines and create an array of the resulting pvalues 
for(i in 1:nrow(COG_prop2)) {
  row <- COG_prop2[i,]
  # do stuff with row
  enrich_pval[[i]]=phyper(as.numeric(row$count)-1, as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = FALSE)
  derich_pval[[i]]=phyper(as.numeric(row$count), as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = TRUE)
}
##add the pvalues to the dataframe
COG_prop3 <- cbind(COG_prop2,enrich_pval)
COG_prop3 <- cbind(COG_prop3,derich_pval)
##now adjust the pvalues and then add that to another column
enrich_pval.adjust=p.adjust(enrich_pval, method = "BH")
derich_pval.adjust=p.adjust(derich_pval, method = "BH")

COG_prop3 <- cbind(COG_prop3,enrich_pval.adjust)
COG_prop3 <- cbind(COG_prop3,derich_pval.adjust)

##now add the adjusted pvalues, without differentiating between enriched and deriched
##we will not differentiate as the difference is clear and we are just showing statistically significant differences

COG_prop3$group1="core"
COG_prop3$group2="starship"

##add stars instead for the significance
COG_prop3$enrich_pval.adjust2 =stars.pval(COG_prop3$enrich_pval.adjust)
COG_prop3$derich_pval.adjust2 =stars.pval(COG_prop3$derich_pval.adjust)

ESE90unfilt=ggplot(data=COG_prop , aes(x=V1, y=proportion, fill=dataset))+
  geom_col(position= position_dodge())+
  theme_pubr()+
  xlab("COG category")+
  ylab("proportion")+theme(legend.title = element_blank())+
  labs(subtitle=expression(~italic(P.fuscoglaucum)~' ESE00090'))+
  scale_fill_manual(values = c("grey","red"))+
  stat_pvalue_manual(COG_prop3, label = "enrich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  stat_pvalue_manual(COG_prop3, label = "derich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  scale_x_discrete(labels=c("Absent"=expression(bold(Absent)), parse=TRUE))


##calculate the proportion of each COG annotation with the core and starship genes seperately
##also re;ove any genes with COGs S (unknown), and L/K (transposable element associated). The latter is because the annotated regions were not softmasked before annotation due to size and databases
##also removed all the 'Absent' (genes with no COG annotation)
COG_prop4= subset(COG, V1 != "Absent" & V1 != "S" & V1 != "L" & V1 != "K") %>%  
  group_by(dataset, V1) %>%
  summarise(count = sum(V2)) %>% 
  mutate(proportion = count / sum(count))

##plot a stacked column with all the COG proportions for core vs starships
ggplot()+
  geom_col(data=COG_prop4 , aes(x=dataset, y=proportion, fill=V1), position= position_stack())+
  theme_pubr(x.text.angle = 35)

##can plot side by side proportions for each COG seperately
ggplot()+
  geom_col(data=COG_prop4 , aes(x=V1, y=proportion, fill=dataset), position= position_dodge())+
  theme_pubr()+
  xlab("COG category")+
  ylab("proportion")+theme(legend.title = element_blank())+
  scale_fill_manual(values = c("grey","red"))

##now need to calculate the hypergeometric COG p-value for all COG groups and then correct for all the COGs
##the four features needed for the hypergeometric COG are q, m, n, k = number genes in starships with the COG, the total number of genes (in and out of starships) with COG, the total count of all other genes (not with the same COG), the total of all starship genes 
##we will actually test for enrichment and deenrichment too

##get the sum of all genes with each COG (m)
xx = COG_prop4 %>%  
  group_by(V1) %>% summarise(m = sum(count))
##merge with the proportion dataset
COG_prop4=merge(COG_prop4, xx, by='V1')

##add the total of all other genes (not the same COG) (n)
xx=COG_prop4 %>%  
  reframe(n = sum(count)-m)
COG_prop4$n=xx

##add the total of all starship genes (k)
##can now subset for just those lines with the starship counts; and then sum all
COG_prop5=COG_prop4 %>%
  filter(stringr::str_detect(dataset, 'starship'))
COG_prop5$k=sum(COG_prop5$count)


##now loop over each line and calculate the pvalue
##perform the hypergeometric COG for each COG using phyper(lower.tail = FALSE), and to COG both sides we COG q-1

##first create an empty dataframe
COG_prop6 = data.frame()
enrich_pval=""
derich_pval=""
##now loop through the lines and create an array of the resulting pvalues 
for(i in 1:nrow(COG_prop5)) {
  row <- COG_prop5[i,]
  # do stuff with row
  enrich_pval[[i]]=phyper(as.numeric(row$count)-1, as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = FALSE)
  derich_pval[[i]]=phyper(as.numeric(row$count), as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = TRUE)
}
##add the pvalues to the dataframe
COG_prop6 <- cbind(COG_prop5,enrich_pval)
COG_prop6 <- cbind(COG_prop6,derich_pval)
##now adjust the pvalues and then add that to another column
enrich_pval.adjust=p.adjust(enrich_pval, method = "BH")
derich_pval.adjust=p.adjust(derich_pval, method = "BH")

COG_prop6 <- cbind(COG_prop6,enrich_pval.adjust)
COG_prop6 <- cbind(COG_prop6,derich_pval.adjust)


##now add the adjusted pvalues, without differentiating between enriched and deriched
##we will not differentiate as the difference is clear and we are just showing statistically significant differences

COG_prop6$group1="core"
COG_prop6$group2="starship"

##add stars instead for the significance
COG_prop6$enrich_pval.adjust2 =stars.pval(COG_prop6$enrich_pval.adjust)
COG_prop6$derich_pval.adjust2 =stars.pval(COG_prop6$derich_pval.adjust)

ESE90=ggplot(data=COG_prop4 , aes(x=V1, y=proportion, fill=dataset))+
  geom_col(position= position_dodge())+
  theme_pubr()+
  xlab("COG category")+
  ylab("proportion")+
  theme(legend.title = element_blank())+
  labs(subtitle=expression(~italic(P.fuscoglaucum)~' ESE00090'))+
  scale_fill_manual(values = c("grey","red"))+
  stat_pvalue_manual(COG_prop6, label = "enrich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  stat_pvalue_manual(COG_prop6, label = "derich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)



###now combine them all together
ggarrange(ESE19unfilt, LCP05531unfilt, LCP06093unfilt, ESE90unfilt, ncol = 1, common.legend = T)
##saved as 1000x1000 svg 'core_vs_starship.unfiltered.combined'
ggarrange(ESE19, LCP05531, LCP06093, ESE90, ncol = 1, common.legend = T)
##saved as 1000x1000 svg 'core_vs_starship.combined'
combined=ggarrange(ESE19, LCP05531, LCP06093, ESE90, ncol = 1, common.legend = T)

##try to add the meanings of the COG in a table
##only for those with significant differences
COGkey=read.csv(file="~/cluster2/projects/Penicillium/functional_enrichment/COGs_key_meaning.tsv", sep='\t', header=T)

COGkey2 <- cbind(COGkey)
ggplot()+theme_void()+ annotation_custom(tableGrob(COGkey, rows=NULL))
table=ggplot()+theme_void() +annotation_custom(tableGrob(COGkey, rows=NULL))

ggarrange(combined, table, ncol = 2, widths = c(2,1))


##can also combine all the cheese strains as they are very similar to look at general trends
COG=read.csv(file="~/cluster2/projects/Penicillium/functional_enrichment/all_strains.combined_mRNA.COG.tsv", header=F, sep ='\t')

##rewrite another column to siplify core and starship categories
COG$dataset<- ifelse(grepl("starship", COG$V3), "starship", "core")
##have a look at a stacked column
ggplot()+geom_col(data=COG , aes(x=dataset, y=V2, fill=V1), position= position_stack())+theme_pubr(x.text.angle = 35)
##calculate the proportion per group (starship, core)
COG_prop= COG %>%  
  group_by(dataset, V1) %>%
  summarise(count = sum(V2)) %>% 
  mutate(proportion = count / sum(count))
##check out a stacked column of the proportions now
ggplot()+
  geom_col(data=COG_prop , aes(x=dataset, y=proportion, fill=V1), position= position_stack())+
  theme_pubr(x.text.angle = 35)


##calculate the proportion of each COG annotation with the core and starship genes seperately
##also re;ove any genes with COGs S (unknown), and L/K (transposable element associated). The latter is because the annotated regions were not softmasked before annotation due to size and databases
##also removed all the 'Absent' (genes with no COG annotation)
COG_prop4= subset(COG, V1 != "Absent" & V1 != "S" & V1 != "L" & V1 != "K") %>%  
  group_by(dataset, V1) %>%
  summarise(count = sum(V2)) %>% 
  mutate(proportion = count / sum(count))

##plot a stacked column with all the COG proportions for core vs starships
ggplot()+
  geom_col(data=COG_prop4 , aes(x=dataset, y=proportion, fill=V1), position= position_stack())+
  theme_pubr(x.text.angle = 35)

##can plot side by side proportions for each COG seperately
ggplot()+
  geom_col(data=COG_prop4 , aes(x=V1, y=proportion, fill=dataset), position= position_dodge())+
  theme_pubr()+
  xlab("COG category")+
  ylab("proportion")+theme(legend.title = element_blank())+
  scale_fill_manual(values = c("grey","red"))

##now need to calculate the hypergeometric COG p-value for all COG groups and then correct for all the COGs
##the four features needed for the hypergeometric COG are q, m, n, k = number genes in starships with the COG, the total number of genes (in and out of starships) with COG, the total count of all other genes (not with the same COG), the total of all starship genes 
##we will actually test for enrichment and deenrichment too

##get the sum of all genes with each COG (m)
xx = COG_prop4 %>%  
  group_by(V1) %>% summarise(m = sum(count))
##merge with the proportion dataset
COG_prop4=merge(COG_prop4, xx, by='V1')

##add the total of all other genes (not the same COG) (n)
xx=COG_prop4 %>%  
  reframe(n = sum(count)-m)
COG_prop4$n=xx

##add the total of all starship genes (k)
##can now subset for just those lines with the starship counts; and then sum all
COG_prop5=COG_prop4 %>%
  filter(stringr::str_detect(dataset, 'starship'))
COG_prop5$k=sum(COG_prop5$count)


##now loop over each line and calculate the pvalue
##perform the hypergeometric COG for each COG using phyper(lower.tail = FALSE), and to COG both sides we COG q-1

##first create an empty dataframe
COG_prop6 = data.frame()
enrich_pval=""
derich_pval=""
##now loop through the lines and create an array of the resulting pvalues 
for(i in 1:nrow(COG_prop5)) {
  row <- COG_prop5[i,]
  # do stuff with row
  enrich_pval[[i]]=phyper(as.numeric(row$count)-1, as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = FALSE)
  derich_pval[[i]]=phyper(as.numeric(row$count), as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = TRUE)
}
##add the pvalues to the dataframe
COG_prop6 <- cbind(COG_prop5,enrich_pval)
COG_prop6 <- cbind(COG_prop6,derich_pval)
##now adjust the pvalues and then add that to another column
enrich_pval.adjust=p.adjust(enrich_pval, method = "BH")
derich_pval.adjust=p.adjust(derich_pval, method = "BH")

COG_prop6 <- cbind(COG_prop6,enrich_pval.adjust)
COG_prop6 <- cbind(COG_prop6,derich_pval.adjust)


##now add the adjusted pvalues, without differentiating between enriched and deriched
##we will not differentiate as the difference is clear and we are just showing statistically significant differences

COG_prop6$group1="core"
COG_prop6$group2="starship"

##add stars instead for the significance
COG_prop6$enrich_pval.adjust2 =stars.pval(COG_prop6$enrich_pval.adjust)
COG_prop6$derich_pval.adjust2 =stars.pval(COG_prop6$derich_pval.adjust)

ggplot(data=COG_prop4 , aes(x=V1, y=proportion, fill=dataset))+
  geom_col(position= position_dodge())+
  theme_pubr()+
  xlab("COG category")+
  ylab("proportion")+
  theme(legend.title = element_blank())+
  labs(subtitle=expression(~italic(P.biforme)~' LCP05531; '~italic(P.caseifulvum)~' ESE00019; '~italic(P.camemberti)~' LCP06093'))+
  scale_fill_manual(values = c("grey","red"))+
  stat_pvalue_manual(COG_prop6, label = "enrich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  stat_pvalue_manual(COG_prop6, label = "derich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  scale_x_discrete(labels=c("A"=expression(bold(A)), "C"=expression(bold(C)),"E"=expression(bold(E)),"H"=expression(bold(H)),"I"=expression(bold(I)),"J"=expression(bold(J)),"M"=expression(bold(M)),"P"=expression(bold(P)),"U"=expression(bold(U)),"Y"=expression(bold(Y)),"Z"=expression(bold(Z)), parse=TRUE))

cheesefunc=ggplot(data=COG_prop4 , aes(x=V1, y=proportion, fill=dataset))+
  geom_col(position= position_dodge())+
  theme_pubr()+
  xlab("COG category")+
  ylab("proportion")+
  theme(legend.title = element_blank())+
  labs(subtitle=expression(~italic(P.biforme)~' LCP05531; '~italic(P.caseifulvum)~' ESE00019; '~italic(P.camemberti)~' LCP06093'))+
  scale_fill_manual(values = c("grey","red"))+
  stat_pvalue_manual(COG_prop6, label = "enrich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  stat_pvalue_manual(COG_prop6, label = "derich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  scale_x_discrete(labels=c("A"=expression(bold(A)), "C"=expression(bold(C)),"E"=expression(bold(E)),"H"=expression(bold(H)),"I"=expression(bold(I)),"J"=expression(bold(J)),"M"=expression(bold(M)),"P"=expression(bold(P)),"U"=expression(bold(U)),"Y"=expression(bold(Y)),"Z"=expression(bold(Z)), parse=TRUE))




##########
###we can also compare the COGs between starships found in environmetnal, food-related and clinical strains


##read in data
##data contains a factor for starship or core, the COG and the count of COG for each group
COGisol=read.csv(file="~/cluster2/projects/Penicillium/functional_enrichment/isolation_comparison/all_isolations.COGs.tsv", header=F, sep ='\t')

##rewrite another column to siplify core and starship categories
COGisol$dataset = factor(COGisol$V3, levels = c("environment", "clinical", "food"))
##have a look at a stacked column
ggplot()+geom_col(data=COGisol , aes(x=dataset, y=V2, fill=V1), position= position_stack())+theme_pubr(x.text.angle = 35)
## do pairwise comparisons of datasets
COGfoodenv=subset(COGisol, dataset == "food" | dataset == "environment")
##calculate the proportion per isolation category
COGfoodenv_prop= COGfoodenv %>%  
  group_by(dataset, V1) %>%
  summarise(count = sum(V2)) %>% 
  mutate(proportion = count / sum(count))

##check out a stacked column of the proportions now
ggplot()+
  geom_col(data=COGfoodenv_prop , aes(x=dataset, y=proportion, fill=V1), position= position_stack())+
  theme_pubr(x.text.angle = 35)


##calculate the proportion of each COG annotation with the core and starship genes seperately
##also re;ove any genes with COG S (unknown), and L/K (transposable element associated). The latter is because the annotated regions were not softmasked before annotation due to size and databases
##also removed all the 'Absent' (genes with no annotation)
COGfoodenv_prop4= subset(COGfoodenv, V1 != "Absent" & V1 != "S" & V1 != "L" & V1 != "K") %>%  
  group_by(dataset, V1) %>%
  summarise(count = sum(V2)) %>% 
  mutate(proportion = count / sum(count))

##plot a stacked column with all the COG proportions for core vs starships
ggplot()+
  geom_col(data=COGfoodenv_prop4 , aes(x=dataset, y=proportion, fill=V1), position= position_stack())+
  theme_pubr(x.text.angle = 35)

##can plot side by side proportions for each COGisol seperately
ggplot()+
  geom_col(data=COGfoodenv_prop4 , aes(x=V1, y=proportion, fill=dataset), position= position_dodge())+
  theme_pubr()+
  xlab("COGisol category")+
  ylab("proportion")+theme(legend.title = element_blank())+
  scale_fill_manual(values = c("#009E73", "#0072B2"))

##now need to calculate the hypergeometric COGisol p-value for all COGisol groups and then correct for all the COGisols
##the four features needed for the hypergeometric COGisol are q, m, n, k = number genes in starships with the COGisol, the total number of genes (in and out of starships) with COGisol, the total count of all other genes (not with the same COGisol), the total of all starship genes 
##we will actually test for enrichment and deenrichment too

##get the sum of all genes with each COGisol (m)
xx = COGfoodenv_prop4 %>%  
  group_by(V1) %>% summarise(m = sum(count))
##merge with the proportion dataset
COGfoodenv_prop4=merge(COGfoodenv_prop4, xx, by='V1')

##add the total of all other genes (not the same COGisol) (n)
xx=COGfoodenv_prop4 %>%  
  reframe(n = sum(count)-m)
COGfoodenv_prop4$n=xx

##add the total of all food genes (k)
##can now subset for just those lines with the food counts; and then sum all
COGfoodenv_prop5=COGfoodenv_prop4 %>%
  filter(stringr::str_detect(dataset, 'food'))
COGfoodenv_prop5$k=sum(COGfoodenv_prop5$count)


##now loop over each line and calculate the pvalue
##perform the hypergeometric COGisol for each COGisol using phyper(lower.tail = FALSE), and to COGisol both sides we COGisol q-1

##first create an empty dataframe
COGfoodenv_prop6 = data.frame()
enrich_pval=""
derich_pval=""
##now loop through the lines and create an array of the resulting pvalues 
for(i in 1:nrow(COGfoodenv_prop5)) {
  row <- COGfoodenv_prop5[i,]
  # do stuff with row
  enrich_pval[[i]]=phyper(as.numeric(row$count)-1, as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = FALSE)
  derich_pval[[i]]=phyper(as.numeric(row$count), as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = TRUE)
}
##add the pvalues to the dataframe
COGfoodenv_prop6 <- cbind(COGfoodenv_prop5,enrich_pval)
COGfoodenv_prop6 <- cbind(COGfoodenv_prop6,derich_pval)
##now adjust the pvalues and then add that to another column
enrich_pval.adjust=p.adjust(enrich_pval, method = "BH")
derich_pval.adjust=p.adjust(derich_pval, method = "BH")

COGfoodenv_prop6 <- cbind(COGfoodenv_prop6,enrich_pval.adjust)
COGfoodenv_prop6 <- cbind(COGfoodenv_prop6,derich_pval.adjust)


##now add the adjusted pvalues, without differentiating between enriched and deriched
##we will not differentiate as the difference is clear and we are just showing statistically significant differences

COGfoodenv_prop6$group1="environment"
COGfoodenv_prop6$group2="food"

##add stars instead for the significance
COGfoodenv_prop6$enrich_pval.adjust2 =stars.pval(COGfoodenv_prop6$enrich_pval.adjust)
COGfoodenv_prop6$derich_pval.adjust2 =stars.pval(COGfoodenv_prop6$derich_pval.adjust)

COGfoodenv_prop_plot=ggplot(data=COGfoodenv_prop4 , aes(x=V1, y=proportion, fill=dataset))+
  geom_col(position= position_dodge())+
  theme_pubr()+
  xlab("COG category")+
  ylab("proportion")+
  theme(legend.title = element_blank())+
  scale_fill_manual(values = c("#009E73", "#0072B2"))+
  stat_pvalue_manual(COGfoodenv_prop6, label = "enrich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  stat_pvalue_manual(COGfoodenv_prop6, label = "derich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  scale_x_discrete(labels=c("B"=expression(bold(B)), "C"=expression(bold(C)),"D"=expression(bold(D)),"F"=expression(bold(F)),"G"=expression(bold(G)),"I"=expression(bold(I)),"J"=expression(bold(J)),"O"=expression(bold(O)),"M"=expression(bold(M)),"P"=expression(bold(P)),"Q"=expression(bold(Q)),"V"=expression(bold(V)), parse=TRUE))




##########
###now the same but with clinical vs env
COGclinenv=subset(COGisol, dataset == "clinical" | dataset == "environment")
##calculate the proportion per isolation category
COGclinenv_prop= COGclinenv %>%  
  group_by(dataset, V1) %>%
  summarise(count = sum(V2)) %>% 
  mutate(proportion = count / sum(count))

##check out a stacked column of the proportions now
ggplot()+
  geom_col(data=COGclinenv_prop , aes(x=dataset, y=proportion, fill=V1), position= position_stack())+
  theme_pubr(x.text.angle = 35)


##calculate the proportion of each COG annotation with the core and starship genes seperately
##also re;ove any genes with COG S (unknown), and L/K (transposable element associated). The latter is because the annotated regions were not softmasked before annotation due to size and databases
##also removed all the 'Absent' (genes with no annotation)
COGclinenv_prop4= subset(COGclinenv, V1 != "Absent" & V1 != "S" & V1 != "L" & V1 != "K") %>%  
  group_by(dataset, V1) %>%
  summarise(count = sum(V2)) %>% 
  mutate(proportion = count / sum(count))

##plot a stacked column with all the COG proportions for core vs starships
ggplot()+
  geom_col(data=COGclinenv_prop4 , aes(x=dataset, y=proportion, fill=V1), position= position_stack())+
  theme_pubr(x.text.angle = 35)

##can plot side by side proportions for each COGisol seperately
ggplot()+
  geom_col(data=COGclinenv_prop4 , aes(x=V1, y=proportion, fill=dataset), position= position_dodge())+
  theme_pubr()+
  xlab("COGisol category")+
  ylab("proportion")+theme(legend.title = element_blank())+
  scale_fill_manual(values = c("#009E73", "#E69F00"))

##now need to calculate the hypergeometric COGisol p-value for all COGisol groups and then correct for all the COGisols
##the four features needed for the hypergeometric COGisol are q, m, n, k = number genes in starships with the COGisol, the total number of genes (in and out of starships) with COGisol, the total count of all other genes (not with the same COGisol), the total of all starship genes 
##we will actually test for enrichment and deenrichment too

##get the sum of all genes with each COGisol (m)
xx = COGclinenv_prop4 %>%  
  group_by(V1) %>% summarise(m = sum(count))
##merge with the proportion dataset
COGclinenv_prop4=merge(COGclinenv_prop4, xx, by='V1')

##add the total of all other genes (not the same COGisol) (n)
xx=COGclinenv_prop4 %>%  
  reframe(n = sum(count)-m)
COGclinenv_prop4$n=xx

##add the total of all food genes (k)
##can now subset for just those lines with the food counts; and then sum all
COGclinenv_prop5=COGclinenv_prop4 %>%
  filter(stringr::str_detect(dataset, 'clinical'))
COGclinenv_prop5$k=sum(COGclinenv_prop5$count)


##now loop over each line and calculate the pvalue
##perform the hypergeometric COGisol for each COGisol using phyper(lower.tail = FALSE), and to COGisol both sides we COGisol q-1

##first create an empty dataframe
COGclinenv_prop6 = data.frame()
enrich_pval=""
derich_pval=""
##now loop through the lines and create an array of the resulting pvalues 
for(i in 1:nrow(COGclinenv_prop5)) {
  row <- COGclinenv_prop5[i,]
  # do stuff with row
  enrich_pval[[i]]=phyper(as.numeric(row$count)-1, as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = FALSE)
  derich_pval[[i]]=phyper(as.numeric(row$count), as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = TRUE)
}
##add the pvalues to the dataframe
COGclinenv_prop6 <- cbind(COGclinenv_prop5,enrich_pval)
COGclinenv_prop6 <- cbind(COGclinenv_prop6,derich_pval)
##now adjust the pvalues and then add that to another column
enrich_pval.adjust=p.adjust(enrich_pval, method = "BH")
derich_pval.adjust=p.adjust(derich_pval, method = "BH")

COGclinenv_prop6 <- cbind(COGclinenv_prop6,enrich_pval.adjust)
COGclinenv_prop6 <- cbind(COGclinenv_prop6,derich_pval.adjust)


##now add the adjusted pvalues, without differentiating between enriched and deriched
##we will not differentiate as the difference is clear and we are just showing statistically significant differences

COGclinenv_prop6$group1="environment"
COGclinenv_prop6$group2="clinical"

##add stars instead for the significance
COGclinenv_prop6$enrich_pval.adjust2 =stars.pval(COGclinenv_prop6$enrich_pval.adjust)
COGclinenv_prop6$derich_pval.adjust2 =stars.pval(COGclinenv_prop6$derich_pval.adjust)

COGclinenv_prop_plot=ggplot(data=COGclinenv_prop4 , aes(x=V1, y=proportion, fill=dataset))+
  geom_col(position= position_dodge())+
  theme_pubr()+
  xlab("COG category")+
  ylab("proportion")+
  theme(legend.title = element_blank())+
  scale_fill_manual(values = c("#009E73", "#E69F00"))+
  stat_pvalue_manual(COGclinenv_prop6, label = "enrich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  stat_pvalue_manual(COGclinenv_prop6, label = "derich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)



##########
###once more but clinical vs food
COGclinfood=subset(COGisol, dataset == "clinical" | dataset == "food")
##calculate the proportion per isolation category
COGclinfood_prop= COGclinfood %>%  
  group_by(dataset, V1) %>%
  summarise(count = sum(V2)) %>% 
  mutate(proportion = count / sum(count))

##check out a stacked column of the proportions now
ggplot()+
  geom_col(data=COGclinfood_prop , aes(x=dataset, y=proportion, fill=V1), position= position_stack())+
  theme_pubr(x.text.angle = 35)


##calculate the proportion of each COG annotation with the core and starship genes seperately
##also re;ove any genes with COG S (unknown), and L/K (transposable element associated). The latter is because the annotated regions were not softmasked before annotation due to size and databases
##also removed all the 'Absent' (genes with no annotation)
COGclinfood_prop4= subset(COGclinfood, V1 != "Absent" & V1 != "S" & V1 != "L" & V1 != "K") %>%  
  group_by(dataset, V1) %>%
  summarise(count = sum(V2)) %>% 
  mutate(proportion = count / sum(count))

##plot a stacked column with all the COG proportions for core vs starships
ggplot()+
  geom_col(data=COGclinfood_prop4 , aes(x=dataset, y=proportion, fill=V1), position= position_stack())+
  theme_pubr(x.text.angle = 35)

##can plot side by side proportions for each COGisol seperately
ggplot()+
  geom_col(data=COGclinfood_prop4 , aes(x=V1, y=proportion, fill=dataset), position= position_dodge())+
  theme_pubr()+
  xlab("COGisol category")+
  ylab("proportion")+theme(legend.title = element_blank())+
  scale_fill_manual(values = c("#E69F00", "#0072B2"))

##now need to calculate the hypergeometric COGisol p-value for all COGisol groups and then correct for all the COGisols
##the four features needed for the hypergeometric COGisol are q, m, n, k = number genes in starships with the COGisol, the total number of genes (in and out of starships) with COGisol, the total count of all other genes (not with the same COGisol), the total of all starship genes 
##we will actually test for enrichment and deenrichment too

##get the sum of all genes with each COGisol (m)
xx = COGclinfood_prop4 %>%  
  group_by(V1) %>% summarise(m = sum(count))
##merge with the proportion dataset
COGclinfood_prop4=merge(COGclinfood_prop4, xx, by='V1')

##add the total of all other genes (not the same COGisol) (n)
xx=COGclinfood_prop4 %>%  
  reframe(n = sum(count)-m)
COGclinfood_prop4$n=xx

##add the total of all food genes (k)
##can now subset for just those lines with the food counts; and then sum all
COGclinfood_prop5=COGclinfood_prop4 %>%
  filter(stringr::str_detect(dataset, 'clinical'))
COGclinfood_prop5$k=sum(COGclinfood_prop5$count)


##now loop over each line and calculate the pvalue
##perform the hypergeometric COGisol for each COGisol using phyper(lower.tail = FALSE), and to COGisol both sides we COGisol q-1

##first create an empty dataframe
COGclinfood_prop6 = data.frame()
enrich_pval=""
derich_pval=""
##now loop through the lines and create an array of the resulting pvalues 
for(i in 1:nrow(COGclinfood_prop5)) {
  row <- COGclinfood_prop5[i,]
  # do stuff with row
  enrich_pval[[i]]=phyper(as.numeric(row$count)-1, as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = FALSE)
  derich_pval[[i]]=phyper(as.numeric(row$count), as.numeric(row$m), as.numeric(row$n), as.numeric(row$k), lower.tail = TRUE)
}
##add the pvalues to the dataframe
COGclinfood_prop6 <- cbind(COGclinfood_prop5,enrich_pval)
COGclinfood_prop6 <- cbind(COGclinfood_prop6,derich_pval)
##now adjust the pvalues and then add that to another column
enrich_pval.adjust=p.adjust(enrich_pval, method = "BH")
derich_pval.adjust=p.adjust(derich_pval, method = "BH")

COGclinfood_prop6 <- cbind(COGclinfood_prop6,enrich_pval.adjust)
COGclinfood_prop6 <- cbind(COGclinfood_prop6,derich_pval.adjust)


##now add the adjusted pvalues, without differentiating between enriched and deriched
##we will not differentiate as the difference is clear and we are just showing statistically significant differences

COGclinfood_prop6$group1="food"
COGclinfood_prop6$group2="clinical"

##add stars instead for the significance
COGclinfood_prop6$enrich_pval.adjust2 =stars.pval(COGclinfood_prop6$enrich_pval.adjust)
COGclinfood_prop6$derich_pval.adjust2 =stars.pval(COGclinfood_prop6$derich_pval.adjust)

COGclinfood_prop_plot=ggplot(data=COGclinfood_prop4 , aes(x=V1, y=proportion, fill=dataset))+
  geom_col(position= position_dodge())+
  theme_pubr()+
  xlab("COG category")+
  ylab("proportion")+
  theme(legend.title = element_blank())+
  scale_fill_manual(values = c("#E69F00", "#0072B2"))+
  stat_pvalue_manual(COGclinfood_prop6, label = "enrich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)+
  stat_pvalue_manual(COGclinfood_prop6, label = "derich_pval.adjust2", y.position = "proportion", x = "V1", size = 5)


ggarrange(COGfoodenv_prop_plot, COGclinenv_prop_plot)
temp1=ggarrange(COGfoodenv_prop_plot, COGclinenv_prop_plot)



##combine the core-vs-starship and starship isolate comparisons
ggarrange(cheesefunc, temp1, ncol=1)
##saved as svg 1000x500 'functionalenrichment_combined.svg'
ggarrange(combined, COGfoodenv_prop_plot, COGclinenv_prop_plot, ncol=1, heights = c(4,1,1))
##saved as svg 700x1000 'functionalenrichment.combined_full.svg'
