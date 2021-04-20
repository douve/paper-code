








# Exces of mortality -----------------------------------------------------------

## exitus pob vs residents
vlines = data.frame(vals = c(min(resi.anl$Date_first_test,na.rm=T),
                             median(resi.anl$Date_first_test,na.rm=T)+4,
                             median(resi.anl$Date_first_test,na.rm=T)+18),
                    label = c('First PCR','Date of intervention','14 days post intervention'))
vlines$dif = c(0,diff(vlines$vals)) %>% as.numeric %>% cumsum
vlines$week = lubridate::week(vlines$vals)

xitus_long$week = exitus_long$Data_exitus %>% lubridate::week()

data = exitus_long %>% filter(Data_exitus>=start.date) %>% filter(Data_exitus<=end.date) %>%
  group_by(week) %>% summarise(pob =sum(Poblacio),
                               resis = sum(Residents),
                               covid19 = sum(Resi_covid)) %>%
  mutate(nocovid_residents = resis-covid19,
         nocovid_pob = pob - resis) %>%
  left_join(exitus_old %>% group_by(year,week) %>% 
              summarise(pob19 = n()) %>% data.frame() %>%
              group_by(week) %>%
              summarise(pob19=median(pob19))) %>%
  mutate(excess = ifelse(pob19>resis,0,resis-pob19),
         excess = ifelse(covid19>excess,0,excess-covid19)
  ) %>% data.frame

bpdat = data %>% select(week,pob,pob19,excess,covid19) %>%
  reshape2::melt(id.vars='week') %>% 
  mutate(variable=factor(variable,levels=rev(c('pob19','excess','covid19')),
                         labels=rev(c('Median exitus 2016-2019','Exitus 2020','Exitus Covid 2020'))))

bpdat = left_join(bpdat,data%>% select(week,pob)) %>% mutate(pob=pob/1.5)
bpdat$date = lubridate::ymd( "2020-01-01" ) + lubridate::weeks( bpdat$week - 1 )

bp = ggplot(bpdat %>% filter(variable!='Poblacio') %>% mutate(variable=droplevels(variable)),
) + 
  geom_line(mapping= aes(x=date,y= pob),size=0.5,col='gray30')+
  geom_point(mapping= aes(x=date,y= pob),size=2,col='black',alpha=0.5)+
  geom_bar(mapping=aes(x=date,y=value,fill=variable),stat = 'identity',col='black',size=0.1) + 
  scale_fill_manual(values = c('brown3','salmon','gray65'))+
  scale_y_continuous(name = "Excess mortality in LTFC",
                     sec.axis = sec_axis(~.*1.5, name = "Mortality in general population"))+
  scale_x_date(date_breaks ='2 weeks',date_labels = '%b %d' )+
  geom_vline(xintercept=vlines$vals[1],size=0.5,linetype='dotted')+
  annotate("text",x=vlines$vals[1],y=535,label=vlines$label[1],
           hjust=-0.1,vjust=2,size=3)+
  theme_minimal()+
  theme(axis.ticks.x=element_blank(),
        legend.title=element_blank(),
        panel.grid.major.x =element_blank(),
        panel.grid.minor.x =element_blank())









# Heatmap -------------------------------------------------------
data = resi.anl %>% 
  select(Codi,Desc_residencia,
         louvain = louvain2,
         median_Age,
         Perc_Dones,
         median_comorbi,
         Residents_alfa,
         Perc_PCC,
         Perc_MACA,
         N_residents_actius,
         SNQ12,
         Perc_leave,
         renta_llar,
         N_pob,
         incidencia
         ) %>% na.omit 

# scale
data = data %>% select(-(1:3)) %>% 
  mutate_if(is.Date,as.numeric) %>% 
  mutate_if(is.factor,as.numeric) %>% 
  scale %>% data.frame %>%
  cbind.data.frame(data %>% select(1:3),.)
params = names(data)[-c(1:3)]
data$Residents_alfa = -data$Residents_alfa
data$renta_llar = -data$renta_llar

classic.colors = colorRampPalette(c(rgb(75,145,35,max=255),
                                    rgb(75,145,35,max=255),
                                    'white',
                                    rgb(200,30,125,max=255),
                                    rgb(200,30,125,max=255)
                                    ))

aux = data %>%
  dplyr::select(louvain,all_of(params)) %>%
  na.omit() %>% group_by(louvain)
aux = aux %>% dplyr::summarise_all(median) %>%
  select(-louvain) %>% t() %>% data.frame
colnames(aux) = table(data[,'louvain']) %>% names
aux = t(aux)

# Calculate cluster frequencies
clustering_table <- as.numeric(table(data %>% select(louvain)))
labels_row <- paste0(rownames(aux), " (",formatC(clustering_table,format="g"),', ',
                     round(clustering_table / sum(clustering_table) * 100, 2), "%)")
labels_col = c(
  'Elderliness',
  '% of males',
  'Number of comorbidities',
  '% of Barthel score <50 recipients',
  '% CCPs recipients',
  '% ACDs recipients',
  'No. of residents',
  'SNQ12',
  '% of residents who leave and go home',
  'Low mean household income',
  'Density of population',
  'Covid incidence at the municipality'
)

# heatmap:
pheatmap(aux,
         color = classic.colors(28),
         scale = "column",
         breaks = seq(-2.5,2.5,5/28),
         cellwidth =25,cellheight=20,
         kmeans_k = NA,
         show_rownames = T, show_colnames = T,
         main = "",
         clustering_method = "ward.D2",
         cluster_rows = F,
         cluster_cols = F,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         labels_row = labels_row,
         labels_col = labels_col,
         display_numbers = TRUE, number_color = "black",
         fontsize = 8,
         fontsize_number = 0,
         angle_col = 45,
         legend=T
)






# Variable importance at each cluster -----------------------
varImportance = lapply(1:n_distinct(resi.anl$louvain2),function(l){
  
  data =resi.anl %>% mutate(cluster = ifelse(louvain2==l,1,0) %>% to_dummy())
  rf = randomForest::randomForest(cluster~median_Age+Perc_Dones+median_comorbi+Residents_alfa+
                                    Perc_PCC+Perc_MACA+N_residents_actius+
                                    SNQ12+Perc_leave+renta_llar+N_pob+incidencia,
                                  data=data)
  labels_col = c(
                 'Elderliness',
                 '% of males',
                 'Number of comorbidities',
                 '% of Barthel score <50 recipients',
                 '% CCPs recipients',
                 '% ACDs recipients',
                 'No. of residents',
                 'SNQ12',
                 '% of residents who leave and go home',
                 'Low mean household income',
                 'Density of population',
                 'Covid incidence at the municipality'
  )
  aux = data.frame(vars=labels_col,importance=randomForest::importance(rf),
                   cluster=paste0('cluster ',l))
  row.names(aux)=1:nrow(aux)
  
  return(aux)
}) %>% bind_rows()

varImportance = varImportance %>% group_by(cluster) %>% 
  mutate_at('MeanDecreaseGini',fn) %>% data.frame
varImportance = varImportance %>% group_by(cluster) %>% 
  arrange(cluster,desc(MeanDecreaseGini)) %>%
  mutate(N=1:n())
varImportance$importance = ifelse(varImportance$N<=3,
                                  paste0(varImportance$MeanDecreaseGini,'%'),NA)
varImportance$col = ifelse(is.na(varImportance$importance),'black','black')
varImportance$size = ifelse(varImportance$N<=1,1,0.25)
varImportance$cluster = factor(varImportance$cluster,
                               levels = rev(levels(varImportance$cluster)))
varImportance$vars = factor(varImportance$vars,levels = rev(labels_col))

ggplot(varImportance,aes(x=cluster,y=MeanDecreaseGini,fill=vars))+
  geom_bar(stat='identity',size=varImportance$size,col=varImportance$col,alpha=.6)+
  scale_fill_brewer(palette='Set3')+
  geom_text(aes(label = importance),
            position=position_stack(0.5))+
  coord_flip()+
  blank_theme+ 
  guides(fill = guide_legend(reverse = TRUE))+
  theme(legend.title = element_blank())



# Barplot mortality at each cluster ---------------------------------------------
jet.colors = colorRampPalette(c(rgb(75,145,35,max=255),
                                    rgb(185,225,135,max=255),
                                    rgb(240,180,220,max=255),
                                    rgb(200,30,125,max=255)
))

resi.anl = pacients.anl %>% group_by(Codi) %>%
  summarise(N=n(),
            Perc_Exitus_covid_pre = perc(sum(Exitus_pre_post=='Pre' & type_exitus=='Covid'),N)
  ) %>% left_join(resi.anl)

data = resi.anl %>% 
  select(louvain=louvain2,
         Perc_Exitus,Perc_Exitus_covid) %>% na.omit %>%
  group_by(louvain) %>% reshape2::melt()

aux = data %>% group_by(louvain,variable) %>% 
  summarise(median = mean(value),sd=sd(value)/sqrt(n()))  %>%
  group_by(variable) %>% mutate(threshold = mean(median)) %>%
  ungroup() %>%
  mutate(intensity = ifelse((median-sd)>=threshold,'high',
                            ifelse(median>=threshold,'mid-high',
                                   ifelse(median+sd>=threshold,'mid-low','low'))),
         intensity = factor(intensity,levels=c('low','mid-low','mid-high','high')),
         colour = factor(intensity,levels=c('low','mid-low','mid-high','high'),
                         labels = jet.colors(4)),
         louvain = factor(louvain))
data = data %>% left_join(aux)
data$title = factor(data$variable,
                    levels = c('Perc_Exitus','Perc_Exitus_covid','Perc_Exitus_covid_pre'),
                    labels = c('% exitus','% exitus covid-19','% exitus covid-19 at baseline'))
aux$title = factor(aux$variable,
                    levels = c('Perc_Exitus','Perc_Exitus_covid','Perc_Exitus_covid_pre'),
                    labels = c('% exitus','% exitus covid-19','% exitus covid-19 at baseline'))


ggplot(data, aes(x=louvain, y=value,fill=intensity)) + 
  geom_boxplot(color="black",position=position_dodge()) +
  geom_hline(aes(yintercept = threshold),
             col='black',linetype='dotted',size=.5)+
  scale_color_manual(values=jet.colors(4))+
  scale_fill_manual(values=jet.colors(4))+
  facet_wrap(.~title,scales='free')+
  blank_theme+
  theme(legend.position = 'none')



# Barplot all variables included in clusterization -----------------------------------------------
data = resi.anl %>% 
  select(louvain=louvain2,median_Age,median_comorbi,
         Perc_Dones,Perc_PCC,Perc_MACA,
         Residents_alfa,incidencia,renta_llar,N_pob,
         N_residents_actius,Perc_leave,SNQ12) %>% na.omit %>%
  select(louvain,N_pob,everything()) %>%
  group_by(louvain) %>% 
  mutate(Residents_beta = 100-Residents_alfa,
         N_pob = log10(N_pob)) %>% 
  select(-Residents_alfa) %>%
  reshape2::melt()
data = data %>% mutate(variable=factor(variable,levels=c( 'median_Age',
                                                   'median_comorbi',
                                                   'Perc_Dones',
                                                   'Perc_PCC',
                                                   'Perc_MACA',
                                                   'Residents_beta',
                                                   'incidencia',
                                                   'renta_llar',
                                                   'N_pob',
                                                   'N_residents_actius',
                                                   'Perc_leave',
                                                   'SNQ12'),
                                       labels = c('Elderliness',
                                                  'Number of comorbidities',
                                                  '% of males',
                                                  '% CCPs recipients',
                                                  '% ACDs recipients',
                                                  '% of Barthel score <50 recipients',
                                                  'Covid incidence at the municipality',
                                                  'lack of household income at the municipality',
                                                  'Density of population at the municipality',
                                                  'No. of residents',
                                                  '% of residents who leave and go home',
                                                  'SNQ12'
                                                   )))

aux = data %>% group_by(louvain,variable) %>% 
  summarise(median = median(value),
            sd=sd(value)/sqrt(n()))  %>%
  group_by(variable) %>% mutate(threshold = mean(median)) %>%
  ungroup() %>%
  mutate(intensity = ifelse((median-sd)>=threshold,'high',
                            ifelse(median>=threshold,'mid-high',
                                   ifelse(median+sd>=threshold,'mid-low','low'))),
         intensity = factor(intensity,levels=c('low','mid-low','mid-high','high')),
         colour = factor(intensity,levels=c('low','mid-low','mid-high','high'),
                         labels = jet.colors(4)),
         louvain = factor(louvain))
data = data %>% left_join(aux)

p = ggplot(data, aes(x=louvain, y=value,fill=intensity)) + 
  geom_boxplot(color="black",position=position_dodge()) +
  geom_hline(aes(yintercept = threshold),
             col='black',linetype='dotted',size=.5)+
  scale_color_manual(values=jet.colors(4))+
  scale_fill_manual(values=jet.colors(4))+
  facet_wrap(~variable,ncol=3,scales='free')+
  blank_theme+
  theme(legend.position = 'none')
  





