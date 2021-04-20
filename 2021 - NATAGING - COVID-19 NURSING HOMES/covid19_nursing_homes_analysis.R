

# Nursing homes clusterization ----------------------------------------------------------------

data = dataset %>% 
  dplyr::select(Codi,Perc_Males=Perc_Dones,
                N_active_residents=N_residents_actius,
                median_Age,Perc_leave,n_comorb = median_comorbi,dependent_residents=Residents_alfa,
                CCP = Perc_PCC,ACP = Perc_MACA,Population_density=N_pob,household_income=renta_llar,
                snq12=SNQ12,incidence_area=incidencia) %>% na.omit 
data$Population_density = log10(data$Population_density)


# data pre-processing:

# 0,1 norm:
nor <-function(x){(x -min(x))/(max(x)-min(x))}
norm.data = data %>% dplyr::select(-(1:3)) %>% 
  mutate_if(is.Date,as.numeric) %>% 
  mutate_if(is.factor,as.numeric) %>% 
  mutate_all(nor) %>% data.frame %>%
  cbind.data.frame(data %>% dplyr::select(1:3),.)

# scaled:
data = data %>% dplyr::select(-(1:3)) %>% 
  mutate_if(is.Date,as.numeric) %>% 
  mutate_if(is.factor,as.numeric) %>% 
  scale %>% data.frame %>%
  cbind.data.frame(data %>% dplyr::select(1:3),.)
params = names(data)[-c(1:3)]


# correlation between variables:
corr = pairwise_simpleLM(data[-c(1:3)])[,c(1,2,9,11)]
uniexport(corr,'excel',path.output,'correlation_numeric',date=T)

pdf(paste0(path.output,'/Correlogram variables clustering.pdf'),height=15,width=15)
viRievac::correlogram(data[-c(1:3)],method = 'pearson',save = F)
dev.off()


# clustering:
Clusters = lapply(seq(4,30,1),function(i){
  pheno = Rphenograph::Rphenograph(data=data %>% dplyr::select(all_of(params)),k=i) # 8 clusters  
  
  pdf(paste0(path.output,'/network plot k_',i,'.pdf'),height=20,width=20)
  plot(pheno[[2]],pheno[[1]],vertex.label='')
  dev.off()
  
  data.frame(k=i,
             mean_distance=mean_distance(pheno[[1]]),
             modularity = modularity(pheno[[2]]),
             n_clusters=length(pheno[[2]]))
}) %>% bind_rows()

pdf(paste0(path.output,'/Number of clusters selection.pdf'),height=10,width=5)
par(mfrow=c(3,1))
plot(Clusters$k,Clusters$modularity,pch=20,ylab='Modularity',xlab='k value')
lines(Clusters$k,Clusters$modularity)
abline(v=12,lty=2,col='steelblue')

plot(Clusters$k,Clusters$mean_distance,pch=20,ylab='mean distance',xlab='k value')
lines(Clusters$k,Clusters$mean_distance)
abline(v=12,lty=2,col='steelblue')

plot(Clusters$k,Clusters$n_clusters+3,ylab='Number of clusters',xlab='k value',pch=19)
abline(v=12,lty=2,col='steelblue')
dev.off()







# Risk factors zero inflated poisson regression ---------------------------

## All-cause exitus: ----
# Univariate:
vars = dataset %>% 
  dplyr::select(median_Age,
                median_comorbi,
                Perc_Dones,
                Perc_PCC,
                Perc_MACA,
                Residents_alfa,
                incidencia,
                renta_llar_10,
                N_pob_log,
                N_residents_actius,
                Perc_leave,
                # Time_to_First_PCR,
                SNQ12
  ) %>% names 

out_zeroinflpoiss = lapply(1:length(vars),function(i){
  
  cat(vars[i],'\n')
  
  formula = as.formula(paste0('Exitus~',vars[i],'+offset(log(N_residents))|1'))
  model = zeroinfl(formula, data=dataset)
  cmodel = model %>% summary
  or = cbind.data.frame(exp(cmodel$coefficients$count[,1]),exp(confint(model))[1:2,],
                        cmodel$coefficients$count[,4])
  or = or %>% round(4)
  or = or[-1,]
  names(or) = c('OR','Lower','Upper','P-value')
  out = cbind.data.frame(Variable = vars[i],or)
  
  return(out)
}) %>% bind_rows()


# Multivariate:
m = as.formula(Exitus~Perc_Dones+Perc_PCC+incidencia+N_residents_actius+SNQ12+offset(log(N_residents))|1)
mod1 <- zeroinfl(m, data=dataset)
smod1 = summary(mod1)

or = cbind.data.frame(exp(smod1$coefficients$count[,1]),exp(confint(mod1))[1:6,],
                      smod1$coefficients$count[,4])
or = or %>% round(4)
or = or[-1,]
names(or) = c('adjOR','adjLower','adjUpper','adj P-value')
out = cbind.data.frame(Variable = rownames(or),or)


allcause_output = left_join(out_zeroinflpoiss,out)





## covid-19 related exitus: ----
# Univariate:
vars = dataset %>% 
  dplyr::select(median_Age,
                median_comorbi,
                Perc_Dones,
                Perc_PCC,
                Perc_MACA,
                Residents_alfa,
                incidencia,
                renta_llar_10,
                N_pob_log,
                N_residents_actius,
                Perc_leave,
                # Time_to_First_PCR,
                SNQ12
  ) %>% names 

out_zeroinflpoiss = lapply(1:length(vars),function(i){
  
  cat(vars[i],'\n')
  
  formula = as.formula(paste0('Exitus_covid~',vars[i],'+offset(log(N_residents))|1'))
  model = zeroinfl(formula, data=dataset)
  cmodel = model %>% summary
  or = cbind.data.frame(exp(cmodel$coefficients$count[,1]),exp(confint(model))[1:2,],
                        cmodel$coefficients$count[,4])
  or = or %>% round(4)
  or = or[-1,]
  names(or) = c('OR','Lower','Upper','P-value')
  out = cbind.data.frame(Variable = vars[i],or)
  
  return(out)
}) %>% bind_rows()

# Multivariate:
m = as.formula(Exitus_covid~Perc_PCC+incidencia+SNQ12+offset(log(N_residents))|1)
mod1 <- zeroinfl(m, data=dataset)
smod1 = summary(mod1)

or = cbind.data.frame(exp(smod1$coefficients$count[,1]),exp(confint(mod1))[1:4,],
                      smod1$coefficients$count[,4])
or = or %>% round(4)
or = or[-1,]
names(or) = c('adjOR','adjLower','adjUpper','adj P-value')
out = cbind.data.frame(Variable = rownames(or),or)


covid19_output = left_join(out_zeroinflpoiss,out)










