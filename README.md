# A spatio-temporal analysis of NO2 concentrations during the Italian 2020 COVID-19 lockdown.

The following is the code for running INLA. The input dataset is available as a `.csv` file on this repository in the `data` folder.

```
rm(list=objects())
library("tidyverse")
library("INLA")
library("sf")
library("sp")
options(warn=0)


YEARS<-c(2019,2020)
MONTHS<-3:4

#read input data, this includes both March and April daily data
read_delim("data.csv",delim=";",col_names = TRUE)->mydata

#run INLA for each month
purrr::walk(MONTHS,.f=function(MM){

  #filter the input dataset in order to select only the month of interest
  mydata %>%
    filter(mm==MM)->dati
  
  #number of days
  nrow(dati[!duplicated(dati[,c("week","wday")]),])->n_giorni

  ######################################
  #Coordinates for each single station and daily observation: this is for the stack
  ######################################
  st_as_sf(dati,coords = c("coordx","coordy"),crs=32632)->sfDati
  
  #from meters to kilometers
  sf::st_transform(sfDati, crs = sp::CRS("+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs"))->sfDati
  as_tibble(st_coordinates(sfDati))->coordinateOsservazioni  #<- per lo stack 
  names(coordinateOsservazioni)<-c("coordx_km","coordy_km")
  bind_cols(sfDati[,c("station_eu_code","wday","week")],coordinateOsservazioni)->coordinateOsservazioni
  st_geometry(coordinateOsservazioni)<-NULL
  rm(sfDati)

  
  #Utilizzo coordinateOsservazioni per costruire la mesh: in questo caso ho bisogno di un punto per stazione (non per stazione E per osservazione)
  coordinateOsservazioni[!duplicated(coordinateOsservazioni$station_eu_code),c("coordx_km","coordy_km")]->stazioni
  as.matrix(stazioni)->stazioni
  
  ### Quale mesh?
  mesh<-inla.mesh.2d(loc =stazioni, max.edge = c(20,25),cutoff=10,min.angle = 30,offset=c(55,30)) 
  

  ######################## SPDE: Priors & more
  inla.spde2.pcmatern(mesh=mesh,alpha=2,constr=FALSE,prior.range = c(150,0.8),prior.sigma = c(0.5,0.2))->spde

  inla.spde.make.index(name="i",n.spde=spde$n.spde,n.group = n_giorni)->iset

  #training
  inla.spde.make.A(mesh=mesh,loc=as.matrix(coordinateOsservazioni[,c("coordx_km","coordy_km")]),group =dati$banda,n.spde=spde$n.spde,n.group =n_giorni )->A.training
  
  #Effects
  EFFETTI<-c("sp","t2m","tp","dtr","wspeed","pblmax","pblmin","pwspeed","rh","nirradiance","altitudedem","d_a2","clc_arable_agri","Intercept","station_eu_code","day","weekend")
  
  inla.stack(data=list(value=dati$value),A=list(A.training,1),effects=list(iset,dati[c(EFFETTI)]),tag="training")->stack.training

  stack.training->mystack

  ########################
  #PC Priors
  ########################
  list(theta = list(prior="pc.prec", param=c(1,0.1)))->prec_hyper  
  list(prior="pc.cor1",param=c(0.8,0.318))->theta_hyper #AR1

  formula(value ~ t2m + tp + dtr + wspeed + pblmax + pblmin + pwspeed + 
    rh + nirradiance + sp + altitudedem + d_a2 + clc_arable_agri + 
    Intercept + f(station_eu_code, model = "iid", hyper = prec_hyper) + 
    day + weekend +f(i, model = spde, group = i.group, control.group = list(model = "ar1",hyper = list(theta = theta_hyper))) - 1)->myformula

  ######################## INLA Vai!
  inla(myformula,
       data=inla.stack.data(mystack,spde=spde),
       family ="gaussian",
       verbose=TRUE,
       control.inla=list(int.strategy="eb"),
       control.compute = list(openmp.strategy="pardiso.parallel",cpo=F,waic=F,dic=F,config=TRUE),
       control.fixed = list(prec.intercept = 1, prec=0.01,mean.intercept=0),
       control.predictor =list(A=inla.stack.A(mystack),compute=TRUE) )->inla.out

  saveRDS(inla.out,glue::glue("result{MM}_{REGIONE}.RDS"))


}) #su MM
```
