## purl("GenExpr_Bioiberica_maig2019_aovambNA.Rnw",
## output = "GenExpr_Bioiberica_maig2019_aovambNA.R", documentation=1 )

  

## ----libs 
library(Biobase)
library(multtest)
library(genefilter)
library(limma)
library(DMwR)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(MASS)
library(gplots)
library(xtable)
require(FactoMineR)
require(missMDA)
require(readxl)
require(mice)
library(RColorBrewer)


## ----readdata 
dades<-read_excel("Exprgen_bioiberica_ili_jeju.xlsx",sheet = "Maindata")
dat<-dades

## Define colors associated to genetic functions, for heatmap, pca(),...
functions<-read_excel("Exprgen_bioiberica_ili_jeju.xlsx",sheet = "Funcions")
levels_func <- data.frame(geneid=functions$Gens,funcio=functions$Funcions, colors=0)
rownames(levels_func)<-functions$Gens

## assing colors to gen functions
for(i in 1:nrow(levels_func)){
   if(levels_func[i,2]=="BF"){
        levels_func[i,3]="dodgerblue"
   }else if(levels_func[i,2]=="OX"){
          levels_func[i,3]="firebrick1"
   }else if(levels_func[i,2]=="IR"){
          levels_func[i,3]="chartreuse1"
   }else if(levels_func[i,2]=="NT"){
          levels_func[i,3]="orange1"
   }else if(levels_func[i,2]=="EH"){
          levels_func[i,3]="gray48"
   }else{
          levels_func[i,3]="aquamarine1"   ## funció "S"
        }
}


## ----manipulateData 
 
## factors
dat$Mostra <- as.factor(dat$Mostra)
dat$Treat <- as.factor(dat$Treat)
dat$Bloc <- as.factor(dat$Bloc)
dat$Teixit <- as.factor(dat$Teixit)


# ORDER columns
## order gens by functions
id2 <- order(levels_func$funcio)
levels_func <- levels_func[id2,]
dat <- cbind(dat[,1:4],dat[,rownames(levels_func)])
## ordenem files de dat: teixit, tract, bloc  
id<-order(dat$Teixit,dat$Treat,dat$Bloc)   # index de files ordenat 
dat<-dat[id,]   

# subsets, by tissue
dat2 <- dat
datIli <- dat2[which(dat2$Teixit=="Ileum"),]      
datJeju <- dat2[which(dat2$Teixit=="Jejunum"),]   
rownames(datIli)<-paste0("Sample",datIli$Mostra)
rownames(datJeju)<-paste0("Sample",datJeju$Mostra)

## --------------------- function to manage NA's   # -----------------------#

gestioNA <- function( data, remove0=TRUE, del.badRows=TRUE, noNaMin=NULL)
{  
  ## columns containing factors
  classe<-character()
  factor <- numeric()
  for (i in 1:ncol(data)){ classe[i]<-class(data[,i])
  if (classe[i]=="factor") factor<-c(factor,i)
  } 
  ## change 0 to NA, if exist
  if (remove0==TRUE) data[!is.na(data)&data==0] <- NA  ## .............!!!
  
  ## delete columns all NA's!
  nacol <- apply(data,2,function(x) sum(is.na(x)))
  if (sum(nacol==nrow(data))>0)
  { 
    elim<-which(nacol==nrow(data),useNames=TRUE) 
    data<-data[-elim]
  }  
  ## delete rows all NA's!
  narow <- apply(data[-factor],1,function(x) sum(is.na(x)))  # suma de NA per files
  if (sum(narow==ncol(data[-factor]))>0)
  { 
    delrow<-which(narow==ncol(data[-factor]),useNames=TRUE) 
    data<-data[-delrow,]
  }  
  ## delete rows with 50% or more NA's!  (opcional: del.badRows==TRUE)
  if (del.badRows==TRUE)
  { 
    condition  <- narow >= floor(0.50*ncol(data[-factor]))   # narrow >= numgens/2
    whichcond  <- which(condition==TRUE,useNames=TRUE)  
    if(length(whichcond)>0) data<-subset(data,condition==FALSE)   
  }
  ## delete columns (gens) with a too small number of valid replicates
  ### in one or more treatments
  ## default threshold: integer part (number of replicates /2)
  if (is.null(noNaMin)) noNaMin<- floor(max(table(data$Treat))/2)        
  
  ## delete rows s.t. number of valid replicates in some treatment is < noNaMin
  llista<-split(data,data$Treat)
  noNASxTrac<-t(sapply(llista,function(df)apply(df,2,function(v)sum(!is.na(v)))))
  eliminar<-  which( apply(noNASxTrac< noNaMin,2,sum) > 0, useNames=TRUE)
  data2<-data[-eliminar]
  
  return(data2)
}

# --------------------------------- end function --------------------------#


datIlina <- gestioNA(datIli)     
datJejuna <- gestioNA(datJeju)  


## ----log, 

classe<-character()
factors <- numeric()
for (i in 1:ncol(dat)){ classe[i]<-class(dat[,i])
                           if (classe[i]=="factor") factors<-c(factors,i)
} 
iden<-which(colnames(dat)=="Mostra")   ## preserving the original samples numbers  
 
datIlinal <- cbind(datIlina[,factors],log(datIlina[,-factors],10))    
rownames(datIlinal) <- paste0("Sample", datIlinal[,iden]) ## keeping original samples' names

datJejunal <- cbind(datJejuna[,factors],log(datJejuna[,-factors],10))   
rownames(datJejunal) <- paste0("Sample", datJejunal[,iden])


## ----ExpressionSets  

## using BioBase library to create an object of class ExpressionSet, allowing to implement
### the package functions

## database=datIlinal
## sepVar="Treat"
## returning object of class expressionset  
 
GetDataExpression <- function(database, sepVar){
  df <- data.frame(Treat=database[,sepVar],
                 row.names=rownames(database))
  dataExpression <- ExpressionSet(as.matrix(t(database[,sapply(database,function(x) is.numeric(x))])),
                    phenoData = AnnotatedDataFrame(data=df))
  return(dataExpression)
}

# apply GetDataExpression to the files
datIliExpression  <- GetDataExpression(datIlinal,"Treat")
datJejuExpression <- GetDataExpression(datJejunal,"Treat")

## Functions ro compute ANOVA statistic and p-value and FDR corrected p-value
## return a 3-colum dataframe

# -------------------  start functions -----------------------·

## gives up output with NAs !

getF<-function(x,trac,dades,alfa)
{ summary(aov(x~trac, data=dades))[[1]]['F value'][1,1]}

getFpval<-function(x,trac,dades,alfa) 
  {summary(aov(x~trac, data=dades))[[1]]['Pr(>F)'][1,1]}

## gives up output with NAs
GetPvalProva <- function(dades=datJejunal,sepvar='Treat',fac=factors,alpha=0.1)
 {
  dades0<-dades[-factors]
  dades1<- data.frame(dades[sepvar],dades0)
  statistic<-apply(dades1[-1],2, function(x){ 
                       trac <- dades1[,1]
                       summary(aov(x~trac))[[1]]['F value'][1,1]
                       })
  p.value<-apply(dades1[-1],2, function(x){
                       trac <- dades1[,1]
                       summary(aov(x~trac))[[1]]['Pr(>F)'][1,1]
                       })
  tt<-data.frame(statistic,p.value)          
  p.BH = p.adjust(tt[,"p.value"], "BH" )
  tt <- cbind(tt,p.BH)
  g <- which( tt[,2] <= alpha)
  if( length(g) > 0){
  tt[which( tt[,2] <= alpha),] 
  }else{
    print("No hi ha cap significatiu")
  }
}

#--------------------------------- end functions --------------#

require(xtable)
tab1<-GetPvalProva(dades=datJejunal,sepvar='Treat',fac=factors,alpha=0.1)
tab3<-GetPvalProva(dades=datIlinal,sepvar='Treat',fac=factors,alpha=0.1)

tab1.tot <- GetPvalProva(dades=datJejunal,sepvar='Treat',fac=factors,alpha=1) 
tab3.tot <- GetPvalProva(dades=datIlinal,sepvar='Treat',fac=factors,alpha=1)

rnames1<-rownames(tab1)    ## significant genes, jejunum
rnames3<-rownames(tab3)    ## significant genes, ilium 

asig.funcio<-function(nomgens)
{## nomgens<-rnames1 
n<-length(nomgens)
nomfuncio<-character(n)  # i=1
for (i in 1:n) {
   num<-which(levels_func$geneid==nomgens[i])
   nom<-as.character(levels_func$funcio[num])
   nomfuncio[i]<-nom}
nomfuncio
}

fnames1<-asig.funcio(rnames1)
fnames3<-asig.funcio(rnames3)

tab1   <-data.frame(fnames1,tab1)
tab3   <-data.frame(fnames3,tab3)

colnames(tab1)   <- c("gen-func","statistic", "p-value", "p-value FDR") 
colnames(tab3)   <- c("gen-func","statistic", "p-value", "p-value FDR")



## ----xtables 

bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}
gray <- function(x) {paste('{\\textcolor{gray}{',x,'}}', sep ='')}
print(xtable(tab1, digits= 6, caption ="P-valors de l'ANOVA entre tractaments pel teixit \\textbf{Jejunum}, es mostren només els gens amb diferències significatives. Les dades faltants (NA) només afecten al gen on hi ha el valor perdut."),
      sanitize.rownames.function=gray, 
      sanitize.colnames.function=bold, 
      booktabs=T)
print(xtable(tab3, digits= 6, caption ="P-valors de l'ANOVA entre tractaments pel teixit \\textbf{Ileum}, es mostren només els gens amb diferències significatives.  Les dades faltants (NA) només afecten al gen on hi ha el valor perdut."), 
      sanitize.rownames.function=gray, 
      sanitize.colnames.function=bold, 
      booktabs=T)


## ---- adding the genes names to the tables and export to excel

## 
## tab1.ambnoms<-data.frame(gen=rownames(tab1),tab1)
## tab1.totambnoms<-data.frame(gen=rownames(tab1.tot),tab1.tot)
## tab3.ambnoms<-data.frame(gen=rownames(tab3),tab3)
## tab3.totambnoms<-data.frame(gen=rownames(tab3.tot),tab3.tot)
## 
## library(writexl)
## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019/aovNA")
## write_xlsx(x = tab1.ambnoms, path = "tabAnovaJeju_aovNA.xlsx", col_names = TRUE)
## write_xlsx(x = tab1.totambnoms, path = "tabAnovaJeju_tot_aovNA.xlsx", col_names = TRUE)
## write_xlsx(x = tab3.ambnoms, path = "tabAnovaIli_aovNA.xlsx", col_names = TRUE)
## write_xlsx(x = tab3.totambnoms, path = "tabAnovaIli_tot_aovNA.xlsx", col_names = TRUE)
## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019")


## ----heatmaps.jeju 

rwb      <- colorRampPalette(colors = c("darkred", "white", "darkgreen"))
BrBG     <- colorRampPalette(brewer.pal(11, "BrBG"))
Spectral <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

da0<-datJejunal
da.clus<-t(da0[-factors])

color.map <- function(x){ if(x=="1") "black" else if(x=="2") "orange" else "red" }  
patientcolors <- unlist(lapply(da0$Treat, color.map))
funcColors    <- levels_func[rownames(da.clus),3]
uniquecolors  <- unique(funcColors)
uniquefunc    <- unique(as.character(levels_func[rownames(da.clus),2]))

## select distances and clustering method
## for rows
r.row<-cor(t(da.clus),use="pairwise.complete.obs")
d.row <- as.dist(0.5*(1-r.row))
distance.row = d.row
cluster.row = hclust(distance.row, method = "complete")
## for columns
distance.col = dist(t(da.clus), method = "euclidean")
cluster.col = hclust(distance.col, method = "ward.D2")

## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019/aovNA")
## tiff("heatmapJeju_aovNA.tiff",height =7.5, width = 7.5, units = 'in', type="cairo", res=300)
  #width = 480, height = 480, units = "px", pointsize = 12,type="tifflzw", res=300)

heatmap.2(da.clus, scale="none", trace="none", dendrogram="both", 
           Rowv=as.dendrogram(cluster.row), Colv=as.dendrogram(cluster.col),
          col=Spectral,
          na.color="dimgray",na.rm=FALSE,cexRow=.75,cexCol=.75,
          ColSideColors=patientcolors, RowSideColors = funcColors, 
          #hclustfun=function(x) hclust(x, method="ward.D2"),
          symkey=FALSE, density.info="none")
coords <- list(x=-0.13,y=-0.08)  
legend(coords,  horiz = T,xpd=T,
       legend = uniquefunc,
       col = uniquecolors, 
       lty= 1,             
       lwd = 4,           
       cex= 0.41
)
coords2 <- list(x=-0.09,y=0.80) ## COMPILAR rnw
legend(coords2,  horiz = T,xpd=T,
       #legend = c("T1","T2","T3"),
       legend = c("CON2","SDP2","SDP-PDP"),
       col = c("black","orange","red"),  
       lty= 1,             
       lwd = 4,           
       cex= 0.41
)

## dev.off()
## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019")


## ----heatmaps.ili

rwb      <- colorRampPalette(colors = c("darkred", "white", "darkgreen"))
BrBG     <- colorRampPalette(brewer.pal(11, "BrBG"))
Spectral <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

da0<-datIlinal
da.clus<-t(da0[-factors])

color.map <- function(x){ if(x=="1") "black" else if(x=="2") "orange" else "red" }  
patientcolors <- unlist(lapply(da0$Treat, color.map))
funcColors    <- levels_func[rownames(da.clus),3]
uniquecolors  <- unique(funcColors)
uniquefunc    <- unique(as.character(levels_func[rownames(da.clus),2]))

## select distances and clustering method
## for rows
r.row<-cor(t(da.clus),use="pairwise.complete.obs")
d.row <- as.dist(0.5*(1-r.row))
distance.row = d.row
cluster.row = hclust(distance.row, method = "complete")
## for columns
distance.col = dist(t(da.clus), method = "euclidean")
cluster.col = hclust(distance.col, method = "ward.D2")

## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019/aovNA")
## tiff("heatmapIli_aovNA.tiff",height =7.5, width = 7.5, units = 'in', type="cairo", res=300)
  #width = 480, height = 480, units = "px", pointsize = 12,type="tifflzw", res=300)

heatmap.2(da.clus, scale="none", trace="none", dendrogram="both", 
           Rowv=as.dendrogram(cluster.row), Colv=as.dendrogram(cluster.col),
          col=Spectral,
          na.color="dimgray",na.rm=FALSE,cexRow=.75,cexCol=.75,
          ColSideColors=patientcolors, RowSideColors = funcColors, 
          #hclustfun=function(x) hclust(x, method="ward.D2"),
          symkey=FALSE, density.info="none")
coords <- list(x=-0.13,y=-0.08)  
legend(coords,  horiz = T,xpd=T,
       legend = uniquefunc,
       col = uniquecolors, 
       lty= 1,             
       lwd = 4,           
       cex= 0.41
)
coords2 <- list(x=-0.13,y=0.80)  
legend(coords2,  horiz = T,xpd=T,
       #legend = c("T1","T2","T3"),
       legend = c("CON2","SDP2","SDP-PDP"),
       col = c("black","orange","red"),     ## canviats colors
       lty= 1,             
       lwd = 4,           
       cex= 0.41
)
## dev.off()
## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019")


## ----tukey functions
 
Tukey_test<- function(dataExpression){
  Tuk <- list()
  for(i in 1:nrow(exprs(dataExpression))){
    if( sum(is.na(exprs(dataExpression)[i,])) < length(exprs(dataExpression)[i,])-2){
      Tuk[[i]] <- TukeyHSD(aov(exprs(dataExpression)[i,]~ pData(dataExpression)[,1]))
               rownames(exprs(dataExpression))[i]
      names(Tuk)[[i]] <- rownames(exprs(dataExpression))[i]
    }
  }
  return(na.omit(Tuk))
}

GetTukeyPlot <- function(DataTuk, noms){
  ###
  pd <- position_dodge(0.78)
  Tractaments <- factor(rownames(DataTuk$`pData(dataExpression)[, 1]`),labels=c(rownames(DataTuk$`pData(dataExpression)[, 1]`)))
  ggplot(as.data.frame(JejuTuk[[1]]$`pData(dataExpression)[, 1]`), aes(x=Tractaments, y = diff, group = Tractaments)) +
   #draws the means
      geom_point(position=pd) +
   #draws the CI error bars
      geom_errorbar(data=as.data.frame(DataTuk$`pData(dataExpression)[, 1]`), aes(ymin=lwr, ymax=upr, 
      color=1), width=.1, position=pd)+ggtitle(DataTuk) + theme(plot.title = element_text(size = 0.5, face = "bold"))+
  xlab("Tractaments") + ylab("Diff means")
  
}



## ----tukeyJeju 
 
## datJejuExpression<- GetDataExpression(datJejunal,"Treat")  ## fet abans
JejuTuk    <- Tukey_test(datJejuExpression)
JejuNames  <- rownames(tab1)
JejuTukGood  <- list()
for(i in 1:length(JejuNames)){ JejuTukGood[[i]] <- JejuTuk[[JejuNames[i]]]}
dataTukJeju <- as.data.frame(t(sapply(JejuTukGood,function(x) x$`pData(dataExpression)[, 1]`[,4], USE.NAMES = F)))

nomfunc<-asig.funcio(JejuNames)
Tukjeju_nom<-data.frame(gen=JejuNames,func=nomfunc,dataTukJeju)  # nom gen sense nom funció
jejunomscomplets<-paste(nomfunc,'.',JejuNames,sep="")
rownames(dataTukJeju) <- jejunomscomplets                       

## mitjanes i sd de cada gen, entre els significatius

taulaTukey<-function(dades=datJejunal, taula=tab1,  # "taula" conté els significatius
                     dataTuk=dataTukJeju,nomscomplets=jejunomscomplets)
{
 aux<- dades
 taulanames<-rownames(taula)
 aux_sig<- aux[,colnames(aux) %in% c("Treat",taulanames)]
  
 aux_llista<- split(aux_sig[-1],aux_sig$Treat)  
 aux_means <- as.data.frame(sapply(aux_llista, function(df)colMeans(df,na.rm=TRUE)))
 aux_sd    <- as.data.frame(sapply(aux_llista, function(df) apply(df,2,sd,na.rm=TRUE))) 
 colnames(aux_means)<-paste0("T",1:ncol(aux_means),".mean")
 colnames(aux_sd)<-paste0("T",1:ncol(aux_sd),".sd")
 
 aux_descr<-data.frame(aux_means,aux_sd)
 nc<-ncol(aux_descr)
 reord<-c()
 for (i in 1:(nc/2)) reord<-c(reord,i,i+3)
 aux_descr<-aux_descr[reord]   ## reordenat
 datatable<-cbind(aux_descr,dataTuk)
 rownames(datatable) <- nomscomplets   
 datatable
 }

datatableJ<-taulaTukey() 

print(xtable(datatableJ,digits=3,caption="Tukey:comparacions múltiples Gene-Tractament"),
      sanitize.rownames.function=gray, 
      sanitize.colnames.function=bold, 
      booktabs=T,
      scalebox='0.5')

## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019/aovNA")
## png('TukeyJeju_aovNA.## png')
plotrs <- list()
for(i in 1:length(JejuTukGood)){
  DataTuk <- JejuTukGood[[i]]
  pd <- position_dodge()
  Tractaments <- factor(rownames(DataTuk$`pData(dataExpression)[, 1]`),labels=c(rownames(DataTuk$`pData(dataExpression)[, 1]`)))
  plotrs[[i]]<- ggplot(as.data.frame(JejuTuk[[1]]$`pData(dataExpression)[, 1]`), aes(x=Tractaments, y = diff, group = Tractaments)) +
   # horizontal line at 0
      geom_hline(yintercept=0, linetype="dashed", color = "red") +
   #draws the CI error bars
      geom_errorbar(data=as.data.frame(DataTuk$`pData(dataExpression)[, 1]`), aes(ymin=lwr, ymax=upr, 
      color=1), width=.1, position=pd)+ggtitle(JejuNames[i]) + theme(axis.line=element_blank(),,
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none")
}

do.call(grid.arrange,plotrs)
## dev.off()
## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019")


## ----tukeyIli 

IliTuk <- Tukey_test(datIliExpression)
IliNames <- rownames(tab3)
IliTukGood <- list()
for(i in 1:length(IliNames)){IliTukGood[[i]] <- IliTuk[[IliNames[i]]]}
dataTukIli <- as.data.frame(t(sapply(IliTukGood,function(x) x$`pData(dataExpression)[, 1]`[,4], USE.NAMES = F)))
## afegit de nou  -------------
nomfunc<-asig.funcio(IliNames)
ilinomscomplets<-paste(nomfunc,'.',IliNames,sep="")
rownames(dataTukIli) <- ilinomscomplets
 

datatableI<-taulaTukey(dades=datIlinal, taula=tab3,  # "taula" conté els significatius
                     dataTuk=dataTukIli,nomscomplets=ilinomscomplets)

print(xtable(datatableI,digits=3,caption="Tukey: comparacions múltiples Gene-Tractament"),
      sanitize.rownames.function=gray, 
      sanitize.colnames.function=bold, 
      booktabs=T,
      scalebox='0.5')

## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019/aovNA")
## png('TukeyIli_aovNA.## png')
plotrs <- list()
for(i in 1:length(IliTukGood)){ # i=1
  DataTuk <- IliTukGood[[i]]  
  pd <- position_dodge()
  Tractaments <- factor(rownames(DataTuk$`pData(dataExpression)[, 1]`),labels=c(rownames(DataTuk$`pData(dataExpression)[, 1]`)))
  plotrs[[i]]<- ggplot(as.data.frame(IliTuk[[1]]$`pData(dataExpression)[, 1]`), aes(x=Tractaments, y = diff, group = Tractaments)) +
   # horizontal line at 0
      geom_hline(yintercept=0, linetype="dashed", color = "red") +
   #draws the CI error bars
      geom_errorbar(data=as.data.frame(DataTuk$`pData(dataExpression)[, 1]`), aes(ymin=lwr, ymax=upr, 
      color=1), width=.1, position=pd)+ggtitle(IliNames[i]) + theme(axis.line=element_blank(),,
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none")
}

do.call(grid.arrange,plotrs)
## dev.off()
## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019")



## ---- export tukey tables to excel
## ##install.packages("writexl")
## 
## datatableJ.ambnoms<-data.frame(gen=rownames(datatableJ),datatableJ)
## datatableI.ambnoms<-data.frame(gen=rownames(datatableI),datatableI)
## 
## library(writexl)
## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019/aovNA")
## write_xlsx(x = datatableJ.ambnoms, path = "tabTukeyJeju_aovNA.xlsx", col_names = TRUE)
## write_xlsx(x = datatableI.ambnoms, path = "tabTukeyIli_aovNA.xlsx", col_names = TRUE)
## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019")
## 


## ----  starting Lines plots -------------------------------------------------# 

## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019/aovNA")
## tiff("Linies_ordreFamTracJeju_aovNA.tiff",height =7.5, width = 7.5, units = 'in', type="cairo", res=300)


### ----jejunum

Tr <- datJejunal$Treat
ntr<- length(unique(Tr))

bytratJeju<-split(datJejunal, Tr)
names(bytratJeju)<-paste0("T", 1:ntr)

de<- data.frame(sapply(bytratJeju,function(df){colMeans(df[-factors],na.rm=T)}))
deJeju<-de

nomfunc<-asig.funcio(rownames(de))
denomcompl<-paste(nomfunc,rownames(de),sep=".")
rownames(de)<-denomcompl
de$funcio<-nomfunc
de<-de[order(de$funcio,-de$T1),]
 
rang.y <- range(de[1:ntr]) 
plot(de$T1,type="o",pch=19,col="black",xaxt='n',xlab=NA, ylab="Treatments' mean-expression",
     ylim=rang.y+c(-0.2,+0.2))
axis(side=1, at=1:nrow(de), labels=rownames(de),las=2, cex.axis=0.5)
lines(de$T2,type="o",pch=19, col="orange")
lines(de$T3, type="o",pch=19,col="red")
 
id<-1:nrow(de)
tab1nomcompl<-paste(tab1[,1],'.',rownames(tab1),sep="")
gens_sig<-id[rownames(de)%in% tab1nomcompl]  ## significant genes
for (i in 1:length(gens_sig)) abline(v=gens_sig[i],col="gray40",lty="dotted")
legend("topright",c("CON2","SDP2","SDP-PDP"), cex=0.8, 
       col=c("black","orange","red"),lty=1, title="Treatment",bg="white")
## dev.off()
## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019/")


## ---- 
## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019/aovNA")
## tiff("Linies_ordreTracJeju_aovNA.tiff",height =7.5, width = 7.5, units = 'in', type="cairo", res=300)

de <- de[order(de$T1,decreasing = T),]    
rang.y <- range(de[1:ntr]) 
plot(de$T1,type="o",pch=19,col="black",xaxt='n',xlab=NA, ylab="Treatments' mean-expression",
     ylim=rang.y+c(-0.2,+0.2))
axis(side=1, at=1:nrow(de), labels=rownames(de),las=2, cex.axis=0.5)
lines(de$T2,type="o",pch=19, col="orange")
lines(de$T3, type="o",pch=19,col="red")
 # legend
id<-1:nrow(de)
tab1nomcompl<-paste(tab1[,1],'.',rownames(tab1),sep="")
gens_sig<-id[rownames(de)%in% tab1nomcompl]  ## significant genes 
for (i in 1:length(gens_sig)) abline(v=gens_sig[i],col="gray40",lty="dotted")
legend("topright",c("CON2","SDP2","SDP-PDP"), cex=0.8, 
       col=c("black","orange","red"),lty=1, title="Treatment",bg="white")
## dev.off()
## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019")


## ---- ileum

## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019/aovNA")
## tiff("Linies_ordreFamTracIli_aovNA.tiff",height =7.5, width = 7.5, units = 'in', type="cairo", res=300)

Tr <- datIlinal$Treat
ntr<- length(unique(Tr))

bytratIli<-split(datIlinal, Tr)
names(bytratIli)<-paste0("T", 1:ntr)

de<- data.frame(sapply(bytratIli,function(df){colMeans(df[-factors],na.rm=T)}))
deIli<-de

nomfunc<-asig.funcio(rownames(de))
denomcompl<-paste(nomfunc,rownames(de),sep=".")
rownames(de)<-denomcompl
de$funcio<-nomfunc
de<-de[order(de$funcio,-de$T1),]
 

rang.y <- range(de[1:ntr]) 
plot(de$T1,type="o",pch=19,col="black",xaxt='n',xlab=NA, ylab="Treatments' mean-expression",
     ylim=rang.y+c(-0.2,+0.2))
axis(side=1, at=1:nrow(de), labels=rownames(de),las=2, cex.axis=0.5)
lines(de$T2,type="o",pch=19, col="orange")
lines(de$T3, type="o",pch=19,col="red")
 
id<-1:nrow(de)
tab3nomcompl<-paste(tab3[,1],'.',rownames(tab3),sep="")
gens_sig<-id[rownames(de)%in% tab3nomcompl]  ## significant genes
for (i in 1:length(gens_sig)) abline(v=gens_sig[i],col="gray40",lty="dotted")
legend("topright",c("CON2","SDP2","SDP-PDP"), cex=0.8, 
       col=c("black","orange","red"),lty=1, title="Treatment",bg="white")
## dev.off()
## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019")
 

## ---- 
## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019/aovNA")
## tiff("Linies_ordreTracIli_aovNA.tiff",height =7.5, width = 7.5, units = 'in', type="cairo", res=300)

de<-de[order(-de$T1),]

rang.y <- range(de[1:ntr]) 
plot(de$T1,type="o",pch=19,col="black",xaxt='n',xlab=NA, ylab="Treatments' mean-expression",
     ylim=rang.y+c(-0.2,+0.2))
axis(side=1, at=1:nrow(de), labels=rownames(de),las=2, cex.axis=0.5)
lines(de$T2,type="o",pch=19, col="orange")
lines(de$T3, type="o",pch=19,col="red")
 
id<-1:nrow(de)
tab3nomcompl<-paste(tab3[,1],'.',rownames(tab3),sep="")
gens_sig<-id[rownames(de)%in% tab3nomcompl]  ## significant genes
for (i in 1:length(gens_sig)) abline(v=gens_sig[i],col="gray40",lty="dotted")
legend("topright",c("CON2","SDP2","SDP-PDP"), cex=0.8, 
       col=c("black","orange","red"),lty=1, title="Treatment",bg="white")
## dev.off()
## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019")


## ----  ending Lines plots   -------------------------------------------------# 


## ---- functions to apply PCA() procedures

## assigning functions names to each gene
## a 2-parameter version of asig.funcio
asig.func<-function(genes,funs)
{ 
 n<-length(genes)
 nomfuncio<-character(n)  # i=1
 for (i in 1:n) 
     {
      num<-which(funs$geneid==genes[i])
      nom<-as.character(funs$funcio[num])
      nomfuncio[i]<-nom
     }
nomfuncio
}


## assigning colors to gens
asig.col<-function(nomfunc,lfuns) 
{
 n<-length(nomfunc)
 colorf<-character(n)  # i=1
 for (i in 1:n) 
   {
    num<-which(lfuns$funcio==nomfunc[i])[1]
    color<-as.character(lfuns$colors[num])
    colorf[i]<-color
   }
 colorf
}


## ---- PCA jejunum variables=treatmeans

dades<-deJeju   
nomfunc<-asig.funcio(rownames(dades))
dades$func <- as.factor(nomfunc)
colnames(dades)<-c("CON2","SDP2","SDP-PDP")
require(FactoMineR)
pcajetr<-PCA(dades,quali.sup=4,graph=F)

## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019/aovNA")
##  tiff("PCAJeju_ambNA.tiff",height =3.8, width = 7.5, units = 'in', type="cairo", res=300)

par(mfrow = c(1,2),
          oma = c(0,0,0,0) + 0.7,
          mar = c(4,4,3,2) + 0.2,
    lwd=1.9, cex.main=.8, cex=.7,  cex.lab=.8,cex.axis=.8)  

## variables plot (treatments means plot)
plot(pcajetr,choix="var",col.var="blue")

nomgenaux<-rownames(pcajetr$ind$coord)
funcColInd <- as.character(levels_func[nomgenaux,3])
nomsquali<-rownames(pcajetr$quali.sup$coord)
cols<-asig.col(nomsquali,levels_func) 

## individuals plot
plot.PCA(pcajetr,choix="ind",col.ind=funcColInd, invisible="quali",label="none")
legend("topright",nomsquali, cex=0.6, col=cols,lty=1,bg="white")

## dev.off()
## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019")


## ---- PCA ileum variables=treatmeans

dades<-deIli  ## té valors imputats
nomfunc<-asig.funcio(rownames(dades))  # one parameter version
dades$func <- as.factor(nomfunc)
colnames(dades)<-c("CON2","SDP2","SDP-PDP")

require(FactoMineR)
pcailtr<-PCA(dades,quali.sup=4,graph=F)

## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019/aovNA")
## tiff("PCAIli_ambNA.tiff",height =3.8, width = 7.5, units = 'in', type="cairo", res=300)

par(mfrow = c(1,2),
         oma = c(0,0,0,0) + 0.7,
         mar = c(4,4,3,2) + 0.2,
    lwd=1.9, cex.main=.9, cex=.7,  cex.lab=.8,cex.axis=.8) 

## variables plot (treatments means plot)
plot(pcailtr,choix="var",col.var="blue") 

nomgenaux<-rownames(pcailtr$ind$coord)
funcColInd <- as.character(levels_func[nomgenaux,3])
nomsquali<-rownames(pcailtr$quali.sup$coord)
cols<-asig.col(nomsquali,levels_func) 

## individuals plot
plot.PCA(pcailtr,choix="ind",col.ind=funcColInd,invisible="quali",label="none") 
legend("topright",nomsquali, cex=0.6, col=cols,lty=1,bg="white")
 
## dev.off()
## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019")


####################   standard PCA -----------------------------------#

## ---- PCA variables=gens, individuals=samples -----------------------
 

###  define function

funcpca<-function(data,eixos=c(1,2),nomstr=nomstracs,gsignif=tab1, nivellsfunc=levels_func,
                  coltract=c("black","orange","red"), gruix=1.7, gris=4, limcos2=0.5,
                  cextit=.7,cexleg=.6,cexlab=.6,cexax=.6,cexlletra=.7)
{ # data=datJejunal[,-c(1:3)]
  # gsignif=tab1  
  # eixos=c(1,3)
 datas<-data
 colnames(datas)[1] <- "Treatment"                   
 
 require(FactoMineR)
 require(missMDA)
 ## impute NAs with EM-algorithm to improve default imputation by mean value
 
 if (sum(is.na(datas))!=0){
   datas.i_1<-imputePCA(datas[-1])$completeObs
   # sum(is.na(datas.i_1)) ## must be 0
   datas.i<-data.frame(datas[1],datas.i_1)
   }
 else datas.i<-datas
   
 pcaout<-PCA(datas.i,quali.sup=1,graph=F)
 
 nomgen<-colnames(datas.i[-1])  ## only genes
 nomfun<-asig.func(genes=nomgen,funs=nivellsfunc)
 
 par(oma = c(0,0,0,0) + 0.5,
         mar = c(4,4,1,2) + 0.2)
 
 ## individuals graph
 nlev<-length(nomstr)
 plot(pcaout,choix="ind",axes=eixos, invisible="quali", habillage=1, col.hab=coltract,
     lwd=gruix, cex.main=cextit, cex=cexlletra, cex.lab=cexlab, cex.axis=cexax)  
 legend("topleft", nomstr, cex=cexleg, col=coltract,lty=1, title="Treatment",bg="white")
 
 ## variables graph
 ## only genes with a quality over threshold in limcos2
 ###  and being significant, that is, belong to the table gsignif
 
 pcaux<-pcaout
 pcaqual2<-apply(pcaux$var$cos2[,eixos],1,sum) 
 ## subset of pcaout
 ### significant genes and quality over threshold
 pcaux$var$coord<-subset(pcaux$var$coord,rownames(pcaux$var$coord)%in% rownames(gsignif)& pcaqual2 > limcos2) 
 pcaux$var$cos2 <-subset(pcaux$var$cos2,rownames(pcaux$var$cos2)%in% rownames(gsignif)& pcaqual2 > limcos2)
 
# gray color is the baseline color for variables, overprinted with functions colors afterwards
 colorbase<-gray.colors(1,0.3+gris*0.05,0.3+gris*0.1)
 plot(pcaux, axes=eixos, choix="var", col.var=colorbase, 
      lwd=gruix, cex.main=cextit, cex=cexlletra*.8,  cex.lab=cexlab,cex.axis=cexax)  

# same colors as in heatmap
 ng<-nrow(pcaux$var$coord)
 nomgenaux<-rownames(pcaux$var$coord)
 nomfuncaux<-asig.func(nomgenaux,nivellsfunc)
 colaux<-character()
 for (i in 1:ng){## i<-1
  fila<-which(rownames(levels_func)==nomgenaux[i])
  colaux[i]<-levels_func[fila,3]           
  }
 for (i in 1:ng)
  arrows(x0=0,y0=0,x1=pcaux$var$coord[i,eixos[1]],
      y1=pcaux$var$coord[i,eixos[2]],col=colaux[i], angle = 14,length=.1)
# legend
 colors<-unique(colaux)
 funcaux<-unique(nomfuncaux)
 legend("topleft",legend=funcaux, col=colors,cex=cexleg,lty=1,title="Functions")
} 
 

## ---- PCA jejunum variables=gens ...

## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019/aovNA")

dades<-datJejunal[,-c(1:3)]  ## delete factors columns
factortrac<-dades[,1]
## class(factortrac); levels(factortrac)
nlev<-length(levels(factortrac))
nomstracs<-c("CON2","SDP2","SDP-PDP")

## tab1: contains the significant genes in Jenunum

## tiff("nouPCAJeju_aovNA_12.tiff",height =4.4, width = 7.3, units = 'in', type="cairo", res=300)

par(mfrow=c(1,2))
funcpca(data=dades,gsignif=tab1,gruix=1.9, gris=2,limcos2=0.5,cextit=.6,cexax=.6,cexleg=.65,cexlletra=.4)

## dev.off()
## tiff("nouPCAJeju_aovNA_13.tiff",height =4.4, width = 7.3, units = 'in', type="cairo", res=300)

par(mfrow=c(1,2))
funcpca(data=dades,gsignif=tab1,gruix=1.9, gris=2,limcos2=0.5,eixos=c(1,3),cextit=.6,cexax=.6,cexleg=.65,cexlletra=.4)

## dev.off()
## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019")


## ---- PCA ileum variables=gens ...

## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019/aovNA")

dades<-datIlinal[,-c(1:3)]  ## trec columnes dels factors, només queda tract i gen
factortrac<-dades[,1]
nlev<-length(levels(factortrac))
nomstracs<-c("CON2","SDP2","SDP-PDP")
## tab3: taula anova+pvalor+FDR de gens significatius imputats

##tiff("nouPCAiLI_aovNA_12.tiff",height =4.5, width = 7.3, units = 'in', type="cairo", res=300)
par(mfrow=c(1,2))
funcpca(data=dades,gsignif=tab3,gruix=1.9, gris=2,limcos2=0.5,cextit=.6,cexax=.6,cexleg=.65,cexlletra=.4)
## dev.off()
## tiff("nouPCAiLI_aovNA_13.tiff",height =4.5, width = 7.3, units = 'in', type="cairo", res=300)
par(mfrow=c(1,2))
funcpca(data=dades,gsignif=tab3,gruix=1.9, gris=2,limcos2=0.5,eixos=c(1,3),cextit=.6,cexax=.6,cexleg=.64,cexlletra=.4)
## dev.off()

## setwd("C:/Users/farre/Desktop/PINSO/ExprGENETICA/DOCS.Finals/bioiberica/fitxers14maig2019")

