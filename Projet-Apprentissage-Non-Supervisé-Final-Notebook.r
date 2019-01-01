# -*- coding: utf-8 -*-
#' ---
#' jupyter:
#'   jupytext:
#'     text_representation:
#'       extension: .R
#'       format_name: spin
#'       format_version: '1.0'
#'       jupytext_version: 0.8.5
#'   kernelspec:
#'     display_name: R
#'     language: R
#'     name: ir
#'   language_info:
#'     codemirror_mode: r
#'     file_extension: .r
#'     mimetype: text/x-r-source
#'     name: R
#'     pygments_lexer: r
#'     version: 3.5.1
#' ---

#' # Apprentisage Non-Supervisé (Notebook Version)

#' ## Bibliothèques R

#install.packages(c("bigstatsr",
#                   "doParallel",
#                    "parallel",
#                    "factoextra",
#                    "FactoMineR",
#                    "NbClust",
#                    "ggplot2",
#                    "dplyr",
#                    "ggrepel",
#                    "foreach",
#                    "MASS",
#                    "matrixStats",
#                    "devtools",
#                    "Rcpp"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if(!require("RDRToolbox", quietly = TRUE)){
    BiocManager::install("RDRToolbox", version = "3.8")
}else 
    library("RDRToolbox")


library("bigstatsr")
library("doParallel")
library("parallel")
library("factoextra")
library("FactoMineR")
library("NbClust")
library('ggplot2')
library('dplyr')
library("ggrepel")
library("foreach")
library("MASS")
library("matrixStats")
library(caret)
library(e1071)
library(scales)

#' ## Chargement des jeux de données


sequencesimu = "sequencesimu.txt"
aiguillage = "aiguillage.txt"


data_sequencesimu <- read.table(sequencesimu,header = FALSE)
data_sequencesimu_mat <- as.matrix(data_sequencesimu)
nbSeqSimu = length(data_sequencesimu_mat)
meanNbSeqSimu = mean(1:nbSeqSimu)
ecNbSeqSimu = sqrt(var(1:nbSeqSimu))


data_aiguillage <- read.table(aiguillage,header = FALSE,sep=",")
data_aiguillage_mat <- as.matrix(data_aiguillage)
data_aiguillage_variable = as.matrix(data_aiguillage[,-553]) 
data_aiguillage_classe = as.matrix(data_aiguillage[,553])
data_aiguillage_classeDF = data_aiguillage[,-553]
NbAiguillage=nrow(data_aiguillage_variable)

#' # Implémentation de l'Algorithme de Fisher

#' ### La fonction Diamètre

#' Version 1D


diam <- function(donnees){ 
  n = nrow(donnees)
  D <- matrix(data = 0, nrow = n, ncol = n)    
  for (a in 1:(n-1)){ 
    for (b in (a+1):n){
      D[a,b] <- var(donnees[a:b,])*(b-a+1) 
    }
  }
  D
}

#' Version nD

diamND_slow <- function(donnees){ 
    n = nrow(donnees)
    D <- matrix(0,n,n)
    foreach(a=1:(n-1)) %do% {
        foreach(b=(a+1):n) %do% {
          D[a,b]= sum(var(donnees[a:b,])*(b-a+1))
        }
    }
    D
}

#' version optimisé avec matrixStats:colVars


diamND <- function(donnees){ 
    n = nrow(donnees)
    D <- matrix(0,n,n)
    foreach(a=1:(n-1)) %do% {
        foreach(b=(a+1):n) %do% {
          D[a,b]= sum(matrixStats::colVars(donnees[a:b,])*(b-a+1))
        }
    }
    D
}

#' ### La fonction clustFisher


clustFisher = function(data,kmax=2,diamF=diamND){
    n=nrow(data)

    D=diamF(data)
    M1=matrix(data=0,nrow=n,ncol=kmax)
    M2=matrix(data=0,nrow=n,ncol=kmax)
    M3=matrix(data=0,nrow=n,ncol=kmax)
    M1[,1] = D[1,]

    for(k in 2:kmax){
        for(i in k:n){
            # skip la recherche des sous groupes kmax + 1
            if(k == kmax && i != n) next
            # skip la derniere boucle des sous groupes  quand k<kmax 
            #if(k<kmax && i == n)next
            val=c()

            for( t in k:i){
                val=c(val,M1[t-1, k-1] + D[t, i])
            }
            M1[i,k] = min(val)
            M2[i,k] = which.min(val) + k - 1
            M3[i,k] = M1[i,k]
        }

    }
    
    op=getVarCluster(list(M2=M2,M3=M3,data=data),kmax)
    t=op$t
    vari = op$var

    cluster = getCluster(list(data=data,t=t))
    
    res = list()
    res$t=t
    res$M1=M1
    res$M2=M2
    res$M3=M3
    res$cluster=cluster
    res$data=data
    res$diam=D
    res$var = vari
    res
}


#' get cluster from data and t

getCluster = function(cf){
    o=0
    curr=0
    n=nrow(cf$data)
    cluster=c()
    for (i in cf$t){
        cluster=c(cluster,rep(o,i-curr))
        curr=i
        o=o+1
    }
    cluster=c(cluster,rep(o,n-curr))
    cluster
}

#' get wss and t for k (var,t) <br> need only the result of clustFisher with the param k <br>
#' with param withClusters = T  -> Call the function `getCluster`

getVarCluster = function(cf,k,withClusters=FALSE){
    n=nrow(cf$data)
    t=rep(0,time=(k-1))
        k=k-1
        m=n
        vari = NULL
        while(k >=1) {
           t[k] =  cf$M2[m, k +1 ] - 1
           if(is.null(vari)) vari = cf$M3[m, k  +1]
           m = t[k] 
           k=k-1
        }
    res=list()
    res$var=vari
    res$t=t
    if(withClusters) res$clusters= getCluster(list(t=t,data=cf$data))
    res
}

#' ### Fonction pour trouver les variances pour toutes les classes jusqu'a n_clust


getVarAndFisher=function(data,n_clust=10){
    clust10 = clustFisher(data,n_clust)
    frz= sapply(2:n_clust,function(i)getVarCluster(clust10,i)$var)
    df=data.frame(unlist(matrix(frz)))
    df=cbind(2:n_clust,df)
    colnames(df) = c("index","V1")
    rownames(df)=df$index
    list(clust10=clust10,df=df)
}

#' # SequenceSimu

#' ## Comparaison des 3 algorithmes CAH Ward, Kmeans, et l'algorithme de Fisher

#' Il y a 210 lignes et 1 colonnes


dim(data_sequencesimu_mat)

#' ## Graphiques Du jeux de données  

#' On arrive deja à avoir une certaine idée des clusters (5 clusters) avec le graphique


data_sequencesimu %>% ggplot(aes(x=1:nbSeqSimu,y=V1)) +
                                                    geom_point() +
                                                    xlab("t") + 
                                                    ylab("value") + 
                                                    ggtitle("SequenceSimu") + 
                                                    theme(plot.title = element_text(hjust = 0.5))

#' ## Application de l'algo de Fisher

#' ### Données Brutes


n_clust=10


getBru = function(force=F,autoSave=T){
    if(file.exists("bru.Rdata")){
        if(!force){
            load("bru.Rdata")
            return(bru)
        }
    }
    bru=getVarAndFisher(data_sequencesimu_mat,n_clust)
    if(autoSave)save(bru,file="bru.Rdata")
    return(bru)
}
bru = getBru()


bru


df=bru$df
clustBru=bru$clust10

#' La variance de chaque cluster allant de 2 à 10 


df[1:10,]

#' ### Méthode du coude

#' ### Données Brutes


df %>% ggplot(aes(x=2:n_clust,y=V1)) + 
                                    ggtitle("Méthode du Coude avec Fisher") + 
                                    geom_line(col="blue")  + 
                                    geom_point() + 
                                    xlab("Nombre de Clusters") +
                                    ylab("wss") +
                                    scale_x_continuous(breaks = 2:n_clust) +
                                    geom_vline(xintercept = 5,col="red") + 
                                    geom_label_repel(aes(label=round(V1,2))) + 
                                    theme(plot.title = element_text(hjust = 0.5))

#' On voit qu'il faut choisir 5 pour le choix du nombre de classes

#' ### Representation des clusters

#' Fonction pour l'affichage


addLine = function(ggplotFig,meanV=meanNbSeqSimu,ecV=ecNbSeqSimu,title=""){
    plotFig=ggplotFig + 
    geom_point() + 
    xlab("row") + 
    ylab("value")
    for(i in cl$t){
        plotFig=plotFig+geom_vline(xintercept = (i-meanV)/ecV,col="red")
    }
     
    plotFig+
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5)) 
                                                        
}

#' Recuperation des clusters pour 5 


cl=getVarCluster(clustBru,5,withClusters = T)


addLine(fviz_cluster(show.clust.cent = T,labelsize = 0,list(data = as.data.frame(
                                        cbind(1:nrow(data_sequencesimu_mat),data_sequencesimu_mat)
                                      ), 
                  cluster = cl$cluster
                ))
        ,title="Nuage de point avec les différentes séquences")

#' ### Kmeans sur le jeux de données de Simulation 


nbc=fviz_nbclust(data_sequencesimu,kmeans,k.max = n_clust,method = "wss")


nbc$data


nbc + 
    geom_vline(xintercept = 4,col="red") + 
    ggtitle("Méthode du Coude avec Kmeans sur le jeux de données simulation") + 
    geom_label_repel(aes(label=round(nbc$data$y,2))) +
    theme(plot.title = element_text(hjust = 0.5)) 

#' Nous utiliserons le nombre de cluster fixé à 4

cluster_kmeans=kmeans(data_sequencesimu,4)

fviz_cluster(list(data = as.data.frame(cbind(1:nbSeqSimu,data_sequencesimu_mat)), 
                  cluster = cluster_kmeans$cluster),
             main="Nuage de points et cluster associés (Kmeans et jeux de données simulation)")

#' ### La classification ascendante hiérarchique (Ward)


nbc=fviz_nbclust(data_sequencesimu,hcut,k.max = n_clust,method = "wss")


nbc + 
    geom_vline(xintercept = 4,col="red") + 
    ggtitle("Méthode du Coude avec CAH Ward sur les données de simulation") + 
    geom_label_repel(aes(label=round(nbc$data$y,2))) +
    theme(plot.title = element_text(hjust = 0.5)) 

#' Nous utiliserons le nombre de cluster fixé à 4


hv=hcut(data_sequencesimu,4)


fviz_dend(hv,4,main="Dendogramme avec une CAH Ward sur le jeux de données simulation")


fviz_cluster(list(data = as.data.frame(cbind(1:length(data_sequencesimu_mat),data_sequencesimu_mat)), 
                  cluster = hv$cluster),
             main="Nuage de points et cluster associés (CAH Ward et jeux de données simulation)")

#' # Aiguillage

#' ## Comparaison des 3 algorithmes CAH Ward, Kmeans, et l'algorithme de Fisher

#' Il y a 140 lignes et 553 colonnes


dim(data_aiguillage_mat)

#' Il y a 4 classes


as.data.frame(table(data_aiguillage_classe))

qplot(data_aiguillage_classe, geom="histogram") +
ggtitle("Histogramme des classes d'aiguillage") + 
theme(plot.title = element_text(hjust = 0.5))

#' ## Graphiques

#' V1 au cours du temps


as.data.frame(data_aiguillage_variable) %>% ggplot(aes(x=1:NbAiguillage,y=V1)) + 
                                                                geom_point() + 
                                                                xlab("row") + 
                                                                ylab("V1") + 
                                                                ggtitle("Aiguillage") + 
                                                                theme(plot.title = element_text(hjust = 0.5))

#' Visuellement nous pouvons constater qu'il y a 4 séquences.

#' ## Algorithme de Fisher


n_clust = 10

test = clustFisher(as.matrix(data_aiguillage_variable),4)

test$cluster

opF = getVarAndFisher(as.matrix(data_aiguillage_variable),n_clust)

df = opF$df
clustFish = opF$clust10

cl=getVarCluster(clustFish,4,withClusters = T)

# Table de confusion 
tab <- table(test$cluster+1,data_aiguillage_classe)
print(tab)
prop.table(tab)

d=data.frame(x=1:nrow(data_aiguillage_classe),y=cl$clusters+1)
d %>% ggplot() + geom_point(aes(x=x,y=y),color=data_aiguillage_classe) + ggtitle("Séquences d'Aiguillage (Vrais partitions)") + 
                                                                theme(plot.title = element_text(hjust = 0.5))  

d=data.frame(x=1:nrow(data_aiguillage_classe),y=data_aiguillage_classe)
d %>% ggplot() + geom_point(aes(x=x,y=y),color=data_aiguillage_classe)

#' ### Kmeans avec 4 clusters sur les données Aiguillage

#' D'après la méthode du coude, nous aurions tendance à dire qu'il y 2 groupement.

cluster_kmeans=kmeans(data_aiguillage_variable,4)

cluster_kmeans$cluster

clustFish

d = table(cluster_kmeans$cluster,data_aiguillage_classe)
d

# Table de confusion 
tab <- table(cluster_kmeans$cluster,data_aiguillage_classe)
print(tab)
prop.table(tab)

d=data.frame(x=1:nrow(data_aiguillage_classe),y=cluster_kmeans$cluster)
d %>% ggplot() + geom_point(aes(x=x,y=y),color=data_aiguillage_classe)  + ggtitle("Séquences d'Aiguillage (K-means)") + 
                                                                theme(plot.title = element_text(hjust = 0.5)) 

hv=hcut(data_aiguillage_variable,4)

fviz_dend(hv,main="Dendogramme découpage (Ward) sur le jeux de données Aiguillage")

table(hv$cluster)

tab <- table(hv$cluster,data_aiguillage_classe)
print(tab)
prop.table(tab)

d=data.frame(x=1:nrow(data_aiguillage_classe),y=hv$cluster)
d %>% ggplot() + geom_point(aes(x=x,y=y),color=data_aiguillage_classe) + ggtitle("Séquences d'Aiguillage (CAH Ward)") + 
                                                                theme(plot.title = element_text(hjust = 0.5)) 

valid_actual <- as.factor(data_aiguillage_classe)

valid_pred   <- as.factor(cluster_kmeans$cluster)

cfm <- confusionMatrix(valid_actual, valid_pred)

ggplotConfusionMatrix <- function(m){
  mytitle <- paste("Confusion Matrix for Kmeans:","Accuracy", percent_format()(m$overall[1]),
                   "Kappa", percent_format()(m$overall[2]))
  data_c <-  mutate(group_by(as.data.frame(m$table), Reference ), percentage = 
percent(Freq/sum(Freq)))
  p <-
    ggplot(data = data_c,
           aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = log(Freq)), colour = "white") +
    scale_fill_gradient(low = "white", high = "green") +
    geom_text(aes(x = Reference, y = Prediction, label = percentage)) +
    theme(legend.position = "none") +
    ggtitle(mytitle)
  return(p)
}

ggplotConfusionMatrix(cfm)

valid_actual <- as.factor(data_aiguillage_classe)

valid_pred   <- as.factor(hv$cluster)

cfm <- confusionMatrix(valid_actual, valid_pred)

ggplotConfusionMatrix <- function(m){
  mytitle <- paste("Confusion Matrix for CAH Ward:","Accuracy", percent_format()(m$overall[1]),
                   "Kappa", percent_format()(m$overall[2]))
  data_c <-  mutate(group_by(as.data.frame(m$table), Reference ), percentage = 
percent(Freq/sum(Freq)))
  p <-
    ggplot(data = data_c,
           aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = log(Freq)), colour = "white") +
    scale_fill_gradient(low = "white", high = "green") +
    geom_text(aes(x = Reference, y = Prediction, label = percentage)) +
    theme(legend.position = "none") +
    ggtitle(mytitle)
  return(p)
}

ggplotConfusionMatrix(cfm)

valid_actual <- as.factor(data_aiguillage_classe)

valid_pred   <- as.factor(cl$clusters+1)

cfm <- confusionMatrix(valid_actual, valid_pred)

ggplotConfusionMatrix <- function(m){
  mytitle <- paste("Confusion Matrix for Fisher:","Accuracy", percent_format()(m$overall[1]),
                   "Kappa", percent_format()(m$overall[2]))
  data_c <-  mutate(group_by(as.data.frame(m$table), Reference ), percentage = 
percent(Freq/sum(Freq)))
  p <-
    ggplot(data = data_c,
           aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = log(Freq)), colour = "white") +
    scale_fill_gradient(low = "white", high = "green") +
    geom_text(aes(x = Reference, y = Prediction, label = percentage)) +
    theme(legend.position = "none") +
    ggtitle(mytitle)
  return(p)
}

ggplotConfusionMatrix(cfm)

c_fast = vector(length=5)
for (a in 2:6){ 
      c_fast[a-1] <- system.time(clustFisher(as.matrix(data_aiguillage_variable),a,diamF=diamND))
  }

c_fast

c_slow = vector(length=5)
for (a in 2:6){ 
      c_slow[a-1] <- system.time(clustFisher(as.matrix(data_aiguillage_variable),a,diamF=diamND_slow))
  }

c_slow

mean(c_slow/c_fast)

row_names = c(2,3,4,5,6)

row_names

r <- data.frame(matrixStats=c_fast,var=c_slow,number_clusters=c(2:6))

r

r %>% ggplot() + geom_line(aes(x=r$number_clusters,y=r$var,colour = "var"))  + 
                 geom_line(aes(x=r$number_clusters,y=r$matrixStats,colour = "MatrixStats")) +
                 ggtitle("Benchmarks") + 
                 xlab("Number of clusters") + 
                 ylab("Seconds") + 
                 theme(plot.title = element_text(hjust = 0.5)) 

ggplot(r, aes(number_clusters)) + 
  geom_line(aes(y = matrixStats, colour = "MatrixStats")) + 
  geom_line(aes(y = var, colour = "var"))
