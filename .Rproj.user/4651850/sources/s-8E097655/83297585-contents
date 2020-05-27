# GalgoR results functions

#' Title
#'
#' @param output
#' @param data
#' @param Surv
#' @param distancetype
#' @param usegpu
#'
#' @return
#' @export
#'
#' @examples
non.dominated.summary <- function(output,data, Surv, distancetype= "pearson",usegpu=FALSE){
  if (!methods::is(output, "galgo.Obj")) {
    stop("object must be of class 'galgo.Obj'")
  }

  require(survcomp)
  require(cluster)


  output_df<- toDataFrame(output)

  NonDom_solutions<- output_df[output_df$Rank==1,]

  select_distance(distancetype = distancetype,usegpu = usegpu)

  RESULT=data.frame(solution=as.character(),k=as.numeric(),ngenes=as.numeric(),mean.Silhouette=as.numeric(),C.Index=as.numeric(),stringsAsFactors = FALSE)

  for(i in 1:nrow(NonDom_solutions)){


    name<- rownames(NonDom_solutions)[i]
    genes<- NonDom_solutions[i,"Genes"][[1]]
    k<- NonDom_solutions[i,"k"]
    Sub_matrix<- data[genes,]

    D<- calculate_distance(Sub_matrix)

    true_class<- cluster_algorithm(D,k)

    Centroids <- kcentroid(Sub_matrix,true_class$cluster)

    predicted_class <- classify(Sub_matrix,Centroids,method= distancetype)

    predicted_class <- as.factor(predicted_class)
    predicted_classdf <- as.data.frame(predicted_class)


    surv_formula <- as.formula("Surv~ predicted_class")
    tumortotal <- survfit(surv_formula)
    totalsdf <- survdiff(surv_formula)
    tumortotalpval <- 1 - pchisq(totalsdf$chisq, length(totalsdf$n) - 1)
    tumortotalpval <- format(tumortotalpval, digits=4)

    coxsimple=coxph(surv_formula,data=predicted_classdf)

    CI=concordance.index(predict(coxsimple),surv.time=Surv[,1],surv.event=Surv[,2],outx=FALSE)$c.index
    mean_Sil =mean(silhouette(as.numeric(predicted_class),D)[,3])

    row= c(name,k,length(genes),mean_Sil,CI)
    RESULT[nrow(RESULT)+1,]=row
    #print(row)

  }
  return(RESULT)
}


create.centroids<- function(output,solution.names,train.set){
  CentroidsList<-list()
  output_df<- toDataFrame(output)
  for(j in solution.names){


    genes<- output_df[j,"Genes"][[1]]
    k<- output_df[j,"k"]
    name<- j
    Sub_matrix<- train.set[genes,]

    D<- calculate_distance(Sub_matrix)

    true_class<- cluster_algorithm(D,k)

    Centroids <- kcentroid(Sub_matrix,true_class$cluster)
    CentroidsList[[name]]=Centroids
  }
  return(CentroidsList)
}


classify.multiple<- function(data, centroid.list,distancetype="pearson"){

  classes<- matrix(rep(NA,ncol(data)*length(centroid.list)),ncol=length(centroid.list))
  as.data.frame<- classes
  colnames(classes)<- names(centroid.list)
  rownames(classes)<- colnames(data)

  for( i in 1:length(centroid.list)){

    centroids<- centroid.list[[i]]
    name<- names(centroid.list)[i]
    genes<- rownames(centroids)
    k<- ncol(centroids)
    Sub_matrix<- data[genes,]

    predicted_class <- classify(Sub_matrix,centroids,method= distancetype)

    predicted_class <- as.factor(predicted_class)
    classes[,name]=predicted_class

  }

  return(classes)
}


