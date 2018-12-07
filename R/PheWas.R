#' Phewas
#' @param input1  The list of different clusters
#' @param input2 the data frame subset to define the prevalence of cluster subject
#' @param depth vector of disease symptoms that are used fro clustering
#' @param threshold the minimum prevalence threshold

#' @examples
#' example1 <- cluster_report(
#'                       clusterList     =clusterList,
#'                       input_subset    = GIAut,
#'                       clusterSubject  =clusterSubject
#'                       n               =5
#'                     )
#'
#'  @export DataInit,Phwas






DataInit=function(input1,input2,depth){
  DataMatrix=function(input)
  {
    D=unique(input[,c("PATIENT_NUM",depth,"START_DATE")])
    D=D[,c("PATIENT_NUM",depth)]
    D1=as.data.frame(table(D))
    D2=D1[D1[,depth]!="-",]
    D3=reshape2::acast(D2,D2$PATIENT_NUM~D2[,paste(depth)],new.row.names=FALSE)
    D4=as.data.frame(D3)
    return(D4)
  }

  DM1 = DataMatrix(input1)
  DM2=  DataMatrix(input2)
  Pheno=intersect(colnames(DM1),colnames(DM2))
  input1=DM1[,Pheno]
  input2=DM2[,Pheno]
  input1$label="A+"
  input2$label="A-"
  DT=rbind(input1,input2)
}













Phwas=function(DT,threshold){

  b=which(colnames(DT)=="label")
  input=data.frame(DT[,b],check.names = FALSE)
  output=data.frame(DT[,-b],check.names = FALSE)
  Pheno=colnames(output)
  V=parallel::mclapply(mc.cores=parallel::detectCores(),
             Pheno, function(x){
               row = Reg(x,input,output,threshold)
               return(row)
             })
  WW=data.table::rbindlist(lapply(V,function(x) data.frame(x)))
  # WW=WW[WW$A.OR.95..CI.!="(0-Inf)"]
}








# DT=DataInit(AD,AnD,"Phenotype")


Reg=function(x,input,output,threshold){
  newTable=data.frame()
  newTable=data.frame(cbind(output[,x],input),check.names = FALSE)
  # newTable=cbind(output[,x],input)
  names(newTable)[1]="x"
  names(newTable)[2]="condition"
  newTable=newTable[!(newTable[,1]>0 & newTable[,1]<threshold),]
  newTable[,1]=ifelse(newTable[,1] >=threshold,1,0)
  newTable=apply(newTable,2,function(x) as.factor(x))
  newTable=as.data.frame(newTable)
  fit <- glm(x~condition , data=newTable,family="binomial")
  pval=data.frame(summary(fit)$coefficients[,4])  #p value
  or=data.frame(round(exp(coef(fit)),2))          # odd ratio
  cI=round(exp(confint.default(fit, level=0.95)),2)
  a= nrow(newTable[newTable$condition=="A+" & newTable[,1]==1,])      # total ASD patients that have phenotype
  b= nrow(newTable[newTable$condition=="A-" & newTable[,1]==1,])      # total non-ASD patients that have phenotype
  c= nrow(newTable[newTable$condition=="A+" & newTable[,1]==0,])      # total ASD patients that doesnt have phenotype
  d= nrow(newTable[newTable$condition=="A-" & newTable[,1]==0,])
  row=list()
  row["Phenotype"]=paste(x)
  row["A+OR"]=or["conditionA+",]
  row["A+OR(95% CI)"]=paste0("(",cI["conditionA+",1],"-",cI["conditionA+",2],")")#OR for Condition
  row["Pval"]=pval["conditionA+",] #Pval for Condition
  row["tPheno+"]=paste0("(",a,"/",b,")")
  row["Pheno+"]=a+b
  row["tPheno-"]=paste0("(",c,"/",d,")")
  row["Pheno-"]=c+d
  return(row)}


Phecode=read.csv("./phecode_icd9_rolled.csv")
load("./pheinfo.rda")
# pheinfo[1815,]=c('563',"Constipation","10","digestive","chartreuse4")
# pheinfo[1816,]=c('561.1',"Diarrhea","10","digestive","chartreuse4")
# pheinfo[1817,]=c('279.2',"Immune disorders","3","endocrine/metabolic","brown")
# pheinfo[1818,]=c('300.3',"Obssesive compulsive disorder","5","mental disorders","magenta")
# pheinfo[1819,]=c('783',"symptoms concerning nutrition,mtblsm and Dvlpmnt","17","symptoms","purple")
# save(pheinfo,file="pheinfo (1).rda")

pheinfo$groupnum=as.numeric(pheinfo$groupnum)


library(devtools)
library(PheWAS)

PhePlot=function(input,refTable){
  names(input)[1]="description"
  regphe=merge(refTable,input,by="description")
  regphe1=regphe[,c(1:6,8)]
  names(regphe1)[1]="description"
  names(regphe1)[2]="phenotype"
  names(regphe1)[6]="OR"
  names(regphe1)[7]="p"
  regphe1$description=paste0(regphe1$description,"(",regphe1$OR,")")
  phenotypeManhattan(regphe1, suggestive.line = 1e-2, significant.line=1e-2/nrow(input),sort.by.category.value=FALSE,x.group.labels=TRUE,
                     OR.size = T, OR.direction = T, x.axis.label="Phenotype",annotate.angle=90,
                     y.axis.interval = 50, annotate.size=3,annotate.phenotype.description=TRUE,annotate.only.largest=TRUE,
                     switch.axis=FALSE) +theme(axis.text.x = element_text(size = 12, angle = 75, hjust =1.25, vjust = -0.5))


}






# SS=Phwas(DT=DT,threshold=3)
