library(lattice)


#set this to change RRHO version
#packageToUse <- "RRHO2"
#packageToUse <- "Bioconductor"
packageToUse <- "LeonVersion"

if(packageToUse=="Bioconductor" ) {
  detach("package:RRHO", unload=TRUE)
  #Bio conductor RRHO:
  source("https://bioconductor.org/biocLite.R")
  biocLite("RRHO")
  library(RRHO)
  RRHO <- RRHO::RRHO
} else if(packageToUse=="LeonVersion" ) {
  detach("package:RRHO", unload=TRUE)
  #Corrected version from github
  library(devtools)
  install_github("leonfrench/RRHO-1")
  library(RRHO)
  RRHO <- RRHO::RRHO
}  else if(packageToUse=="RRHO2" ) {
  detach("package:RRHO", unload=TRUE)
  detach("package:RRHO2", unload=TRUE)
  library(devtools)
  install_github("Caleb-Huo/RRHO2")
  library(RRHO2)
  RRHO <- RRHO2
}



#two sorted lists of 4
list.length=4
list.names <- paste0("Gene", 1:list.length)
gene.list1<- data.frame(list.names, rank=1:list.length, stringsAsFactors = F)
gene.list2<- data.frame(list.names, rank=1:list.length, stringsAsFactors = F)

stepsize <- 1
alternative = "two.sided"
#alternative = "enrichment"

RRHOResultSameLists <- RRHO(gene.list1,gene.list2,alternative=alternative, stepsize = stepsize)
pValues <- exp(-RRHOResultSameLists$hypermat)
isSymmetric(pValues)
levelplot(pValues)

#one of the lists is reversed
list.names <- paste0("Gene", 1:list.length)
gene.list1<- data.frame(list.names, rank=1:list.length, stringsAsFactors = F)
gene.list2<- data.frame(list.names, rank=list.length:1, stringsAsFactors = F)

RRHOResultReverseList <- RRHO(gene.list1,gene.list2, alternative = alternative, stepsize = stepsize)
pValues <- exp(-RRHOResultReverseList$hypermat)
isSymmetric(pValues)
levelplot(pValues)


all.equal(RRHOResultReverseList$hypermat,RRHOResultSameLists$hypermat)
#below should be true (is not true if using the bioconductor RRHO)
all.equal(RRHOResultReverseList$hypermat[nrow(RRHOResultReverseList$hypermat):1,], RRHOResultSameLists$hypermat) #equal if is mirrored (ignores signs)


########################################
#create two random lists, compare outcomes
#similar to the original RRHO test case

list.length <- 100
list.names <- paste('Gene',1:list.length, sep='')
set.seed(1)
gene.list1<- data.frame(list.names, rank=sample(100),stringsAsFactors = F)
gene.list2<- data.frame(list.names, rank=sample(100),stringsAsFactors = F)
gene.list1 <- gene.list1[order(gene.list1$rank),]
gene.list2 <- gene.list2[order(gene.list2$rank),]

stepsize <- 20
alternative = "two.sided"

RRHOResult <- RRHO(gene.list1,gene.list2,alternative=alternative, stepsize = stepsize)
exp(-RRHOResult$hypermat)
levelplot(RRHOResult$hypermat)

gene.list2$rank <- rev(gene.list2$rank) #reverse the list
RRHOResultReverseList <- RRHO(gene.list1,gene.list2,alternative=alternative, stepsize = stepsize)
exp(-RRHOResultReverseList$hypermat)
levelplot(RRHOResultReverseList$hypermat)


all.equal(RRHOResultReverseList$hypermat[,ncol(RRHOResultReverseList$hypermat):1], RRHOResult$hypermat) #equal if is mirrored (ignores signs)
