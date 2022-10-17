###############################################################################
library(ape)
library(phangorn)
library(kmer)

#creates a panel of 1
layout(matrix(c(1))); layout.show(1);
db6 = read.FASTA("SeqDatabase_6.aln.fasta")
#creates the kdistance matrix for k=5
kd4.5 <- kdistance(db6, 7)
#plots the pcoa graph
plot(pcoa(kd4.5)$vectors[,1:2],main="PCOA of SeqDatabase_6.aln.fasta kdistance k=5");

###############################################################################

#imports necessary libraries
library(ape)
library(phangorn)
#sets working directory
setwd("~/Desktop/School/MCB432_folder/MA2Group")
#reads in db5
db6 = read.FASTA("SeqDatabase_6.aln.fasta")
db6 = as.phyDat(db6)
#tests the GTR model
a = modelTest(db6, model=c("GTR"))
#creates fitGTR
fitGTR <- as.pml(a, "GTR+G(4)+I")
#Optimizes using NNI rearrangement
fitGTR <- optim.pml(fitGTR,rearrangement="NNI")
#Creates bootstrapped tree
bsGTR10 <- bootstrap.pml(fitGTR,bs=10,optNni=TRUE,control = pml.control(trace=0))
layout(matrix(c(1),ncol=1));layout.show(1);par(mar = c(0,1,1,1))
#Graphs bootstrapped tree
plotBS(midpoint(fitGTR$tree),bsGTR10,p=50, type = "p", cex=0.55, main="ML:GTR + G(4) + I, bootstrapping = 10"); add.scale.bar(0,-2,font=2,col="red");







###############################################################################
#sets the working directory to the MA2 Group Folder
setwd("~/Desktop/School/MCB432_folder/MA2Group")
#Reads in the metadata table from producing SeqDatabase_6.fasta
meta = read.csv("SeqData6_meta.csv")
#reads in the clustering data from cd-hit, using line breaks a as a delimeter
clstr = read.table("tmp.99.clstr",sep='\n',header=F)
#converts the table file to a dataframe
clstr = data.frame(clstr)

#sets the column name to clustering center text
colnames(clstr) = "Clustering Center Text"

#initializes a counting variable and an empty list
c = 0;
listAccCent = (NULL);
#loops through all row numbers
for(k in 1:nrow(clstr)){
  #looks for rows that start with a >, this looks for the clustering center #
  if(substr(clstr[k,1],0,1) == ">") {
    #adds one to the counter variable
    c = c+1
  }
  #looks to see if this line conttains a substring ">NR"
  if(grepl(">NR",clstr[k,1],fixed = TRUE)){
    #initializes an empty integer variable
    ind = 0
    #searches through the line looking for substring of size 3 that is >NR
    for(s in 1:nchar(clstr[k,1])-3){
      #prints the substring
      print(substr(clstr[k,1],s,s+2))
      
      if( substr(clstr[k,1],s,s+2) == ">NR"){
        #sets the index variable equal to the index of wher the >NR starts
        ind = s
      }
    }
    #prints the index
    print(ind)
    #adds the substring of the line, the accession number to a list with the 
    #clustering center number
    a=list(c, substr(clstr[k,1],ind+1,ind+11))
    listAccCent = c(listAccCent,list(a))
  }
    
}
#initializes an empty column in the metadata table
meta[,8] = 0
#sets column name equal to "Clustering Numbers"
names(meta)[names(meta)=="V8"] <- "Clustering Numbers"


#loops through each object in the listAccCent
for(j in listAccCent){
  #loops through all rows in meta
  for(i in 1:nrow(meta)){
    print(substr(meta$ACC[[i]],1,11))
    print(j[[2]])
    #checks if substrings match
    if(substr(meta$ACC[[i]],1,11)==j[[2]]){
      #writes the clustering number to its appropriate row
      meta$`Clustering Numbers`[[i]] = j[[1]]
    }
  }

}

write.csv(meta,"metadata_Seq6.csv")

