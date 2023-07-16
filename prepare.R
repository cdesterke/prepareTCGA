### import patient data with removing header comment lines

patient<-read.table("data_clinical_patient.txt",h=T,sep="\t",comment.char="#",na.strings = c("NA" , "[Not Available]" , "[Not Applicable]"))

library(dplyr)

patient %>% mutate(OS_STATUS=recode(OS_STATUS, '1:DECEASED'='1', '0:LIVING'='0'))->patient

patient %>% mutate(DFS_STATUS=recode(DFS_STATUS, '1:Recurred/Progressed'='1', '0:DiseaseFree'='0'))->patient

patient$DFS_STATUS<-as.factor(patient$DFS_STATUS)

patient$OS_STATUS<-as.factor(patient$OS_STATUS)


### import sample data with removing header comment lines

sample<-read.table("data_clinical_sample.txt",h=T,sep="\t",comment.char="#",na.strings = c("NA" , "[Not Available]" , "[Not Applicable]"))


### joined sample and patient information

sample %>% inner_join(patient, by ="PATIENT_ID") -> annot

dim(annot)

### import data

library(data.table)

data<-fread("data.txt")

# remove rows with empty gene symbol
data<-data[!(is.na(data$Hugo_Symbol) | data$Hugo_Symbol==""), ]
data$Entrez_Gene_Id<-NULL

data %>% dplyr::rename(gene="Hugo_Symbol")->data

filtermatrix <- function(data)
{
  #install require R packages if necessary and load them
  if(!require(dplyr)){
    install.packages("dplyr")
    library(dplyr)
  }

  nbc=ncol(data)
  data2<-data[,2:nbc]
  m2 <- apply(data2,1, mean)
  combined2<-data.frame(m2,data$gene,data)
  ord2<-combined2[with(combined2,order(-m2)),]
  ok<-ord2 %>% distinct(gene, .keep_all = T)
  nbcok=ncol(ok)
  row.names(ok)<-ok$gene
  final<-ok[,4:nbcok]
  return(final)
}

ok<-filtermatrix(data)

trans<-as.data.frame(t(ok))

trans$SAMPLE_ID<-row.names(trans)
trans%>%relocate(SAMPLE_ID)->data

data$SAMPLE_ID<-gsub('\\.', '-',data$SAMPLE_ID)


### assemble all

annot %>% inner_join(data,by="SAMPLE_ID")->all

all[all == ''] <- NA

save(all,file="all.rda")

write.csv(all,file="all.csv",row.names=F)
