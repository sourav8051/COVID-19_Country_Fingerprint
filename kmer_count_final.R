## kmer count using kmer and ape package 06.03.20
kmer_count=function(new_dir,metadata)
  #inpud: put all fasta files in a folder named 'dataset' and km=no of k
{

#load libraries
library(ape)
library(readr)
library(kmer)
#library(reshape2)

#Step1: download FASATA files from NCBI https://www.ncbi.nlm.nih.gov/labs/virus/vssi/

#Step2: Go to that folder containing those FASTA files
#system('cd dataset')
setwd(new_dir)
#km=3  #no of k in k-mer 
km <- readline(prompt="Enter value of k ")
all_km=0
for (i in 1:km)
{
  all_km= all_km+ (4^i)
}


files=metadata$Accession
no_of_files=NROW(files)

kmer_result=NULL #matrix(data = NA, nrow = no_of_files, ncol = all_km)

for (i in 1:no_of_files)
{
  
  
  
  x=suppressMessages(read_csv(paste0(files[i],".fasta"), col_names = FALSE))
  val=gsub(".fasta","",files[i])
  
  y=NULL
  for (j in 1: length(x[[1]]))
  {
    y=paste0(y,x[[1]][j])  #merging all rows into one row
  }
  
  y <- t(sapply(strsplit(y,""), tolower)) #split the row into n no of columns
  rownames(y) <- val #save the seq name as rowname
  z=as.DNAbin(y)  #convert to DNAbin using ape 
  
  tmp=NULL
  for (j in 1:km)
  {
    tmp=cbind(tmp,kcount(z,k=j)) #kmer package
  }
tmp=tmp/length(y)

kmer_result=rbind(kmer_result,tmp) 
 #Sys.sleep(0.05)
 #setTxtProgressBar(pb, (i/28)*100)   # update progress bar

}

#close(pb)
#browser()

##check the ordering of kmer_result and metadata
View(rownames(kmer_result)==metadata$Accession)


View(table(metadata$Country)) #use table or prop.table function to get the %
View(prop.table(table(metadata$Country)))

#merge classes with k_mer result
kmer_result=cbind(kmer_result,as.matrix(metadata$Country))
kmer_result=cbind(kmer_result,as.matrix(metadata$Collection_Date))
kmer_result=cbind(kmer_result,as.matrix(metadata$Months))
kmer_result=cbind(kmer_result,as.matrix(metadata$Weeks))
kmer_result=cbind(kmer_result,as.matrix(metadata$Continent))
kmer_result=cbind(kmer_result,as.matrix(metadata$Col_Rel_week_diff))
 


#write.csv(as.character(country_class),"country_class.csv")
write.csv(kmer_result,paste0("kmer_result_",km,".csv"))
#write.csv(as.character(release_class),"release_class.csv")
View(kmer_result)
return(kmer_result)

}

read_multiple_fasta=function()

{
#open multiple seq file from one FASTA 08.03.20
#read the file
#sequences <- read_csv(<"filename.fasta">, col_names = FALSE)
  
fname <- file.choose()
library(readr)
sequences <- suppressMessages(read_csv(fname, col_names = FALSE))
#read the no of seq files and locations from this multiple seq fasta file using "<... pattern

fasta_seq=which(substr(sequences$X1,1,1)=='>')  #locations from <....
n=length(fasta_seq)    #the number of sequences in one fasta file

#km <- readline(prompt="Enter value of k ")

#split the file into n number of small FASTA files in a folder
#create a new directory in local machine and put all those files there
cur_dir=dirname(fname)
new_dir=paste0(cur_dir,"/Result_",Sys.Date())
dir.create(new_dir)
setwd(new_dir)

for (i in 1:n) {
  seq_name=stringr::str_extract(sequences$X1[fasta_seq[i]],">........")   #Extrace the first 9 characters: ">LR757995"
  seq_name=sub(">","",seq_name)  #remove the first >  : "LR757995"
  filename=paste0(seq_name,".fasta")
  if (i==n) {
    write(sequences$X1[fasta_seq[i]:length(sequences$X1)],filename)
    } else{ 
    write(sequences$X1[fasta_seq[i]:(fasta_seq[i+1]-1)],filename) 
    }
}



#kmer_result=kmer_count(km,new_dir)
#setwd(cur_dir)
#write.csv(kmer_result,paste0(fname,"_",km,".csv"))
#return(kmer_result)

return(new_dir)

}



data_preprocessing2=function()
{
  library(readr)
  library(tidyr)
  
  print('Enter the metadata table from NCBI in xls format')
  fname <- file.choose()
  setwd(dirname(fname))  #set the working dir as the input file dir
  
  metadata<- suppressMessages(read_csv(fname))  
  
  ##separate Geo_location into country and states
  metadata_copy=metadata
  metadata=separate(metadata,Geo_Location,c("Country","States"),sep=":")
  
  
  #change the date column
  metadata$Collection_Date=as.Date.character(metadata$Collection_Date,format = "%d-%m-%y")
  metadata$Collection_Date=format(as.Date(metadata$Collection_Date), "%d %B %Y")
  
  metadata$Release_Date=format(as.Date(metadata$Release_Date), "%d %B %Y")
 
  #add the months class
  Months=as.Date.character(metadata_copy$Collection_Date,format = "%d-%m-%y")
  Months=format(as.Date(Months), "%B %Y")
  metadata=cbind(metadata,Months)
  
  #add the weeks class 
  Weeks=as.Date.character(metadata_copy$Collection_Date,format = "%d-%m-%y")
  Weeks=format(as.Date(Weeks),"%W")
  metadata=cbind(metadata,Weeks)
  
  ###get continent class
  library(countrycode)
  Continent <- data.frame(country = metadata$Country)
  Continent <- countrycode(sourcevar = Continent[, "country"],
                              origin = "country.name",
                              destination = "continent")
  metadata=cbind(metadata,Continent)
  
  #delete rows with NA in dates
  t=which(is.na(metadata$Collection_Date))
  metadata=metadata[-t,]
  
  #difference between collection and release dates 
  t1=strptime(metadata$Release_Date, format = "%d %B %Y")
  t2=strptime(metadata$Collection_Date, format = "%d %B %Y")
  Col_Rel_week_diff=difftime(t1,t2,units = "weeks")
  metadata=cbind(metadata,Col_Rel_week_diff)
  
  
  fn=format(Sys.Date(),"%d_%b")
  write.csv(metadata,paste0("processed_metadata_",fn,".csv"),row.names = FALSE)
  return(metadata)
  
}

#function calls:
#new_dir=read_multiple_fasta()
#kmer_result=kmer_count(new_dir,metadata)
#metadata=data_preprocessing2()