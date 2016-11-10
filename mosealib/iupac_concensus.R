#information : https://jpwendler.wordpress.com/2013/06/11/an-r-function-to-convert-a-dna-alignment-into-an-iupac-consensus-sequence/
#rscript
library(seqinr)
#########################################################################################
# FUNCTION: IUPAC_consensus
# Takes an aligned fasta object (seqinr) and returns the consensus sequence
#
#USAGE: Rscript <script.R> <aligned_fasta_file.fa>
# 
#########################################################################################
IUPAC_consensus<-function(myfasta, deletionTreatment="ignore" )
{
  myconsensus<-NA
 
  # For each position in the alignment
  for(i in 1:length(myfasta[[1]]))
  {
    # Get the unique nucleotides present at this position
    thisbase<-unique(sapply(myfasta,"[", i))
    
    if(length(thisbase)>1)
    {
    	#check if there is a '-' in sequence
    	if("-" %in% thisbase){
    		#count number of '-' in the column, if it is greater that 50%, replace it with N
	    	count <- table(sapply(myfasta,"[", i))   
	    	count <- as.data.frame(count)
	    	s = subset(count, Var1 == '-')
        total_freq = sum(count$Freq)
	    	percent <- (s$Freq / total_freq) 
	    	if(percent >= 0.4){
	    		thisbase <- gsub('-', 'N', thisbase)
    		}
    	}
    	else{
    		thisbase <- toupper(thisbase)
    		# If the position is variable
   		 }
      
      # Get the IUPAC representation
      thisbase<-base2iupac(thisbase, deletionTreatment=deletionTreatment)
    }
     
    myconsensus<-c(myconsensus, thisbase)
  }
   
  # Get rid of the leading NA
  myconsensus<-myconsensus[-1]
  return(myconsensus)
}
 
 
#########################################################################################
# FUNCTION: base2iupac
# 
#########################################################################################
base2iupac<-function(mybase, deletionTreatment="ignore")
{
  # Just in case this wasn't done before the call
  mybase<-unique(mybase)

  # Check for deletion if more than one base sent
  if(length(mybase)>1)
  {
    delpositions<-which(mybase %in% c(".","-"))
     
    # If there are deletions
    if(length(delpositions)>=1)
    {
      # This will return the IUPAC of the sequences that are there (ignoring the deletion)
      if(deletionTreatment=="ignore")
      {
        mybase<-mybase[-delpositions]
      }else
      {
        mybase<-"-" # This returns a "-", ignoring other bases
      }
    }
  }
   
  # Get in proper format to be recognized in the table
  mybase<-paste(sort(toupper(mybase)), collapse="")
  mytable<-data.frame(base=c("A","C","G","T","AG",
  				"CT","CG","AT","GT",
  				"AC","CGT","AGT","ACT",
  				"ACG","ACGT",".", "-"),
  				iupac=c("A","C","G","T","R",
  				"Y","S","W","K","M","B",
  				"D","H","V","N",".","-"))
 
  matchposition<-which(mytable$base %in% mybase)
   
  # If there's an unambiguous match in the IUPAC table, return it.  
  # If the occurence of '-' was >50% of the time, return an N for that positon
  if(length(matchposition)==1)
  {
    myiupac<-as.character(mytable[matchposition, "iupac"])
  }else
  {
    myiupac<-"N"
  }
 
  return(myiupac)
}
 
 
 
##############################################################
# Let's see it in action!
##############################################################
args <- commandArgs(TRUE)
aligned_fasta <- args[1]

if(aligned_fasta == ''){
  print("USAGE: Rscript <script.R> <aligned_fasta_file.fa>")
  stop()
}

df <- read.fasta(file = aligned_fasta)
#df <- unlist(df)

result <- IUPAC_consensus(df)
result <- toupper(result)
result <- paste(result, collapse = '')
cat(result)
# "N" "c" "W" "D"
 
#
