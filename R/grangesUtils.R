# general useful functions for granges

#' Convert DNAStringSet to GenomicRanges
#' @param DNAss A DNAStringSetGR
#' @param defaultStrand One of "+", "-" or "*" for the strand from which the sequences come from (default="*")
#' @return A GenomicRanges object
#' @examples
#' DNAstringsetToGR(Biostrings::DNAStringSet(x=c(chr1="agctgtagct",chr2="agagagagttt"))
#' @export
DNAStringSetToGR<-function(DNAss,defaultStrand="*") {
  #extract only the first word in fasta '>' row to be the seq name
  chrNames<-sapply(strsplit(names(DNAss)," "),'[[',1)
  gr<-GenomicRanges::GRanges(seqnames=chrNames,
                             ranges=IRanges::IRanges(start=1,width=Biostrings::width(DNAss)),
                             strand=defaultStrand)
  return(gr)
}



#' Convert BSgenome to DNAStringSet
#' @param BSgnm A BSgenome object
#' @return A DNAStringSet object
#' @examples
#' BSgenomeToDNAStringSet(BSgenome.Celegans.UCSC.ce11::Celegans)
#' @export
BSgenomeToDNAStringSet<-function(BSgnm){
  listOfChr<-sapply(seqnames(BSgnm),function(x){BSgnm[[x]]})
  DNAss<-methods::as(listOfChr,"DNAStringSet")
  return(DNAss)
}



#' Convert BSgenome to GenomicRanges
#' @param BSgnm A BSgenome object
#' @param defaultStrand One of "+", "-" or "*" for the strand from which the sequences come from (default="*")
#' @return A GenomicRanges object
#' @examples
#' BSgenomeToGR(BSgenome.Celegans.UCSC.ce11::Celegans)
#' @export
BSgenomeToGR<-function(BSgnm,defaultStrand="*") {
  gr<-GenomicRanges::GRanges(seqnames=BSgenome::seqnames(BSgnm),
                             ranges=IRanges::IRanges(start=1,end=GenomeInfoDb::seqlengths(BSgnm)),
                             strand="*")
  return(gr)
}



#' Convert UCSC to wormbase chromosome names
#' @param ucscGR A GRanges object with ucsc chomosome names ("chrI"..."chrM")
#' @return A GenomicRanges object
#' @examples
#' ucscToWbGR(BSgenomeToGR(BSgenome.Celegans.UCSC.ce11::Celegans))
#' @export
ucscToWbGR<-function(ucscGR) {
  wbGR<-ucscGR
  GenomeInfoDb::seqlevels(wbGR)<-gsub("chr","",GenomeInfoDb::seqlevels(wbGR))
  GenomeInfoDb::seqlevels(wbGR)<-gsub("M","MtDNA",GenomeInfoDb::seqlevels(wbGR))
  if(class(ucscGR)=="BSgenome") {
    wbGR@provider<-"Wormbase"
    wbGR@provider_version<-"WS235"
  }
  return(wbGR)
}



#' Convert wormbase to ucsc chromosome names
#' @param wbGR A GRanges object with ucsc chomosome names ("chrI"..."chrM")
#' @return A GenomicRanges object
#' @examples
#' wbToUcscGR(ucscToWbGR(BSgenomeToGR(BSgenome.Celegans.UCSC.ce11::Celegans)))
#' @export
wbToUcscGR<-function(wbGR) {
  ucscGR<-wbGR
  GenomeInfoDb::seqlevels(ucscGR)<-gsub("MtDNA","M",GenomeInfoDb::seqlevels(ucscGR))
  GenomeInfoDb::seqlevels(ucscGR)<-paste0("chr",GenomeInfoDb::seqlevels(ucscGR))
  return(ucscGR)
}




#' Convert MIndex object to Genomic Ranges object
#'
#' Function to convert Mindex object (obtained by matching patterns on DNAstringset) to genomic ranges
#'
#' @param mIdx A MIndex object
#' @return A GenomicRanges object
#' @examples
#' mIdxToGR(Biostrings::vmatchPattern("ATTTAGGGTTTTAGAATACTGCCATTAATTAAAAAT",ttTi5605dna))
#' @export
mIdxToGR<-function(mIdx) {
  allGR<-GenomicRanges::GRanges()
  GenomeInfoDb::seqlevels(allGR)<-names(mIdx)
  for (n in names(mIdx)) {
    if (length(mIdx[[n]])>0) {
      gr<-GenomicRanges::GRanges(seqnames=S4Vectors::Rle(c(n),length(mIdx[[n]])),
                     ranges=mIdx[[n]], strand="*")
      allGR<-append(allGR,gr)
    }
  }
  return(allGR)
}



#' Check if vector is only NAs
#'
#' @param vec Vector of values to be checked
#' @examples
#' allNAs(c(NA, NA, NA))
#' allNAs(c(NA,1,NA))
#' allNAs(c(1,2,3))
#' allNAs(c())
#' @return boolean TRUE or FALSE
allNAs<-function(vec) {
  if (sum(is.na(vec)) == length(vec))  {
    returnVal=TRUE
  } else {
    returnVal=FALSE
  }
  return(returnVal)
}



#' Apply a function on a data from gr2 with windows provided by gr1
#'
#' Apply a function to summarise data from multiple metadata columns of GRanges in gr2 that
#' fall within GRanges provided by gr1
#' @param gr1 A GenomicRanges object with windows of interest
#' @param gr2 A GenomicRanges object with data of interest
#' @param applyTo Vector with names of columns in gr2 on which to apply the function
#' @param fun Function to apply (e.g. sum, mean, paste0)
#' @return GRanges from gr1 with summed values of metadata columns from gr2
#' @examples
#' # single gr
#' gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1,3,7), c(5,6,10),names=paste0("win", letters[1:3])), score=4:6)
#' gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 3, 8), c(1, 3, 8),names=paste0("dataID:", letters[1:3])), score=c(10,20,30))
#' applyGRonGR(gr1,gr2,applyTo="score",fun=sum)
#' gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1,3,7), c(5,6,10),names=paste0("win", letters[1:3])), score=4:6, other=1:3)
#' gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 3, 8), c(1, 3, 8),names=paste0("dataID:", letters[1:3])), score=c(10,20,30), other=c(100,NA,300))
#' applyGRonGR(gr1,gr2,c("score","other"),fun=sum,na.rm=T)
#' @export
applyGRonGR<-function(gr1,gr2,applyTo,fun,...) {
  newGR<-IRanges::subsetByOverlaps(gr1,gr2)
  ol<-IRanges::findOverlaps(gr1,gr2)
  grps<-S4Vectors::queryHits(ol)
  dataCols<-tibble::as_tibble(cbind(S4Vectors::DataFrame(grps),
                                    tibble::as_tibble(GenomicRanges::mcols(gr2)[S4Vectors::subjectHits(ol),applyTo])))
  colnames(dataCols)<-c("grps",applyTo)
  newData<-dataCols %>% dplyr::group_by(grps) %>% dplyr::summarise_all(fun,na.rm=T) %>% dplyr::select(applyTo)
  # mark groups where all values are NA
  NAidx<-dataCols %>% dplyr::group_by(grps) %>% dplyr::summarise_all(allNAs) %>% dplyr::select(applyTo)
  newData[as.matrix(NAidx)]<-NA
  # either replace columns or add columns
  if (sum(!(applyTo %in% colnames(GenomicRanges::mcols(newGR))))==0) {
    GenomicRanges::mcols(newGR)[,applyTo]<-newData[,applyTo]
  } else {
    GenomicRanges::mcols(newGR)<-cbind(GenomicRanges::mcols(newGR),newData)
  }
  return(newGR)
}




