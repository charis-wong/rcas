
#' check if Cas number valid
#' 
#' This function uses the CAS structure see more info on https://www.cas.org/support/documentation/chemical-substances/faqs to check if a string is a valid CAS. 
#' @param cas string to check 
#' @keywords cascheck
#' @return TRUE if string is a valid CAS, FALSE if not 
#' @export


cascheck<-function(cas){
  if(!is.na(cas)){
    len<-str_length(cas)
    checkdigit<-str_sub(cas, len)%>%as.numeric()
    casn<-str_replace_all(cas, "-", "")
    lencasn<-str_length(casn)
    for(i in 1:(lencasn-1)){
      n<-lencasn - i
      producti<-as.numeric(str_sub(casn, n, n))*i 
      if(!exists("checksum")) checksum<-producti else checksum<-producti+checksum
    }
    if(!is.na(checksum)){
      cascheck <- ifelse(checksum %% 10 == checkdigit, TRUE, FALSE)
      rm(checksum)
    }else cascheck <- FALSE
  }else cascheck <- FALSE
  
  
  return(cascheck)
}


#' get CAS number
#' 
#' This function combines different functions to get CAS number for a drug/compound name from DrugBank Vocabulary, Common Chemistry API and PubChem API
#' @param drugname name of drug/compound to parse
#' @export
#' @returns CAS or list of CAS numbers for drug/compound name

getCAS <- function(drugname){
  drugname<-tolower(drugname)
  cas <- getCASFromDB(drugname)
  cas <- unique(cas)
  cas <- cas[!is.na(cas)]
  if(length(cas)==0){
    cas<-getCASFromCC(drugname)
    cas <- unique(cas)
    cas <- cas[!is.na(cas)]
    
    if(length(cas)==0){
      cas<-getCASFromPubChem(drugname)
      cas <- unique(cas)
      cas <- cas[!is.na(cas)]
    }
    
    if(length(cas)==0){
      cas<-getCASFromPubChem2(drugname)
      cas <- unique(cas)
      cas <- cas[!is.na(cas)]
    }
    
    
  }
  if(length(cas)==0){cas<-NA}
  
  return(cas)
  
}

#' get CAS number from sources other than drugbank - common chemistry API and pubchem API
#' @param drugname drug/compound name to parse
#' @return cas or list of cas numbers for drug/compound name
#' @export

getCASNonDB <- function(drugname){
  drugname<-tolower(drugname)
  cas<-getCASFromCC(drugname)
  cas <- unique(cas)
  cas <- cas[!is.na(cas)]
  
  if(length(cas)==0){
    cas<-getCASFromPubChem(drugname)
    cas <- unique(cas)
    cas <- cas[!is.na(cas)]
  }
  
  if(length(cas)==0){
    cas<-getCASFromPubChem2(drugname)
    cas <- unique(cas)
    cas <- cas[!is.na(cas)]
  }
  
  if(length(cas)==0){cas<-NA}
  
  return(cas)
  
}

#' get CAS number from drugbank
#' @param drugname drug/compound name to parse
#' @return cas or list of cas numbers for drug/compound name
#' @export


getCASFromDB <- function(drugname){
  drugname<-tolower(drugname)
  data(drugBank)
    drugBankData<-drugBankVocabulary%>%filter(drugname == tolower(`Common name`))
    if(nrow(drugBankData)>0){
      CAS<-drugBankData$CAS
    }else{
      
      
      drugBankData<-drugBankVocabulary%>%filter(grepl(drugname,tolower(`Common name`)))
      if(nrow(drugBankData)>0){
        CAS <- drugBankData$CAS
      } else {
        drugBankData<-drugBankVocabulary %>% filter(grepl(drugname, tolower(Synonyms)))
        
        if(nrow(drugBankData)>0){
          CAS<-drugBankData$CAS
        } else CAS<- NA
      }}
  return(CAS)
}

#' get CAS number from Common Chemistry API
#' @param drugname drug/compound name to parse
#' @return cas or list of cas numbers for drug/compound name
#' @export

getCASFromCC<-function(drug){
  base<- "commonchemistry.cas.org/api"
  drugname<- tolower(drug)
  drugname <-gsub(" ", "%20", drugname)
  call1<-paste0(base, "/search?q=", drugname, "&size=1")
  call1<-URLencode(call1)
  get_CAS<-GET(call1)
  
  get_CAS_text<-content(get_CAS, "text")
  
  get_CAS_json <- fromJSON(get_CAS_text, flatten = "TRUE")
  
  CAS<-get_CAS_json$results
  CAS<-ifelse(length(CAS) == 0, NA, CAS$rn)
  return(CAS)
}

#' Extract CAS number from string of drug name 
#' @param drugname drug/compound name to parse
#' @return cas or list of cas numbers for drug/compound name
#' @export

extractCASFromName <- function(drugname){
  hasCAS <- str_detect(drugname, "\\d+-\\d\\d-\\d")
  if(hasCAS == TRUE){
    CAS <- str_extract(drugname, "\\d+-\\d\\d-\\d")
  } else CAS <- NA
  return(CAS)
}

#' get CAS number from Pubchem where listed on pubchem as compound
#' @param drugname drug/compound name to parse
#' @return cas or list of cas numbers for drug/compound name
#' @export

getCASFromPubChem<-function(drug){
  drugname<-tolower(drug)
  drugname <-gsub(" ", "", drugname)
  drugname <-gsub("\\?", "*", drugname)
  
  pubchemurl<- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", drugname, "/synonyms/XML")
  pubchemxmldata<-GET(pubchemurl)
  if (pubchemxmldata$status_code ==200) { 
    try({data<-read_xml(pubchemurl)
    synonyms<-xml_children(xml_children(data))%>%xml_text()
    CAS<-synonyms[str_detect(synonyms, "^\\d+-\\d\\d-\\d$")]
    CAS <- unique(CAS)
    for(i in CAS){
      if(cascheck(i) == TRUE) {
        if(!exists("CASlist")) CASlist <- i else CASlist <- append(CASlist, i)
      }
    }
    })
  } 
  
  
  if(!exists("CASlist")){
    drugname <- gsub("\\-", "%20", drug)
    drugname <- gsub(" ", "%20", drugname)
    pubchemurl<- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", drugname, "/synonyms/XML")
    pubchemxmldata<-GET(pubchemurl)
    if (pubchemxmldata$status_code ==200) { 
      try({data<-read_xml(pubchemurl)
      synonyms<-xml_children(xml_children(data))%>%xml_text()
      CAS<-synonyms[str_detect(synonyms, "^\\d+-\\d\\d-\\d$")]
      CAS <- unique(CAS)
      for(i in CAS){
        if(cascheck(i) == TRUE) {
          if(!exists("CASlist")) CASlist <- i else CASlist <- append(CASlist, i)
        }
      }
      })
    }
    
    
  }
  
  if(!exists("CASlist")){
    drugname <- gsub(" ", "%20", drug)
    pubchemurl<- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", drugname, "/synonyms/XML")
    pubchemxmldata<-GET(pubchemurl)
    if (pubchemxmldata$status_code ==200) { 
      try({data<-read_xml(pubchemurl)
      synonyms<-xml_children(xml_children(data))%>%xml_text()
      CAS<-synonyms[str_detect(synonyms, "^\\d+-\\d\\d-\\d$")]
      CAS <- unique(CAS)
      for(i in CAS){
        if(cascheck(i) == TRUE) {
          if(!exists("CASlist")) CASlist <- i else CASlist <- append(CASlist, i)
        }
      }
      })
    } 
    
    
  }
  
  if(exists("CASlist")){
    CASlist <- CASlist[!is.na(CASlist)] %>%unique()
    if (length(CASlist) == 0) CASlist <- NA
  } else {CASlist <- NA}
  
  
  return(CASlist)
}


#' get CAS number from Pubchem where listed on pubchem as substance
#' @param drugname drug/compound name to parse
#' @return cas or list of cas numbers for drug/compound name
#' @export

getCASFromPubChem2<-function(drug){
  drugname<-tolower(drug)
  drugname <-gsub(" ", "", drugname)
  drugname <-gsub("\\?", "*", drugname)
  
  pubchemurl<- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/", drugname, "/synonyms/XML")
  pubchemxmldata<-GET(pubchemurl)
  if (pubchemxmldata$status_code ==200) { 
    try({data<-read_xml(pubchemurl)
    synonyms<-xml_children(xml_children(data))%>%xml_text()
    CAS<-synonyms[str_detect(synonyms, "^\\d+-\\d\\d-\\d$")]
    CAS <- unique(CAS)
    for(i in CAS){
      if(cascheck(i) == TRUE) {
        if(!exists("CASlist")) CASlist <- i else CASlist <- append(CASlist, i)
      }
    }
    })
  } 
  
  if(!exists("CASlist")){
    drugname <- gsub("\\-", "%20", drug)
    drugname <- gsub(" ", "%20", drugname)
    pubchemurl<- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/", drugname, "/synonyms/XML")
    pubchemxmldata<-GET(pubchemurl)
    if (pubchemxmldata$status_code ==200) { 
      try({data<-read_xml(pubchemurl)
      synonyms<-xml_children(xml_children(data))%>%xml_text()
      CAS<-synonyms[str_detect(synonyms, "^\\d+-\\d\\d-\\d$")]
      CAS <- unique(CAS)
      for(i in CAS){
        if(cascheck(i) == TRUE) {
          if(!exists("CASlist")) CASlist <- i else CASlist <- append(CASlist, i)
        }
      }
      })
    } else CASlist <- NA
    
    
  }
  
  if(exists("CASlist")){
    CASlist <- CASlist[!is.na(CASlist)] %>%unique()
    if (length(CASlist) == 0) CASlist <- NA
  } else {CASlist <- NA}
  
  
  return(CASlist)
}