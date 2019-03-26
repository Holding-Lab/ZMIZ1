#' annotate gene
#
#' This function converts a EntrezID to a Gene name
#' @keywords EntrezID Gene
#' @export
#' @examples annotategene()
#' annotategene_function()
annotategene<-function(x){

    library(biomaRt)
    library(org.Mm.eg.db)
    library(org.Hs.eg.db)

    tab<-AnnotationDbi::select(org.Hs.eg.db, keys=x, columns=c("GENENAME"), keytype="ENTREZID")
    tab<-tab[!duplicated(tab[,"ENTREZID"]),]
    out<-setNames(tab[,"GENENAME"],tab[,"ENTREZID"])
    out<-out[x]
    return(out)
}

#' Symbol 2 Entrez
#
#' This function converts a gene symbol to Entrez ID
#' @keywords EntrezID Gene Symbol
#' @export
#' @examples sym2eg('TFF1')
#' sym2eg_function()

sym2eg<-function(ids){
    library(org.Hs.eg.db)
    list_symbol2eg <- as.character(org.Hs.egALIAS2EG[mappedkeys(org.Hs.egALIAS2EG)])
    ids <- as.character(ids)
    outlist <- list_symbol2eg[ids]
    names(outlist) <- ids
    outlist[is.na(outlist)] <- paste("unknown.", ids[is.na(outlist)], sep = "")
    outlist <- gsub("unknown.unknown.", "", outlist)
    return(outlist)
}


#' Any to Entrez
#
#' This function converts a anyhting to Entrez ID
#' @keywords EntrezID Gene
#' @export
#' @examples any2entrez('TFF1')
#' any2entrez_function()
#'
any2entrez<-function(x){
    library(org.Hs.eg.db)
    tab<-AnnotationDbi::select(org.Hs.eg.db, keys=x, columns=c("ENTREZID"), keytype="ALIAS")
    symbols<-tab[,1]
    entrez<-tab[,2]
    dups<-which(duplicated(symbols))
    if(length(dups)>0){
        symbols<-symbols[-dups]
        entrez<-entrez[-dups]
    }
    out<-setNames(entrez,symbols)
    return(out)
}

#' Ensemble 2 Entrez
#
#' This function converts a Ensemble ID to Entrez ID
#' @keywords Ensemble Entrez Convert
#' @export
#' @examples ens2eg('ENSG00000160182')
#' ens2eg_function()

ens2eg<-function(x){
    library(org.Hs.eg.db)
    if(!exists("ens2egmap")){
        ens2egmap<<-as.list(org.Hs.egENSEMBL2EG)
    }
    out<-ens2egmap[x]
    names(out)<-x
    out2<-sapply(out,function(x){
        if(is.null(x)){
            return(NA)
        } else {
            return(x[1])
        }
    })
    out3<-unlist(out2)
    out<-setNames(out3,names(out))
    return(out)
}

#' Entrez 2 Ensemble
#
#' This function converts a Entrez ID to Ensemble ID
#' @keywords Ensemble Entrez Convert
#' @export
#' @examples ens2eg('7031')
#' ens2eg_function()

### eg2ens function
eg2ens<-function(x){
    library(org.Hs.eg.db)
    if(!exists("eg2ensmap")){
        eg2ensmap<<-as.list(org.Hs.egENSEMBL[mappedkeys(org.Hs.egENSEMBL)])
    }
    out<-eg2ensmap[x]
    names(out)<-x
    out2<-sapply(out,function(x){
        if(is.null(x)){
            return(NA)
        } else {
            return(x[1])
        }
    })
    out3<-unlist(out2)
    out<-setNames(out3,names(out))
    return(out)
}


#' Entrez 2 Symbol
#
#' This function converts a Entrez ID to Gene Symbol
#' @keywords Gene Symbol Entrez Convert
#' @export
#' @examples eg2sym('7031')
#' eg2sym_function()
#'
eg2sym<-function(x){
    library(org.Hs.eg.db)
    if(!exists("eg2symmap")){
        eg2symmap<<-as.list( org.Hs.egSYMBOL[mappedkeys( org.Hs.egSYMBOL)])
    }
    x<-as.character(x)
    out<-eg2symmap[x]
    names(out)<-x
    out2<-sapply(out,function(x){
        if(is.null(x)){
            return(NA)
        } else {
            return(x[1])
        }
    })
    out3<-unlist(out2)
    out<-setNames(out3,names(out))
    return(out)
}



#' Entrez 2 RegSeq
#
#' This function converts a Entrez ID to RefSeq
#' @keywords Gene Symbol Entrez Convert
#' @export
#' @examples eg2refseq('7031')
#' eg2refseq_function()
#'
#'
eg2refseq<-function(x){
    library(org.Hs.eg.db)
    if(!exists("ens2egmap")){
        eg2refseqmap<<-as.list(org.Hs.egREFSEQ)
    }
    out<-eg2refseqmap[x]
    names(out)<-x
    out2<-sapply(out,function(x){
        if(is.null(x)){
            return(NA)
        } else {
            return(x[1])
        }
    })
    out3<-unlist(out2)
    out<-setNames(out3,names(out))
    return(out)
}



