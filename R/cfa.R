# S.F. 15-01-2000

# Analysis of configuration frequencies 
# (Konfigurationsfrequenzanalyse nach Lienert und Krauth, KFA)

# Accepts more than two values in each cofiguration. Input may be
# factors or numerical data. 
# Only existing configurations are taken into account.
# Checks for missing data and will ignore cases containig NA's in any
# variable including the counts if ignore.na is TRUE



# function name: cfa

# arguments:
# configmatrix: matrix of factors or numbers
# cntvector   : vector of counts
# descending=TRUE	output in descending order
# sort.on.chisq=TRUE	sort output on chi squared
# sort.on.n=FALSE	sort output on n
# ignore.na=FALSE	ignore NA in input
# binom.test=FALSE	perform binomial test
# binom.test.limit=10   but only if cell contents < this limit
# lehmacher=F		perform Lehmacher's test
# holm.alpha=0.01	alpha for Holm's adjustment
# verbose=FALSE        long output

# returns:
# List consisting of table and list of global parameters


# helper function for binomial coefficients

bin.coeff<-function(n,k) exp(lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1))

# helper function for binomial test

bin.test<-function(ni,ep,ntotal) # inexact for large numbers (>10) !
   {
     sum.p<-0
     for (j in ni:ntotal)
       {
         sum.p<-sum.p+((ep^j)*(1.0-ep)^(ntotal-j)*bin.coeff(ntotal,j) )
       }
     sum.p
   }

# helper function for z-test

z.test<-function(ni,ep,ntotal)
   {
     num<-(ni-(ntotal*ep))
     den<-sqrt(ntotal*ep*(1.0-ep))
     if ((ni*ep)<=10.0) num<-num-0.5 # continuity correction
     num/den
   }

# main function

cfa<-function(configmatrix,cntvector,
              descending=TRUE,       # output in descending order
              sort.on.chisq=TRUE,    # sort output on chi squared
              sort.on.n=FALSE,       # sort output on n
              ignore.na=FALSE,       # ignore NA in input
              binom.test=FALSE,      # perform binomial test
              binom.test.limit=10,   # (but only if cell contents < this limit)
              bonferroni.p.z=T,      # Bonferroni-adjust sig(p(z))
              bonferroni.alpha=0.05, # Alpha for Bonferroni-adjustment
              lehmacher=F,           # perform Lehmacher's test
              holm.alpha=0.01,       # Alpha for Holm's adjustment
              verbose=FALSE)         # long output
  {
    cols<-ncol(configmatrix)
    configs<-cols
    if (ignore.na)
      {
       
        delvec<-vector(mode="numeric",length=0)
        for (i in 1:nrow(configmatrix))
          {
            cfgna<-any(is.na(configmatrix[i,1:cols]))
            cntna<-is.na(cntvector[i,1])
            if (cfgna || cntna) delvec<-c(-i,delvec)            
          }
        if (length(delvec>0)) # important!
              {
                interim<-configmatrix[delvec,1:cols]
                configmatrix<-interim
                interim<-cntvector[delvec,1,drop=T]
                cntvector<-as.data.frame(interim)
              }
       }
    
     if (any(is.na(configmatrix)))
         stop("Missing data in configs!")
     if (any(is.na(cntvector)))
         stop("Missing data in counts!")  

    for (i in 1:cols) 
          {              
            if (!is.factor(configmatrix[[i]]) ) # coerce to factors if necessary 
               configmatrix[[i]]<-factor(configmatrix[[i]])
            levels(configmatrix[[i]])<-format(levels(configmatrix[[i]]))

          }
    
    bivariate<-FALSE
    maxfac<-1
    factorn<-numeric(length=configs)
    for (i in 1:ncol(configmatrix))
      {
        locmax<-max(codes(configmatrix[,i]),na.rm=T)
        factorn[i]<-locmax
        if (locmax>maxfac) maxfac<-locmax
      }
    if (maxfac==2) bivariate<-TRUE
    if (configs<2) stop("Only one configuration specified")
    if (bivariate) df<-2^configs-configs-1
    n<-sum(cntvector)
    colsums<-tapply(unlist(cntvector),configmatrix[,1],sum)
    for (i in 2:ncol(configmatrix))
         {
              colsums<-c(colsums,tapply(unlist(cntvector),configmatrix[,i],sum))
         }    
    colsummatrix<-matrix(data=NA,nrow=configs,ncol=maxfac)
    startpos<-0
    dfprod<-1
    for (i in 1:ncol(configmatrix))
         {
            ncodes<-max(codes(configmatrix[,i]))
            idx<-startpos
            for (j in 1:ncodes)
              {
                colsummatrix[i,j]<-colsums[idx+j]
              }
            startpos<-startpos+ncodes
         }

    cfgsum<-vector(mode="numeric",length=0)
    expected<-vector(mode="numeric",length=0)
    chisq<-vector(mode="numeric",length=0)
    zleh<-vector(mode="numeric",length=0)
    q<-vector(mode="numeric",length=0)
    p.bin<-vector(mode="numeric",length=0)
    z<-vector(mode="numeric",length=0)
    p.z<-vector(mode="numeric",length=0)
    
    for (i in 1:nrow(configmatrix))
         {
           idx<-""
           prodsums<-1
           lprodsum1<-1
           lprodsum2<-1
           for (j in 1:ncol(configmatrix))
             {
               cidx<-(configmatrix[i,j])
               colidx<-codes(cidx)
               lprodsum1<-lprodsum1*(colsummatrix[j,colidx]-1)
               lprodsum2<-lprodsum2*(colsummatrix[j,colidx])
               prodsums<-prodsums*colsummatrix[j,colidx]
               idx<-paste(idx,as.character(configmatrix[i,j]),sep=" ")
             }
           plehmacher1<-lprodsum1/(n-1)^configs
           plehmacher2<-lprodsum2/n^configs
           vlehmacher<-n*plehmacher2*(1.0-plehmacher2-(n-1)*(plehmacher2-plehmacher1 ))
           if(is.na(cfgsum[idx])) cfgsum[idx]<-0
             cfgsum[idx]<-cfgsum[idx]+cntvector[i,1] 
           expected[idx]<-prodsums/n^(configs-1)
#          Lehmacher's test
           zlehmacher<-(cfgsum[idx]-expected[idx])/sqrt(vlehmacher)
#          Lehmacher's test with continuity correction
           if (cfgsum[idx]>expected[idx])
             zlehmacher<-(cfgsum[idx]-expected[idx]-0.5)/sqrt(vlehmacher)
           if (cfgsum[idx]<expected[idx])
             zlehmacher<-(cfgsum[idx]-expected[idx]+0.5)/sqrt(vlehmacher)
           zleh[idx]<-zlehmacher
           chisq[idx]<-((cfgsum[idx]-expected[idx])^2)/expected[idx]
           q[idx]<-abs(cfgsum[idx]-expected[idx])/max(expected[idx],n-expected[idx])
           if (bivariate)
              {
                if (binom.test==TRUE)
                  {
                    if (cfgsum[idx]<=binom.test.limit )
                       p.bin[idx]<-bin.test(cfgsum[idx],expected[idx]/n,n)
                    if (cfgsum[idx]>binom.test.limit)
                       p.bin[idx]<-NA
                  }
                z[idx]<-z.test(cfgsum[idx],expected[idx]/n,n)
                p.z[idx]<-1.0-pnorm(z[idx])
              }
           if (bivariate==FALSE)
             {
                z[idx]<-sqrt(chisq[idx])
                p.z[idx]<-1.0-pnorm(z[idx])
             }
         }

     if (bivariate==TRUE)
      {
        df<-(2^configs-configs-1)    # (KL27&57)
        dfcell<-1                    # (KL 26)
      }
    if (bivariate==FALSE)
     {
       df<-prod(factorn)+sum(-factorn)+(configs-1)  # (LW96)
       dfcell<-prod(factorn-1)                      # richtig?
     }
    nz<-length(z)
    if (bonferroni.p.z==FALSE) nz<-1 # no adjustment
#   pzsig is Bonferroni-adjusted
    pzsig<-((1.0-pnorm(z,0,1))<(bonferroni.alpha/nz) | 
                 (pnorm(z,0,1)<(bonferroni.alpha/nz)))
    rankedzleh<-rank(zleh)
    names(rankedzleh)<-names(zleh)
    alpha.holm<-holm.alpha/rankedzleh # Holm's alpha adjustment 
    lehsig<-((1-pnorm(zleh,0,1)<alpha.holm)|(pnorm(zleh,0,1)<alpha.holm))
    totalxsq<-sum(chisq)
    p<-1.0-pchisq(totalxsq,df)
    sortidx<-seq(1:length(cfgsum))
    if (sort.on.chisq==T)
      sortidx<-order(chisq)
    if (sort.on.n==T)
      sortidx<-order(cfgsum)             
    if( descending==T) sortidx<-rev(sortidx)
    sorted.chisq<-chisq[sortidx]
    sorted.cfgsum<-cfgsum[sortidx]
    sorted.expected<-expected[sortidx]
    sorted.q<-q[sortidx]
    sorted.p.bin<-p.bin[sortidx]
    sorted.z<-z[sortidx]
    sorted.p.z<-p.z[sortidx]
    sorted.pzsig<-pzsig[sortidx]
    sorted.zleh<-zleh[sortidx]
    sorted.lehsig<-lehsig[sortidx]
    restbl<-cbind(sorted.cfgsum,sorted.expected,sorted.chisq)
    dimnames(restbl)[[2]]<-c("n","expected","chisq")
    if (verbose==T)
      {
       restbl<-cbind(sorted.cfgsum,sorted.cfgsum/n*100.0,sorted.expected,
                     sorted.q,sorted.chisq,sorted.z,sorted.p.z)
       dimnames(restbl)[[2]]<-c("n","pct.","expected","Q","chisq.","z","p(z)")
      }
    if(verbose==T)
       {
         restbl<-cbind(restbl,sorted.pzsig)
         colnames(restbl)[ncol(restbl)]<-"sig(p(z))"
       }
    if (bivariate==TRUE & binom.test==TRUE)
       {
         restbl<-cbind(restbl,sorted.p.bin)
         colnames(restbl)[ncol(restbl)]<-"p(bin)"
       }
    if (lehmacher==TRUE)
       {
         restbl<-cbind(restbl,sorted.zleh)
         colnames(restbl)[ncol(restbl)]<-"z(Lehmacher)"
         restbl<-cbind(restbl,sorted.lehsig)
         colnames(restbl)[ncol(restbl)]<-"Lehmacher.sig"
       }
    if (verbose==FALSE)
        res<-list(table=restbl,parms=c(sum.of.cells=n,overall.chisq=totalxsq,overall.p=p,df=df))
    if (verbose==TRUE)
        res<-list(table=restbl,parms=c(sum.of.cells=n,overall.chisq=totalxsq,overall.p=p,df=df,
             sig.limit=0.05/2^configs))
    class(res)<-"cfa"
    res
             
  }


# Bare bones version of the above function which returns only
# the global chi squared (used for hierarchical cfa)

shortcfa<-function(configmatrix,cntvector)                      
  {
    cols<-ncol(configmatrix)
    configs<-cols
    ignore.na<-TRUE
    if (ignore.na)
      {
       
        delvec<-vector(mode="numeric",length=0)
        for (i in 1:nrow(configmatrix))
          {
            cfgna<-any(is.na(configmatrix[i,1:cols]))
            cntna<-is.na(cntvector[i,1])
            if (cfgna || cntna) delvec<-c(-i,delvec)            
          }
        if (length(delvec>0)) # important!
              {
                interim<-configmatrix[delvec,1:cols]
                configmatrix<-interim
                interim<-cntvector[delvec,1,drop=T]
                cntvector<-as.data.frame(interim)
              }
       }
    
     if (any(is.na(configmatrix)))
         stop("Missing data in configs!")
     if (any(is.na(cntvector)))
         stop("Missing data in counts!")  

    maxfac<-1
    locmax<-1
    for (i in 1:cols) 
          {              
            if (!is.factor(configmatrix[[i]]) ) # coerce to factors if necessary
               configmatrix[[i]]<-factor(configmatrix[[i]]) 
            locmax<-max(codes(configmatrix[,i]),na.rm=T)
            if (locmax>maxfac) maxfac<-locmax
           
          }
    if (configs<2) stop("Only one configuration specified")
    n<-sum(cntvector)
    colsums<-tapply(unlist(cntvector),configmatrix[,1],sum)
    for (i in 2:ncol(configmatrix))
         {
              colsums<-c(colsums,tapply(unlist(cntvector),configmatrix[,i],sum))
         }    
    colsummatrix<-matrix(data=NA,nrow=configs,ncol=maxfac)
    startpos<-0
    dfprod<-1
    for (i in 1:ncol(configmatrix))
         {
            ncodes<-max(codes(configmatrix[,i]))
            idx<-startpos
            for (j in 1:ncodes)
              {
                colsummatrix[i,j]<-colsums[idx+j]
              }
            startpos<-startpos+ncodes
         }

    cfgsum<-vector(mode="numeric",length=0)
    expected<-vector(mode="numeric",length=0)
    chisq<-vector(mode="numeric",length=0)  
    for (i in 1:nrow(configmatrix))
         {
           idx<-""
           prodsums<-1
           for (j in 1:ncol(configmatrix))
             {
               cidx<-(configmatrix[i,j])
               colidx<-codes(cidx)
               prodsums<-prodsums*colsummatrix[j,colidx]
               idx<-paste(idx,as.character(configmatrix[i,j]),sep=" ")
             }
           if(is.na(cfgsum[idx])) cfgsum[idx]<-0
             cfgsum[idx]<-cfgsum[idx]+cntvector[i,1] 
           expected[idx]<-prodsums/n^(configs-1)
           chisq[idx]<-((cfgsum[idx]-expected[idx])^2)/expected[idx]
          
         }
    totalxsq<-sum(chisq)
    totalxsq           
  }


# Extension of cfa to two or more samples (counts)

# function name: mcfa

# arguments:
# configmatrix: matrix of factors or numbers
# cntvector   : vector of counts
# descending=TRUE	output in descending order
# sort.on.chisq=TRUE	sort output on chi squared
# ignore.na=FALSE	ignore NA in input
# verbose=FALSE        long output

# returns:
# List consisting of table and list of global parameters


mcfa<-function(configmatrix,cntmatrix,
              descending=TRUE,
              sortonchisq=TRUE,
              ignore.na=FALSE,
              verbose=FALSE)
  {
    if(ncol(cntmatrix)<2) stop("Only one column of counts - no multiple cfa")
    cols<-ncol(configmatrix)
    configs<-cols
    counts<-ncol(cntmatrix)
    if (ignore.na)
      {
       
        delvec<-vector(mode="numeric",length=0)
        for (i in 1:nrow(configmatrix))
          {
            cfgna<-any(is.na(configmatrix[i,1:cols]))
            cntna<-any(is.na(cntmatrix[i,1:counts]))
            if (cfgna || cntna) delvec<-c(-i,delvec)            
          }
        if (length(delvec>0)) # important!
              {
                interim<-configmatrix[delvec,1:cols]
                configmatrix<-interim
                interim<-cntmatrix[delvec,1:counts]
                cntmatrix<-as.data.frame(interim)
              }
       }
     cases<-nrow(configmatrix)
     ntotal<-sum(cntmatrix)
     if (any(is.na(configmatrix)))
         stop("Missing data in configs!")
     if (any(is.na(cntmatrix)))
         stop("Missing data in counts!")  

    for (i in 1:cols) 
          {              
            if (!is.factor(configmatrix[[i]]) ) # coerce to factors if necessary 
            configmatrix[[i]]<-factor(configmatrix[[i]])
            levels(configmatrix[[i]])<-format(levels(configmatrix[[i]]))

          }
    
    bivariate<-FALSE
    maxfac<-1
    factorn<-numeric(length=configs)
    for (i in 1:ncol(configmatrix))
      {
        locmax<-max(codes(configmatrix[,i]),na.rm=T)
        factorn[i]<-locmax
        if (locmax>maxfac) maxfac<-locmax
      }
    if (maxfac==2) bivariate<-TRUE
    if (configs<2) stop("Only one configuration specified")
    n<-apply(cntmatrix,2,sum)
    colsumfield<-array(NA,dim=c(configs,maxfac,counts))

    for (sample in 1:counts)
      {
        colsums<-tapply(unlist(cntmatrix[sample]),configmatrix[,1],sum)
        for (i in 2:ncol(configmatrix))
         {
              colsums<-c(colsums,tapply(unlist(cntmatrix[sample]),configmatrix[,i],sum))
         }    
        colsummatrix<-matrix(data=NA,nrow=configs,ncol=maxfac)
        startpos<-0
        dfprod<-1
        for (i in 1:ncol(configmatrix))
         {
            ncodes<-max(codes(configmatrix[,i]))
            dfprod<-dfprod*ncodes
            idx<-startpos
            for (j in 1:ncodes)
              {
                colsummatrix[i,j]<-colsums[idx+j]
              }
            startpos<-startpos+ncodes
         }
        colsumfield[,,sample]<-colsummatrix
      }

    cfgsum<-vector(mode="numeric",length=0)
    expected<-vector(mode="numeric",length=0)
    sumchisq<-vector(mode="numeric",length=0)
    p.chisq<-vector(mode="numeric",length=0)

    for (sample in 1:counts)
     {
       cfgsum<-vector(mode="numeric",length=0)
                     
       for (i in 1:nrow(configmatrix))
         {
           idx<-""
           for (j in 1:ncol(configmatrix))
             {
               idx<-paste(idx,as.character(configmatrix[i,j]),sep=" ")
             }
           if(is.na(cfgsum[idx])) cfgsum[idx]<-0
           cfgsum[idx]<-cfgsum[idx]+cntmatrix[i,sample] 
         }
       if (sample==1)
         {
            cfgs<-length(cfgsum)
            allcfgsum<-matrix(data=NA,nrow=cfgs,ncol=counts)
            rownames(allcfgsum)<-names(cfgsum)
            allexpected<-matrix(data=NA,nrow=cfgs,ncol=counts)
            rownames(allexpected)<-names(expected)
         }     
       allcfgsum[,sample]<-cfgsum
     }

    cs<-apply(allcfgsum,2,sum)
    rs<-apply(allcfgsum,1,sum)
    for (i in 1:length(rs))
     {
       for (j in 1:length(cs))
         {
            allexpected[i,j]<-rs[i]*cs[j]/ntotal
         }
     }

    if (bivariate==TRUE)
      {
        df<-(2^configs-1)*(counts-1) # (KL82)
        dfcell<-1                    # (KL77)
      }
    if (bivariate==FALSE)
     {
       df<-(prod(factorn)-1)*(configs-1)  # (rows-1)*(cols-1)
       dfcell<-1                          # according to Lautsch
     }                                                 

    chisqsums<-apply((allcfgsum-allexpected)^2/allexpected,1,sum)
    totalxsq<-sum(chisqsums)
    p<-1.0-pchisq(totalxsq,df)

    p.chisq<-pchisq(chisqsums,dfcell)   

    sortidx<-seq(1:length(chisqsums))
    if (sortonchisq==T)
      sortidx<-order(chisqsums)       
    if( descending==T) sortidx<-rev(sortidx)
    sorted.chisq<-chisqsums[sortidx]

    sorted.cfgsum<-allcfgsum[sortidx,]
    sorted.expected<-allexpected[sortidx,]
    sorted.p.chisq<-1.0-p.chisq[sortidx]  
    
    restbl<-cbind(sorted.cfgsum,sorted.expected,sorted.chisq)
    dimnames(restbl)[[2]]<-c(paste("n",1:counts,sep=""),
                             paste("exp.",1:counts,sep=""),
                             "chi.sq")
    if (verbose==T)
      {
         restbl<-cbind(restbl,sorted.p.chisq)
         colnames(restbl)[ncol(restbl)]<-"p(chisq)"
      }
    resparms<-list(total.chi.sq=totalxsq,df=df,p=p)
    reslist<-list(table=restbl,parms=resparms)
    class(reslist)<-"mcfa"             
    reslist
  }


print.cfa<-function(x)
   {
     cat("\n*** Analysis of configuration frequencies ***\n\n")
     print(x$table)
     cat("\n")
     cat("Overall chi squared ",x$parms["overall.chisq"],"\n")
     cat("p(chi squared)      ",x$parms["overall.p"],"\n")
     cat("Degrees of freedom  ",x$parms["df"],"\n")
     cat("Total n             ",x$parms["sum.of.cells"],"\n")
   }

print.mcfa<-function(x)
   {
     cat("\n*** n-sample analysis of configuration frequencies ***\n\n")
     print(x$table)
     cat("\n")
     cat("Overall chi squared ",x$parms[["total.chi.sq"]],"\n")
     cat("p(chi squared)      ",x$parms[["p"]],"\n")
     cat("Degrees of freedom  ",x$parms[["df"]],"\n")
   }

callcfa<-function(configs,counts,...)
  {
     ncfg<-ncol(configs)
     ncnt<-ncol(counts)
     lcfg<-nrow(configs)
     lcnt<-nrow(counts)
     
     if (!is.data.frame(configs))
       stop("Configs must be a data frame")
     if (!is.data.frame(counts)) 
       stop("Counts must be a data frame")   
     if (lcfg!=lcnt)
       stop("Configs and counts must have the same number of rows")
     if (ncfg<2)
       stop("Cannot calculate a CFA with only one configuration")

     if (ncnt==1)
        res<-cfa(configs,counts,...)
     if (ncnt>1)
        res<-mcfa(configs,counts,...)
     res
  }

# delete subconfigurations and run CFA for the remaining combination
perm.cfa<-function(configmatrix,cntvector)
   {    
     for (i in 1:ncol(configmatrix))
      {
         configmatrix.1<-configmatrix[-i]
         res<-shortcfa(configmatrix.1,cntvector)
         cfglbl<-paste(names(configmatrix.1),sep=" ")
         h.cfa.configs<<-c(h.cfa.configs,list(cfglbl))
         h.cfa.chisqs<<-c(h.cfa.chisqs,res)
         if (ncol(configmatrix.1)>=3)
           {
             perm.cfa(configmatrix.1,cntvector)             
           }
      }
     list(h.cfa.configs,h.cfa.chisqs)
   }

# Main function for hierarchical CFA

# function name: hier.cfa

# arguments:
# configmatrix: matrix of factors or numbers
# cntvector   : vector of counts
#
# returns:
# list, consisting of the chisq's and the order of the table

hier.cfa<-function(configmatrix,cntvector)
   {
     lcfg<-nrow(configmatrix) 
     lcnt<-nrow(cntvector)
     if (!is.data.frame(configmatrix))
       stop("Configs must be a data frame")
     if (!is.data.frame(cntvector)) 
       stop("Counts must be a data frame")   
     if (lcfg!=lcnt)
       stop("Configs and counts must have the same number of rows")
      if (ncol(configmatrix)<3) 
          stop("Less than three config variables make no sense")
#     Global variables needed to store results of recursive function perm.cfa
#     Global variables needed to store results of recursive function perm.cfa
      h.cfa.configs<<-list()
      h.cfa.chisqs<<-vector(mode="numeric")

      res<-perm.cfa(configmatrix,cntvector)
      lbls<-res[[1]]
      chisqs<-res[[2]]
      sortidx<-rev(order(chisqs))
      sorted.lbls<-lbls[sortidx]
      sorted.chisqs<-chisqs[sortidx]
      namevec<-""
      for (i in 1:length(sorted.lbls))
        {
          namevec<-c(namevec,paste(unlist(sorted.lbls[i]),collapse=" "))
          
        }
      names(sorted.chisqs)<-namevec[-1]
      orders<- unlist(lapply(strsplit(names(sorted.chisqs)," "),length))
#     Remove global variables       
      remove("h.cfa.configs","h.cfa.chisqs",inherits=TRUE)
      res<-list(chisq=sorted.chisqs,orders=orders)
      class(res)<-"hcfa"
      res
   }

print.hcfa<-function(x)
   {
      cat("\n*** Hierarchical CFA ***\n\n")
      restbl<-matrix(NA,nrow=length(x$chisq),ncol=2)
      restbl[,1]<-x$chisq
      restbl[,2]<-x$orders
      colnames(restbl)<-c("Overall chi squared","order")
      rownames(restbl)<-names(x$chisq)
      print(restbl)
      restbl
   }

# Bootstrap-CFA

# function name: boot.cfa

# arguments:
# configmatrix: matrix of factors or numbers
# cntvector   : vector of counts
# runs        : number of samples to draw
# bonferroni  : Bonferroni-adjust test for significance
# sig.limit   : Alpha limit for test of significance
#
# returns:
# Table with the following colums:
# Number of types
# Number of antitypes
# Percentage of types
# Number of cases found to be a significant type
# Percentage of cases found to be a significant type

boot.cfa<-function(configmatrix,cntvector,runs=100,bonferroni=F,sig.limit=0.05)
   {
     n<-sum(cntvector) 
     for (i in  1:runs)
       {
#        incvector<-cntvector
         incvector<-numeric(len=length(cntvector[,1]))
         for (j in 1:(nrow(cntvector)))
           {
             lim<-cntvector[j,1]/n
             rndnums<-runif(cntvector[[j,1]])
             lres<-(rndnums<=lim)
             inc<-length(lres[lres==TRUE])
#            incvector[[j,1]]<-incvector[[j,1]]+inc
             incvector[j]<-incvector[j]+inc
           }
         if (sum(incvector)==0) stop("All counts are zero")
         res<-cfa(configmatrix,as.data.frame(incvector) ,verbose=T)         
         limvec<-as.vector(res$table[,"p(z)"])
         if (bonferroni==TRUE)
          limit<-res$parms["sig.limit"]
         if (bonferroni==FALSE)
          limit<-sig.limit
         boolres<-(limvec<limit)
         if (i==1)
           {
             cnttype<-rep(0,length(res$table[,"n"]))
             cntleexpected<-rep(0,length(res$table[,"n"]))
             cntgtexpected<-rep(0,length(res$table[,"n"]))
           }
         idxleexpected<-res$table[,"n"]-res$table[,"expected"]<=0
         idxgtexpected<-res$table[,"n"]-res$table[,"expected"]>0
         cntleexpected[idxleexpected]<-cntleexpected[idxleexpected]+1
         cntgtexpected[idxgtexpected]<-cntgtexpected[idxgtexpected]+1
         cnttype[boolres]<-cnttype[boolres]+1
      }
    qcnttype<-cnttype/runs
    qtypes<-cntgtexpected/runs
    restbl<-cbind(cntleexpected,
                  cntgtexpected,
                  qtypes*100.0,
                  cnttype,
                  qcnttype*100.0)
    rownames(restbl)<-rownames(res$table)
    colnames(restbl)<-c("cnt.antitype",
                        "cnt.type",
                        "pct.types",
                        "cnt.sig",
                        "pct.cnt.sig")
    restbl
   } 



