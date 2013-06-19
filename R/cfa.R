# cfa 0.9.3 

# Original part by Funke

# Note: If the argument is a dataframe, call as cfa(df[1:3],df[4])
# if it is a matrix, call as cfa(df[,1:3],df[,4])


scfa<-function(cfg,cnt=NA,sorton="chisq",sort.descending=TRUE,
              format.labels=TRUE)
   {
     if(is.na(cnt)[1]) cnt<-as.vector(rep(1,nrow(cfg)))
     if(!is.data.frame(cfg)) cfg<-as.data.frame(cfg)
     if (!is.data.frame(cnt)) cnt<-as.data.frame(cnt) 
     cnt<-unlist(cnt) 
     d<-aggregate(cnt,by=cfg,sum,na.rm=TRUE) 
     counts<-d[,ncol(d)] 
     n<-sum(cnt)
     configs<-ncol(cfg) 
     if (configs<2) stop("At least two configuration variables required!")
     if (format.labels==TRUE)
       {
         for (i in 1:configs) levels(d[,i])<-format(levels(d[,i]))
       }
     sums<-matrix(data=NA,nrow=nrow(d),ncol=configs)
     n.levels<-vector(mode="numeric",length=configs)
     string<-do.call("paste",d[,1:configs]) 
     strings<-strsplit(string,"\" *\"*")
     for (i in 1:configs)
       {
          if (!is.factor(d[,i])) d[,i]<-as.factor(d[,i]) 
          tbl<-aggregate(counts,d[i],sum)  
          n.levels[i]<-nrow(tbl)          
          dimnames(tbl)[1]<-tbl[1] 
#         handle "":          
          if (dimnames(tbl)[1][[1]][1]=="") dimnames(tbl)[1][[1]][1]<-"NA" 
          cidx<-as.character(d[,i])
          cidx[cidx==""]<-"NA" 
          sums[,i]<-tbl[cidx,2]
       }
     expected<-apply(sums,1,prod)/n^(configs-1)
     chisq<-(counts-expected)^2/expected 
     sortidx<-order(chisq)
     if (sorton=="chisq") sortidx<-order(chisq)
     if (sorton=="n") sortidx<-order(counts) 
     if (sorton=="label") sortidx<-order(unlist(strings))
     if (sort.descending==TRUE) sortidx<-rev(sortidx) 
     colnames(sums)<-names(cfg)
     rownames(sums)<-strings 
     reslist<-list(labels=unlist(strings[sortidx]),
                   n.levels=n.levels,
                   sums=sums[sortidx,],
                   counts=counts[sortidx],     
                   expected=expected[sortidx],
                   chisq=chisq[sortidx])
     reslist 
   }                                          


mcfa<-function(cfg,cnts,sorton="chisq",sort.descending=TRUE,
              format.labels=TRUE)
   {
     if(!is.data.frame(cfg)) cfg<-as.data.frame(cfg)
     configs<-ncol(cfg)
     if (configs<2) stop("At least two configuration variables required!")
     if (ncol(cnts)<2) stop("Only one column of counts - use cfa instead!")
     if (ncol(cnts)>2)  {options("warn"=-1)}
     if (format.labels==TRUE)
       {
         for (i in 1:configs) levels(cfg[,i])<-format(levels(cfg[,i]))
       }
     string<-do.call("paste",cfg[,1:configs]) 
     strings<-strsplit(string,"\" *\"*")
     n<-sum(cnts) 
     colsums<-apply(cnts,1,sum)
     rowsums<-apply(cnts,2,sum) 
     nsmat<-matrix(data=rep(rowsums,ncol(cnts)),nrow=nrow(cfg),ncol=ncol(cnts),byrow=TRUE)
     expected<-colsums*nsmat/n
     chisq<-apply((cnts-expected)^2/expected,1,sum)
     sortidx<-order(chisq)
     if (sorton=="chisq") sortidx<-order(chisq)
     if (sorton=="label") sortidx<-order(unlist(strings))
     if (sort.descending==TRUE) sortidx<-rev(sortidx)
     reslist<-list(labels=unlist(strings[sortidx]),
               counts=cnts[sortidx,],     
               expected=expected[sortidx,],
               chisq=chisq[sortidx]) 
     if (ncol(cnts)>2)  {options("warn"=0)}
     reslist 
   }

cfa<-function(cfg,cnts=NA,
              sorton="chisq",
              sort.descending=TRUE,
              format.labels=TRUE,
              casewise.delete.empty=TRUE,
              binom.test=FALSE,
              exact.binom.test=FALSE,
              exact.binom.limit=10,
              perli.correct=FALSE,
              lehmacher=FALSE,
              lehmacher.corr=TRUE,
              alpha=0.05,
              bonferroni=TRUE
              )
  {
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
         if ((ni[1]*ep[1])<=10.0) num<-num-0.5 # continuity correction
         num/den
       }  
     # The order of the following checks is important
     if (is.null(dim(cfg)))
       stop("cfg is probably part of a matrix! Use cfg[,1:n] rather than cfg[1:n]!") 
     if (is.na(cnts)[1]) cnts<-rep(1,nrow(cfg))
     if (length(cnts)==1 && !is.data.frame(cnts))
       stop("cnts is probably part of a matrix! Use cnts[,n] rather than cnts[n]!")

     if (casewise.delete.empty)  # delete cases with "" and NA in any variable
        {
           cfg<-as.data.frame(cfg)
           idx<-rep(TRUE,nrow(cfg)) 
           for (i in 1:ncol(cfg))
             { 
               interim<-idx & (factor(cfg[,i])!="") & (!is.na(cfg[,i]))
               idx<-interim 
             } 
           cfg<-cfg[idx,]
           if (is.null(dim(cnts))) 
             {
               cnts<-cnts[idx]
               ncnts<-1
             }
           else 
             {
               cnts<-cnts[idx,]
               ncnts<-ncol(cnts) 
             }
        } 
     ncases<-nrow(cfg)
     ntotal<-sum(cnts,na.rm=TRUE)
     if (!is.null(nrow(cnts)) && (nrow(cfg)!=nrow(cnts)))
        stop("Configurations and counts must have the same numer of cases!") 
     nconfigs<-ncol(cfg)
     n<-sum(cnts,na.rm=TRUE)
     multicfa<-FALSE
     if (!is.null(dim(cnts))) multicfa<-TRUE 

     if (multicfa==FALSE)  # ---------- One sample CFA
      {
        res<-scfa(cfg,cnts,sorton,sort.descending,format.labels)
        bivariate<-TRUE        
        if (any(res$n.levels>2)) bivariate<-FALSE
        configs<-length(res$chisq) 
        if (bivariate==TRUE)
          {
            df<-(2^nconfigs-nconfigs-1)    # (KL27&57)
            dfcell<-1                      # (KL 26)
          }
        if (bivariate==FALSE)
         {

           df<-prod(res$n.levels)-sum(res$n.levels)+(nconfigs-1)  # (LW96, Korrektur Gruener)
           dfcell<-prod(res$n.levels-1)                           # richtig?
         }    
       if (bivariate==TRUE)  
              z<-z.test(res$counts,res$expected/n,n)
       else
              z<-sqrt(res$chisq)
       pz<-1.0-pnorm(z)
       nz<-length(z)
       if (bonferroni==FALSE) nz<-1 # no adjustment 
       pzsig<-((1.0-pnorm(z,0,1))<(alpha/nz) | 
                 (pnorm(z,0,1)<(alpha/nz)))  # z: Type
       sig.chisq<-(1.0-pchisq(res$chisq,dfcell))<(alpha/nz)  # Chi^2: Type or antitype
       totalchisq<-sum(res$chisq,na.rm=TRUE)
       p.total.chisq<-1.0-pchisq(totalchisq,dfcell)
       if (binom.test==TRUE)   # LW 43 Approximate test for large n and "not too low" p
         { 
            z.bin<-(res$counts-res$expected)/sqrt(res$expected*(1-res$expected/n))
            low.limit<-qnorm(alpha/2)  # Binom: Type or antitype
            hi.limit<-qnorm(1.0-alpha/2)
            sig.bin<-((z.bin<low.limit) | (z.bin>hi.limit))
         } 
       if (exact.binom.test==TRUE && bivariate==TRUE)
         {
            p.exact.bin<-rep(NA,ncases)
            for (i in 1:ncases)
              {
                if (res$counts[i]<=exact.binom.limit) # Binom exact: Type
                  p.exact.bin[i]<-bin.test(res$counts[i],res$expected[i]/n,n)
              }
         }
       if (perli.correct==TRUE) # LW 42
          {
             numerator.corr<-rep(0,length(res$counts)) 
             idx<-((res$counts>5) & (res$counts<=10) & (res$counts-res$expected>0)) 
             numerator.corr[idx]<-(-0.5) 
             idx<-((res$counts>5) & (res$counts<=10) & (res$counts-res$expected<=0)) 
             numerator.corr[idx]<-0.5 
             x.perli<-(res$counts-res$expected+numerator.corr)/sqrt(res$expected) 
             low.limit<-qnorm(alpha/2)  # Perli z: Type or antitype
             hi.limit<-qnorm(1.0-alpha/2)
             sig.perli<-((x.perli<low.limit) | (x.perli>hi.limit))
          }
       if (lehmacher==TRUE)  # LW 44
         {
           p18<-apply((res$sums-1),1,prod)/((n-1)^nconfigs)
           pn<-apply(res$sums,1,prod)/(n^nconfigs)
           vl<-n*pn*(1-pn-(n-1)*(pn-p18))
           zl<-(res$counts-res$expected)/sqrt(vl) 
           # Kuechhoff's correction
           corr<-rep(0,length(res$counts)) 
           idx<-res$counts>res$expected
           corr[idx]<--0.5
           idx<-res$counts<res$expected
           corr[idx]<-0.5 
           zl.corr<-(res$counts-res$expected+corr)/sqrt(vl)
           low.limit<-qnorm(alpha/2)  # 2 sided test (type OR antitype)
           hi.limit<-qnorm(1.0-alpha/2)
           sig.zl<-((zl<low.limit) | (zl>hi.limit)) 
           sig.zl.corr<-((zl.corr<low.limit) | (zl.corr>hi.limit)) 
         } 
       q<-2*abs(res$counts-res$expected)/(n+abs(2*res$expected-n)) # Q: KL p. 34  
       # Table
       resdataframe<-data.frame(label=res$labels,
                                n=res$counts,
                                expected=res$expected,
                                Q=q, 
                                chisq=res$chisq,
                                p.chisq=1-pchisq(res$chisq,dfcell),
                                sig.chisq=sig.chisq,
                                z=z,
                                p.z=pz,
                                sig.z=pzsig)
      if (perli.correct==TRUE) resdataframe<-cbind(resdataframe,x.perli,sig.perli)                                
      if (lehmacher==TRUE && lehmacher.corr==FALSE) resdataframe<-cbind(resdataframe,zl,sig.zl) 
      if (lehmacher==TRUE && lehmacher.corr==TRUE) resdataframe<-cbind(resdataframe,zl.corr,sig.zl.corr)
      if (binom.test==TRUE) resdataframe<-cbind(resdataframe,z.bin,sig.bin)
      if (exact.binom.test==TRUE && bivariate==TRUE) resdataframe<-cbind(resdataframe,p.exact.bin)
      # Summary stats                                
       resvector<-c(totalchisq, df,p.total.chisq,n)
       names(resvector)<-c("totalchisq","df","p","sum of counts")
       # Levels of all configs
       reslevels<-res$n.levels 
       names(reslevels)<-colnames(cfg)
       res<-list(table=resdataframe,summary.stats=resvector,levels=reslevels)
       class(res)<-"scfa"
     }
      

     if (multicfa==TRUE)   # ----------- Multi sample CFA
      {
        res<-mcfa(cfg,cnts,sorton,sort.descending,format.labels) 
        bivariate<-TRUE        
        if (max(dim(table(cfg)))>2) bivariate<-FALSE 
        configs<-length(res$chisq) 
        counts<-ncol(cnts)
        factorn<-dim(table(cfg)) 
        if (bivariate==TRUE)
          {
            df<-(2^configs-1)*(counts-1) # (KL82)
            dfcell<-1                    # (KL77)
          }
        if (bivariate==FALSE)
          {
            df<-prod(factorn)+sum(-factorn)+(configs-1)*(counts-1) # richtig?
            dfcell<-prod(factorn-1)                                # richtig?
          }
          
       p.chisq<-1.0-pchisq(res$chisq,dfcell) 
       nz<-length(p.chisq)
       if (bonferroni==FALSE) nz<-1 # no adjustment
#      pzsig is Bonferroni-adjusted
       pchisqsig<-(p.chisq<(alpha/nz)) 
       totalchisq<-sum(res$chisq,na.rm=TRUE)
       p.total.chisq<-1.0-pchisq(totalchisq,dfcell) 
#      Table
       resdataframe<-data.frame(label=res$labels,
                                n=res$counts,
                                expected=res$expected,
                                chisq=res$chisq,
                                p.chisq=p.chisq,
                                p.chisq.sig=pchisqsig)
      # Summary stats                                
       resvector<-c(totalchisq, df,p.total.chisq,n)
       names(resvector)<-c("totalchisq","df","p","sum of counts") 
       res<-list(table=resdataframe,summary.stats=resvector)
       class(res)<-"mcfa"
      } 
     res
  }
 
hcfa<-function(configs,cnts)
   {   
      h.cfa.configs<-list()
      h.cfa.chisqs<-vector(mode="numeric") 
      h.cfa.dfs<-vector(mode="numeric")
      # delete subconfigurations and run CFA for the remaining combination
      perm.cfa<-function(configs,cnts)
       { 
         for (i in 1:ncol(configs))
          {
             configs.1<-configs[-i]
             res<-cfa(configs.1,cnts) 
             cfglbl<-paste(names(configs.1),sep=" ")
             h.cfa.configs<<-c(h.cfa.configs,list(cfglbl))
             h.cfa.chisqs<<-c(h.cfa.chisqs,res$summary.stats["totalchisq"])
             h.cfa.dfs<<-c(h.cfa.dfs,res$summary.stats["df"])
             if (ncol(configs.1)>=3)
               {
                 perm.cfa(configs.1,cnts)             
               }
          }
         list(h.cfa.configs,h.cfa.chisqs,h.cfa.dfs)
        } 
     configs<-as.data.frame(configs)
     cnts<-as.data.frame(cnts) 
     lcfg<-nrow(configs) 
     lcnt<-nrow(cnts)
     if (!is.data.frame(configs))
       stop("Configs must be a data frame")
     if (!is.data.frame(cnts)) 
       stop("Counts must be a data frame")   
     if (lcfg!=lcnt)
       stop("Configs and counts must have the same number of rows")
      if (ncol(configs)<3) 
          stop("Less than three config variables make no sense")

      res<-perm.cfa(configs,cnts) 
      lbls<-res[[1]]
      chisqs<-res[[2]]
      dfs<-res[[3]] 
      sortidx<-rev(order(chisqs))
      sorted.lbls<-lbls[sortidx]
      sorted.chisqs<-chisqs[sortidx]
      sorted.dfs<-dfs[sortidx]
      namevec<-""
      for (i in 1:length(sorted.lbls))
        {
          namevec<-c(namevec,paste(unlist(sorted.lbls[i]),collapse=" ")) 
        }
      names(sorted.chisqs)<-namevec[-1]
      orders<- unlist(lapply(strsplit(names(sorted.chisqs)," "),length)) 
      df<- unlist(lapply(strsplit(names(sorted.dfs)," "),length))
      res<-list(chisq=sorted.chisqs,df=sorted.dfs,orders=orders)
      class(res)<-"hcfa"
      res
   }
   
print.scfa<-function(x,...)
   { 
      cat("\n*** Analysis of configuration frequencies (CFA) ***\n\n")
      print(x$table)
      cat("\n\n")
      cat("Summary statistics:\n\n")
      cat("Total Chi squared         = ",x$summary.stats["totalchisq"],"\n") 
      cat("Total degrees of freedom  = ",x$summary.stats["df"],"\n") 
      cat("p                         = ",x$summary.stats["p"],"\n") 
      cat("Sum of counts             = ",x$summary.stats["sum of counts"],"\n")
      cat("\n")
      cat("Levels:\n\n")
      print(x$levels,...)
      cat("\n")
   }
   

print.mcfa<-function(x,...)
   { 
      cat("\n*** Repeated analysis of configuration frequencies (MCFA) ***\n\n")
      print(x$table,...)
      cat("\n\n")
      cat("Summary statistics:\n\n")
      cat("Total Chi squared         = ",x$summary.stats["totalchisq"],"\n") 
      cat("Total degrees of freedom  = ",x$summary.stats["df"],"\n") 
      cat("p                         = ",x$summary.stats["p"],"\n") 
      cat("Sum of counts             = ",x$summary.stats["sum of counts"],"\n")
      cat("\n")
   }

print.hcfa<-function(x,...)
   {
      cat("\n*** Hierarchical CFA ***\n\n")
      restbl<-matrix(NA,nrow=length(x$chisq),ncol=4)
      restbl[,1]<-x$chisq
      restbl[,2]<-x$df
      restbl[,3]<-1.0-pchisq(x$chisq,x$df) 
      restbl[,4]<-x$orders 
      colnames(restbl)<-c("Overall chi squared","df","p","order")
      rownames(restbl)<-names(x$chisq)
      print(restbl,...)
      restbl
   }

plot.scfa<-function(x,...)
   { 
     plot(unlist(x$table["n"]),unlist(x$table["chisq"]),
          main="Chi squared vs. n\n(left click to identify, right click to stop)",
          xlab="n",
          ylab="Chi squared",...
         )
      selected<-identify(unlist(x$table["n"]),unlist(x$table["chisq"]),
         n=nrow(x$table)) 
     if (length(selected)>0) print(x$table[c(selected),"label",drop=FALSE])
   }


plot.mcfa<-function(x,...)
   { 
     s<-colnames(x$table) 
     idx<-grep("n.counts",s)
     nsums<-apply(x$table[,idx],1,sum) 
     plot(nsums,unlist(x$table["chisq"]),
          main="Chi squared vs. sum of n\n(left click to identify, right click to stop)",
          xlab="sum of n",
          ylab="Chi squared",...
         )
      selected<-identify(nsums,unlist(x$table["chisq"]),
         n=nrow(x$table)) 
     if (length(selected)>0) print(x$table[c(selected),"label",drop=FALSE])
   }

plot.hcfa<-function(x,...)
   { 
     dotchart(unlist(x["chisq"]),groups=factor(unlist(x["orders"])),
              main="Chi squared by order",
              xlab="item",
              ylab="Chi squared",...
             )
   }

bcfa<-function(configs,cnts,runs=100,sig.item="sig.z",...)
   {
       if (is.null(dim(cnts))) 
             {
#               incvector<-numeric(len=length(cnts))
                incvector<-numeric(length(cnts)) 
             }
           else 
             { 
#               incvector<-numeric(len=length(cnts[,1]))
               incvector<-numeric(length(cnts[,1]))
             }    
      cnts<-unlist(cnts)       
      n<-sum(cnts) 
      for (i in  1:runs) 
        {
          for (j in 1:(length(cnts)))
           {
             lim<-cnts[j]/n
             rndnums<-runif(cnts[j])
             lres<-(rndnums<=lim)
             inc<-length(lres[lres==TRUE])
             incvector[j]<-incvector[j]+inc
           }
         if (sum(incvector)==0) stop("All counts are zero")
         # must be sorted on label to prevent re-ordering due to different n or chisq
         res<-cfa(configs,as.data.frame(incvector),sorton="label",...) 
         if (!sig.item %in% colnames(res$table)) stop("invalid sig.item specified!")
         boolres<-unlist(res$table[sig.item])
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
    rownames(restbl)<-as.character(res$table[,"label"])
    colnames(restbl)<-c("cnt.antitype",
                        "cnt.type",
                        "pct.types",
                        "cnt.sig",
                        "pct.cnt.sig")
    class(restbl)<-"bcfa"
    restbl
   } 
   
print.bcfa<-function(x,...)
   {
     print(x[,1:5],...) 
   }

plot.bcfa<-function(x,...)
   { 
      plot(unlist(x[,"cnt.type"]),unlist(x[,"cnt.sig"]),
           main="Count(sig.) by count(type)\n(left click to identify, right click to stop)",
           xlab="Count(type)",
           ylab="Count(sig.)",...
          )
      selected<-identify(unlist(x[,"cnt.type"]),unlist(x[,"cnt.sig"]),
                n=length(unlist(x[,"cnt.sig"]))) 
      names(selected)<-row.names(x)[selected]         
     if (length(selected)>0) print(selected)          
   }



