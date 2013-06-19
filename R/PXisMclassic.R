
H<-function(m,a,b,n){
if(n>=b) invisible(choose(a,m)*choose(n-a,b-m)/choose(n,b))
else if(n>=a) invisible(choose(b,m)*choose(n-b,a-m)/choose(n,a))
else return(NaN)
}

# Old algorithm
PXisMclassic<-function(m,n,Nt,k){
k_1<-k-1
k_2<-k-2
r<-integer(k_1)
s<-integer(k_1)
n<-sort(n,decreasing=FALSE)
ret<-0
for(l in 1:k_2) r[l]<-sum(n[1:(l+1)])-l*Nt
r[k_2]<-max(r[k_2],m)
for(l in (k-3):1) r[l]<-max(r[l],r[l+1])
for(l in 1:k_2) s[l]<-min(n[1:(l+1)])
s[k_2]<-min(s[k_2],(Nt-n[k]+m))
for(l in (k-3):1) s[l]<-min(s[l],(Nt-n[l+2]+s[l+1]))
ks<-r
prod2<-1
jj<-k_2
repeat{
if(jj>2){
for(i in 3:jj) Hvalues[i]<-H(ks[i],ks[i-1],n[i+1],Nt)
prod2<-prod(Hvalues[3:k_2])
}
prodC<-prod2*H(ks[1],n[1],n[2],Nt)
if(k_2>1) prodC<-prodC*H(ks[2],ks[1],n[3],Nt)
ret<-ret+H(m,ks[k_2],n[k],Nt)*prodC
ks[1]<-ks[1]+1
jj<-1
while(jj<k_2){
ksa<-min(s[jj],(Nt-n[jj+2]+ks[jj+1]))
if(jj<k_2 && ks[jj]<=ksa){
jj<-jj+1
break
}
if(jj<k_2 && ks[jj]>ksa){
ks[jj+1]<-ks[jj+1]+1
for(j2 in jj:1) ks[j2]<-max(r[j2],ks[j2+1])
}
if(jj==k_2 && ks[jj]<=s[jj]) break
jj<-jj+1
}
if(ks[k_2]>s[k_2]) break
}
return(ret)
} #PXisMclassic








