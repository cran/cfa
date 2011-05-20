H<-function(m,a,b,n){
if(n>=b) invisible(choose(a,m)*choose(n-a,b-m)/choose(n,b))
else if(n>=a) invisible(choose(b,m)*choose(n-b,a-m)/choose(n,a))
else return(NaN)
}

# New algorithm
PXisM<-function(m,n,Nt,k){
k_1<-k-1
k_2<-k-2
r<-integer(k_1)
s<-integer(k_1)
n<-sort(n,decreasing=FALSE)
for(l in 1:k_2) r[l]<-sum(n[1:(l+1)])-l*Nt
r[k_2]<-max(r[k_2],m)
for(l in (k-3):1) r[l]<-max(r[l],r[l+1])
for(l in 1:k_2) s[l]<-min(n[1:(l+1)])
s[k_2]<-min(s[k_2],(Nt-n[k]+m))
for(l in (k-3):1) s[l]<-min(s[l],(Nt-n[l+2]+s[l+1]))
size<-max(s)+1
tmp<-rep(0,size)
for(i in r[k_2]:s[k_2]) tmp[i+1]<-H(m,i,n[k],Nt)
if(k_2>1){
for(l in k_2:2){
tmp2<-rep(0,size)
for(j in r[l]:s[l]){
for(i in max(r[l-1],j): min(s[l-1],(Nt-n[l+1]+j))) tmp2[i+1]<-tmp2[i+1]+tmp[j+1]*H(j,i,n[l+1],Nt)
}
tmp<-tmp2
}
}
ret<-0
for(i in r[1]:s[1]) ret<-ret+tmp[i+1]*H(i,n[2],n[1],Nt)
return(ret)
} #PXisM


