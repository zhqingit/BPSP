args = (commandArgs(trailingOnly = T))
print (args)
for (i in 1:length(args)){
	        eval(parse(text=args[[i]]))
}

##-branch point

polyn = read.table(file=dist_polyn,row.names=1,stringsAsFactors=F)
polyn = polyn+1
colnames(polyn) = c("back","bp","ppt")
totn = sum(polyn[,1])
polyn.t = polyn
polyn.t[,1] = polyn.t[,1]/totn
polyn.t[,3] = polyn.t[,3]/totn

p=apply(polyn.t,1,function(x,totn){
	x = as.numeric(x)
	p1=binom.test(x[2],totn,x[1],alternative="greater")$p.value
	p2=binom.test(x[2],totn,x[3],alternative="greater")$p.value
	#if (p1<=0.01 & p2 <=0.01) return(T)
	#else return(F)
	return(c(p1,p2))
},totn=totn)


pp = t(p)
#fdr[,1] = p.adjust(fdr[,1])
#fdr[,2] = p.adjust(fdr[,2])
colnames(pp) = c("p1","p2")
pp = as.data.frame(pp,stringsAsFactors=F)

pp.filter = subset(pp,p1<=0.05 & p2<=0.05)
polyn.filter = polyn[rownames(pp.filter),]
write.table(file=train_data,polyn.filter,row.names=T,col.names=T,sep="\t",quote=F)

