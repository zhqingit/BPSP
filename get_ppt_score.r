
args = (commandArgs(trailingOnly = T))
print (args)
for (i in 1:length(args)){
       	eval(parse(text=args[[i]]))
}

##-ppt
polyn = read.table(file=dist_polyn,row.names=1,stringsAsFactors=F)
colnames(polyn) = c("back","bp","ppt")
polyn = polyn+1
sc = polyn[,"ppt"]^2/polyn[,"back"]
sc = sc/sum(sc)
sc = log(sc)
polyn.filter = cbind(polyn,sc)
write.table(file=PPTscore,polyn.filter,row.names=T,col.names=T,sep="\t",quote=F)
