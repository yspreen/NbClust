source("NbClust.R")
#source("NbClust_orig.R")
file = "../test/s3-cb_2.job/file.txt"
path = "../test/"
mtb  =  t(read.table(file, header=TRUE, row.names=1, sep="\t")) 
maxnc  =  min(nrow(mtb) / 3, 20)
d="maximum"
m="centroid"
i="silhouette"

res  =  NbClust(data=mtb, diss=NULL, method=m, distance=d, index=i, min.nc=2, max.nc=maxnc)

zf = paste(path, paste(m,d,i, sep='_'), '.tsv', sep="")
zzf = file(zf, open = 'w+')
write.table(res$Best.partition, zzf)
close(zzf)
