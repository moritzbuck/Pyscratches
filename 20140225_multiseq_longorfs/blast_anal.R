pat="/home/moritz/GEFES/views/projects/acI/binning/clustering/kmeans/tetramer_clusters/nb_clusts-30-min_len-5000-max_freq-0.15/"
files.acI=grep("blast",list.files(pat),val=T)
files.acI=paste(pat,files.acI,sep="")
files.acI=files.acI[sapply(files.acI,function(x) file.info(x)$size)!=0]

data.acI=lapply(files.acI,read.table)
nb_contigs.acI=sapply(sub(".blast","",files.acI),function(x) as.numeric(system(paste(c("grep \\> ",x," | wc -l"),collapse=""),intern=T)))
length.acI=sapply(sub(".blast","",files.acI),function(x) as.numeric(system(paste(c("grep -v \\> ",x," | wc -c"),collapse=""),intern=T)))
nb_genes.acI=sapply(sub(".blast",".genes",files.acI),function(x) as.numeric(system(paste(c("grep \\> ",x," | wc -l"),collapse=""),intern=T)))

min.evs.acI=lapply(data.acI,function(dat) sapply(levels(dat$V1),function(x) min(dat[dat$V1==x,]$V11)))

pat="/home/moritz/GEFES/views/projects/alinen/binning/clustering/kmeans/tetramer_clusters/nb_clusts-16-min_len-1000-max_freq-0.12-transform-rank/"
#pat="/home/moritz/GEFES/views/projects/alinen/binning/clustering/kmeans/tetramer_clusters/nb_clusts-8-min_len-500-max_freq-0.12-transform-rank/"
files.alinen=grep("blast",list.files(pat),val=T)
files.alinen=paste(pat,files.alinen,sep="")
files.alinen=files.alinen[sapply(files.alinen,function(x) file.info(x)$size)!=0]

#setwd(pat)
data.alinen=lapply(files.alinen,read.table)
nam=sapply(strsplit(files.alinen,"-rank/"),function(x) x[2])
#max.bit.alinen=lapply(data.alinen,function(dat) sapply(levels(dat$V1),function(x) max(dat[dat$V1==x,]$V12)))


#logs.evs.alinen=sapply(min.evs.alinen,function(x) -log10(x)[-log10(x) > 1 & !is.infinite(-log10(x)) ])
length.alinen=sapply(sub(".blast","",files.alinen),function(x) as.numeric(system(paste(c("grep -v \\> ",x," | wc -c"),collapse=""),intern=T)))
nb_contigs.alinen=sapply(sub(".blast","",files.alinen),function(x) as.numeric(system(paste(c("grep \\> ",x," | wc -l"),collapse=""),intern=T)))
nb_genes.alinen=sapply(sub(".blast",".genes",files.alinen),function(x) as.numeric(system(paste(c("grep \\> ",x," | wc -l"),collapse=""),intern=T)))

contigs=sapply(min.evs.alinen,function(x) factor(sapply(strsplit(names(x),";"),function(x) x[2])))
mean_intracontig_iqr = sapply(1:8,function(x) mean(sapply(levels(contigs[[x]]),function(c) IQR(inf.rm(-log10(min.evs.alinen[[x]][contigs[[x]]==c]))))))
binning_iqr = sapply(1:8,function(x) IQR(-log10(min.evs.alinen[[x]])))


mean_intracontig_sd = sapply(1:8,function(x) mean(sapply(levels(contigs[[x]]),function(c) sd(inf.rm(-log10(min.evs.alinen[[x]][contigs[[x]]==c]))))))
binning_sd = sapply(1:8,function(x) sd(-log10(min.evs.alinen[[x]])))


for j in `ls *fasta.genes`
do
for i in `ls ~/glob/data/blast.dbs/40scg_from_specI/*.fna | cut -f8 -d \/`
do
mkdir $j.blast
blastall -p tblastx -d /home/moritz/glob/data/blast.dbs/40scg_from_specI/$i -i $j -o $j.blast/$i -Y 9.61e+14 -a 16 -m8
done
done 

for i in `ls ~/glob/data/blast.dbs/40scg_from_specI/*.fna | cut -f8 -d \/`
do
blastall -p tblastx -d /home/moritz/glob/data/blast.dbs/40scg_from_specI/$i -i contigs.genes.fasta -o contigs.genes.fasta.blast/$i -Y 2.13e+14 -a 16 -m8
done
