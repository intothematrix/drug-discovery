setwd("~/Documents/Gene")
source('Functions.R')
load.libraries()
df<-load.data()

cor.result.p <- df %>% tidyr::gather(R1,R1v,contains("dependency")) %>%
  tidyr::gather(R2,R2v,contains("log_copy_number"),-cell_line) %>%
  dplyr::group_by(gene_name,R1,R2) %>%
  dplyr::summarize(pearson = cor(x=R1v, y=R2v,method = 'pearson',use = "pairwise.complete.obs")) %>%
  tidyr::unite(Pair, R1, R2, sep="_")

cor.result.s <- df %>% tidyr::gather(R1,R1v,contains("dependency")) %>%
  tidyr::gather(R2,R2v,contains("log_copy_number"),-cell_line) %>%
  dplyr::group_by(gene_name,R1,R2) %>%
  dplyr::summarize(spearman = cor(x=R1v, y=R2v,method = 'spearman',use = "pairwise.complete.obs")) %>%
  tidyr::unite(Pair, R1, R2, sep="_")

cor.result.p %>% dplyr::inner_join(cor.result.s,by=c("gene_name","Pair"))->cor.result
rm(cor.result.p)
rm(cor.result.s)

cor.result<-cor.result %>% dplyr::select(gene_name, pearson, spearman)

# summary
summary(cor.result$pearson)
summary(cor.result$spearman)
hist(cor.result$pearson, br=20)
hist(cor.result$spearman, br=20)
cor.result[cor.result$pearson <= -0.2,]$gene_name ->gene.0.2
cor.result[cor.result$pearson <= -0.3,]$gene_name ->gene.0.3
cor.result[cor.result$pearson <= -0.4,]$gene_name ->gene.0.4
df %>% dplyr::group_by(gene_name) %>%
  dplyr::summarise(n_cl=n_distinct(cell_line),
                   mx.cn=max(log_copy_number),
                   mn.cn=min(log_copy_number),
                   mdn.cn = median(log_copy_number),
                   rng.cn=max(log_copy_number)-min(log_copy_number),
                   avg.cn=mean(log_copy_number),
                   mm.cn=median(log_copy_number)-mean(log_copy_number),
                   mx.dp=max(dependency),
                   mn.dp=min(dependency),
                   mdn.dp = median(dependency),
                   rng.dp=max(dependency)-min(dependency),
                   avg.dp=mean(dependency),
                   mm.dp=median(dependency)-mean(dependency)
  )->gene.features

gene.features %>%dplyr::inner_join(cor.result, by=("gene_name"="gene_name")) ->gene.features
set.seed(123)
clus<-kmeans(scale(gene.features[,2:16]), centers=20)
plotcluster(gene.features[,2:16], clus$cluster)
clus$cluster[gene.features$gene_name == "ARC"]
clus$cluster[is.element(gene.features$gene_name,gene.0.4)]
clus$cluster[is.element(gene.features$gene_name,gene.0.3)]
clus$cluster[is.element(gene.features$gene_name,gene.0.2)]
table(clus$cluster)
view(clus$centers)
write.csv(clus$centers, file="centers.csv")
write.csv(clus$cluster, file="clusters.csv")

# Determine number of clusters
wss <- c(1:50)
for (i in 2:50) wss[i] <- kmeans(scale(gene.features[,2:16]), centers=i)$tot.withinss
plot((2:50), wss[2:50], type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#need to do this with a sampling method so it's not that big
ggpairs(gene.features[,2:16]) 

#figure out how to develop features to categorize copynumber distribution
num.peak<-function(x,br=100,p=0.01){
  hist(x,br=br,plot=FALSE)->y
  l<-length(y$mids[intersect(which(diff(sign(diff(y$counts)))==-2),which(y$counts/sum(y$counts)>=p))])
  return(l)
}

df %>% dplyr::group_by(gene_name) %>%
  dplyr::summarise(npeak=num.peak(log_copy_number,br=50,p=0.05)
  )->junk
table(junk$npeak)
junk[junk$gene_name=="ERBB2",]
for (i in 1:16){
  a=junk[junk$npeak == 4,]$gene_name[i]
  hist(df[df$gene_name==a,]$log_copy_number,br=20)
}


hist(df[df$gene_name=="ERBB2",]$log_copy_number,br=100)
set.seed(234)
kmeans(df[df$gene_name=="ERBB2",]$log_copy_number,centers=5)->junk
table(junk$cluster)
table(junk$centers)
abline(v=junk$centers)
length(df[df$gene_name == "ERBB2"&df$log_copy_number<=sort(junk$centers)[1],]$log_copy_number)/length(df[df$gene_name == "ERBB2",]$log_copy_number)
length(df[df$gene_name == "ERBB2"&df$log_copy_number<=sort(junk$centers)[2]&df$log_copy_number>sort(junk$centers)[1],]$log_copy_number)
length(df[df$gene_name == "ERBB2"&df$log_copy_number<=sort(junk$centers)[3]&df$log_copy_number>sort(junk$centers)[2],]$log_copy_number)
length(df[df$gene_name == "ERBB2"&df$log_copy_number<=sort(junk$centers)[4]&df$log_copy_number>sort(junk$centers)[3],]$log_copy_number)
length(df[df$gene_name == "ERBB2"&df$log_copy_number<=sort(junk$centers)[5]&df$log_copy_number>sort(junk$centers)[4],]$log_copy_number)
length(df[df$gene_name == "ERBB2"&df$log_copy_number>sort(junk$centers)[5],]$log_copy_number)

hist(df[df$gene_name=="ERBB2",]$log_copy_number,br=100,plot=FALSE)->y
y$density[y$density>=0.01]
peakx <- y$mids[which(diff(sign(diff(y$density)))==-2)]



# num cl,
# max range within cl,
# range within cl,
# median range within cl
# mean range within cl,
# mean-median within cl
# overall range

#crispr <- eh[["EH3081"]]
#drug_sensitivity <- eh[["EH3087"]]
#mutationCalls <- eh[["EH3085"]]
#metadata <- eh[["EH3086"]]
#TPM <- eh[["EH3084"]]
#tidyr::pivot_wider(names_from = cell_line, values_from = dependency)
#tidyr::pivot_wider(names_from = cell_line, values_from = log_copy_number)
