###########################################################################################################################################
#SCVA of CCLE and BRCA of which each is added to CCLE distribution
source("/home/users/thomas0915/Research/pc_clustering/function.R")
parCores <- paste0("node", 13:22)
cl <- makeCluster(parCores)
mod2g <- MOD2G %>% filter(GENE_ID %in% as.numeric(rownames(BRCA.exp.tumor)))
MOD_ID.list <- unique(mod2g$MOD_ID) #52919
clusterExport(cl, c("mod2g", "MOD_ID.list"))
clusterEvalQ(cl, {library("multicore"); library("sfsmisc")})

num_seq <- list()
div <- floor(nrow(BRCA.exp.tumor)/12)
for(j in 1:11) num_seq[[j]] <- (div*(j-1)+1):(div*j)
num_seq[[12]] <- (div*11+1):nrow(BRCA.exp.tumor)

kde <<- apply(CCLE.dwd, 1, function(x) density(x, bw=sd(x)/4, kernel="gaussian"))
exp.matrix <<- CCLE.dwd
clusterExport(cl, c("exp.matrix", "kde"))
CCLE.norm <- do.call(rbind, clusterApplyLB(cl, num_seq, normalize.matrix.par))

num_seq2 <- list()
div <- floor(length(MOD_ID.list)/12)
for(k in 1:11) num_seq2[[k]] <- MOD_ID.list[(div*(k-1)+1):(div*k)]
num_seq2[[12]] <- MOD_ID.list[(div*11+1):length(MOD_ID.list)]

normalized.matrix <<- CCLE.norm; clusterExport(cl, "normalized.matrix")
CCLE.scva <- do.call(rbind, clusterApplyLB(cl, num_seq2, module.activation.chunk, method="am"))
rownames(CCLE.scva) <- MOD_ID.list

num_seq3 <- list()
div <- floor(ncol(BRCA.exp.tumor.dwd)/12)
for(l in 1:11) num_seq3[[l]] <-(div*(l-1)+1):(div*l)
num_seq3[[12]] <- (div*11+1):ncol(BRCA.exp.tumor.dwd)

add.matrix <<- BRCA.exp.tumor.dwd; standard.matrix <<- CCLE.dwd
clusterExport(cl, c("add.matrix", "standard.matrix"))
BRCA.exp.tumor.norm.add <- do.call(cbind, clusterApplyLB(cl, num_seq3, normalize.cell.matrix.par))

normalized.matrix <<- BRCA.exp.tumor.norm.add; clusterExport(cl, "normalized.matrix")
BRCA.scva.tumor.add <- do.call(rbind, clusterApplyLB(cl, num_seq2, module.activation.chunk, method="am"))
rownames(BRCA.scva.tumor.add) <- MOD_ID.list

stopCluster(cl)

cc.mod <- read.csv("/home/users/thomas0915/Research/pc_clustering/cancer.selected.mod.txt", sep= "\t")
cc.mod <- as.character(cc.mod$MOD_ID)

#stat of added normalized BRCA tumor expression
cor.add <- cor(BRCA.exp.tumor.norm.add, method="spearman")
cor.origin <- cor(BRCA.exp.tumor, method="spearman")
cor.origin.norm <- cor(BRCA.exp.tumor.norm, method="spearman")
summary(cor.add[lower.tri(cor.add)])
summary(cor.origin[lower.tri(cor.origin)])
summary(cor.origin.norm[lower.tri(cor.origin.norm)])
boxplot(cor.add[lower.tri(cor.add)], cor.origin[lower.tri(cor.origin)])
boxplot(cor.add[lower.tri(cor.add)], cor.origin.norm[lower.tri(cor.origin.norm)])
hist(cor.origin.norm[lower.tri(cor.origin.norm)])
hist(cor.add[lower.tri(cor.add)])

###########################################################################################################################################
##Parameter testing in samp-mod matrix
library("QUBIC")
input.cell.mod <- t(cbind(CCLE.scva, BRCA.scva.tumor.add)[cc.mod, ])
num.sel.bic <- c()
covered.cell <- c()
covered.mod <- c()
covered.bic.med <- c() # number of bicluster that each cell has
covered.bic.mean <- c()
for(r in 21:40){
  res.cell.mod <- biclust::biclust(input.cell.mod, method = BCQU(),
                                   r, q = 0.06,
                                   c = 0.95, o = 100000, f = 1,
                                   k = 10,
                                   type = 'default', P = FALSE, C = FALSE, verbose = FALSE,
                                   weight = NULL, seedbicluster = NULL)
  
  ##Differentially correlated biclusters(samples)
  col.QUBIC <- apply(res.cell.mod@NumberxCol, 1, function(x) (1:ncol(input.cell.mod))[x])
  row.QUBIC <- apply(res.cell.mod@RowxNumber, 2, function(x) (1:nrow(input.cell.mod))[x])
  module_dc_pval <- c()
  for(bicluster_num in 1:length(col.QUBIC)){
    case_mat <- input.cell.mod[row.QUBIC[[bicluster_num]], col.QUBIC[[bicluster_num]]]
    cont_mat <- input.cell.mod[row.QUBIC[[bicluster_num]], -col.QUBIC[[bicluster_num]]]
    
    case_cor <- cor(t(case_mat))
    case_cor_mean <- mean(case_cor[lower.tri(case_cor)])
    cont_cor <- cor(t(cont_mat))
    cont_cor_mean <- mean(cont_cor[lower.tri(cont_cor)])
    
    cor_diff 	<- case_cor_mean - cont_cor_mean 
    n  <- dim(case_mat)[2]
    n2 <- dim(cont_mat)[2]
    pval_fisher	<- paired.r(xy=case_cor_mean, xz=cont_cor_mean, n=n, n2=n2, twotailed=TRUE)
    module_dc_pval <- c(module_dc_pval, pval_fisher$p)
  }
  
  ##Differentially activated biclusters(modules)
  #make population distribution from all samples
  mu.samp_mean <- apply(input.cell.mod, 1, mean)
  sd.samp_mean <- apply(input.cell.mod, 1, sd)
  module_da_pval <- c()
  for(bicluster_num in 1:length(col.QUBIC)){
    idx.row <- row.QUBIC[[bicluster_num]]
    idx.col <- col.QUBIC[[bicluster_num]]
    mu_hat <- input.cell.mod[idx.row, idx.col] %>% apply(1, mean)
    
    z_hat <- (mu_hat-mu.samp_mean[idx.row])/sd.samp_mean[idx.row]
    pval_chisq <- pchisq(sum(z_hat^2), df=length(z_hat), lower.tail = F)
    module_da_pval <- c(module_da_pval, pval_chisq)
  }
  
  ##Both DCM & DAM
  idx.sel.bic <- which(module_da_pval < 0.05)[which(module_da_pval < 0.05) %in% which(module_dc_pval<0.05)]
  num.sel.bic <- c(num.sel.bic, length(idx.sel.bic))
  ##Covered rows(samples) number
  covered.cell <- c(covered.cell, length(unique(unlist(row.QUBIC[idx.sel.bic]))))
  covered.mod <- c(covered.mod, length(unique(unlist(col.QUBIC[idx.sel.bic]))))
  covered.bic.med <- c(covered.bic.med, median(as.numeric(table(unlist(row.QUBIC[idx.sel.bic])))))
  covered.bic.mean <- c(covered.bic.mean, mean(as.numeric(table(unlist(row.QUBIC[idx.sel.bic])))))
  message(r, " done")
}
which.max(num.sel.bic)
which.max(covered.cell)
which.max(covered.mod)
plot(num.sel.bic, type="o", ylab="Number of selected biclusters", xlab="r")
plot(covered.cell, type="o", ylab="Number of covered samples", xlab="r")
plot(covered.mod, type="o", ylab="Number of covered modules", xlab="r")
plot(covered.bic.med, type="o", ylab="Median of covered biclusters", xlab="r")
plot(covered.bic.mean, type="o", ylab="Mean of covered biclusters", xlab="r")

###########################################################################################################################################
## mod-samp matrix (make result2)
#extract Bicluster
input.mod.cellpat <- cbind(CCLE.scva, BRCA.scva.tumor.add)[cc.mod, ]
res <- biclust::biclust(input.mod.cellpat, method = BCQU(),
                        r=1, q = 0.06,
                        c = 0.95, o = 100000, f = 1,
                        k = 10,
                        type = 'default', P = FALSE, C = FALSE, verbose = TRUE,
                        weight = NULL, seedbicluster = NULL)
res
#Differentially correlated biclusters(samples)
col.QUBIC <- apply(res@NumberxCol, 1, function(x) (1:ncol(input.mod.cellpat))[x])
row.QUBIC <- apply(res@RowxNumber, 2, function(x) (1:nrow(input.mod.cellpat))[x])

cor.bicluster <- dc.mc.bicluster(row.QUBIC, col.QUBIC, input.mod.cellpat)
module_dc_pval <- cor.bicluster$module_dc_pval
module_mu_cor <- cor.bicluster$module_mu_cor
module_da_pval <- da.bicluster(row.QUBIC, col.QUBIC, input.mod.cellpat)

#Both HCM & DCM & DAM
tmp <- intersect(which(module_da_pval < 0.05), which(module_dc_pval < 0.05)) #1528
idx.sel.bic <- intersect(tmp, which(module_mu_cor > 0.9)) #736
idx.sel.pat.bic <- idx.sel.bic[unlist(lapply(col.QUBIC[idx.sel.bic], function(x) any(x <=526)))] #704
lapply(col.QUBIC, function(x) x[x > 1037]) %>% unlist %>% unique %>% length #451
lapply(col.QUBIC[idx.sel.pat.bic], function(x) x[x > 1037]) %>% unlist %>% unique %>% length #278

#combine overlapped modules in idx.sel.pat.bic
new.procd.comb <- comb.bicluster(row.QUBIC[idx.sel.pat.bic], col.QUBIC[idx.sel.pat.bic]) ##382##
comb.row.QUBIC <- new.procd.comb$comb.row.QUBIC
comb.col.QUBIC <- new.procd.comb$comb.col.QUBIC

# adapt filter one more time
comb.module_dc_pval <- dc.mc.bicluster(comb.row.QUBIC, comb.col.QUBIC, input.mod.cellpat)$module_dc_pval
comb.module_mu_cor <- dc.mc.bicluster(comb.row.QUBIC, comb.col.QUBIC, input.mod.cellpat)$module_mu_cor
comb.module_da_pval <- da.bicluster(comb.row.QUBIC, comb.col.QUBIC, input.mod.cellpat)
tmp <- intersect(which(comb.module_da_pval < 0.05), which(comb.module_dc_pval < 0.05)) #364
comb.idx.sel.bic <- intersect(tmp, which(comb.module_mu_cor > 0.9)) #354
comb.idx.sel.pat.bic <- comb.idx.sel.bic[unlist(lapply(comb.col.QUBIC[comb.idx.sel.bic], function(x) any(x <=526)))] ##354##
lapply(comb.col.QUBIC, function(x) x[x > 1037]) %>% unlist %>% unique %>% length #278
lapply(comb.col.QUBIC[comb.idx.sel.pat.bic], function(x) x[x > 1037]) %>% unlist %>% unique %>% length #259

hist(unlist(lapply(comb.row.QUBIC[comb.idx.sel.pat.bic], length)), breaks=seq(0,1000,10), freq=F,
     xlab="Number of modules in each bicluster", main="After combining overlapped bicluster")  
hist(unlist(lapply(row.QUBIC[idx.sel.pat.bic], length)), breaks=seq(0,1000,10), freq=F,
     xlab="Number of modules in each bicluster", main="Before combining overlapped bicluster")
hist(unlist(lapply(comb.col.QUBIC[comb.idx.sel.pat.bic], length)), breaks=seq(0,1000,10), freq=F,
     xlab="Number of samples in each bicluster", main="After combining overlapped bicluster")
hist(unlist(lapply(col.QUBIC[idx.sel.pat.bic], length)), breaks=seq(0,1000,10), freq=F,
     xlab="Number of samples in each bicluster", main="Before combining overlapped bicluster")

#bicluster profile of samples
CCLE.BRCA.bic <-
  t(sapply(comb.col.QUBIC[comb.idx.sel.pat.bic], function(x){
    bic.vector <- rep(0, ncol(input.mod.cellpat))
    bic.vector[x] <- 1
    return(bic.vector)
  }))
colnames(CCLE.BRCA.bic) <- colnames(input.mod.cellpat)
rownames(CCLE.BRCA.bic) <- paste0("BIC-", 1:length(comb.idx.sel.pat.bic))

CCLE.BRCA.bic <- CCLE.BRCA.bic[, -which(apply(CCLE.BRCA.bic, 2, sum)==0)]
sel.patient <- colnames(CCLE.BRCA.bic)[which(grepl("tcga", colnames(CCLE.BRCA.bic)))] #259
sel.cell <- colnames(CCLE.BRCA.bic)[-which(grepl("tcga", colnames(CCLE.BRCA.bic)))] #
sel.cell <- sel.cell[-which(duplicated(sel.cell))]
BRCA.sel.tumor.dwd <- BRCA.exp.tumor.dwd[, sel.patient]
CCLE.sel.dwd <- CCLE.dwd[, sel.cell]

apply(CCLE.BRCA.bic, 2, function(x) rownames(CCLE.BRCA.bic)[x==1])


##clustering based on selected bicluster profile of modular expression
#exclude samples not having bicluster
library(vegan)
dist.CCLE.BRCA.bic <- vegdist(t(CCLE.BRCA.bic), method = "jaccard")
nc <- NbClust(data=t(CCLE.BRCA.bic), diss=NULL, distance="binary", method="ward.D") ## optimal k = 5
barplot(table(nc$Best.n[1,]), xlab="Numer of Clusters", ylab="Number of Criteria", main="Number of Clusters Chosen")
plot(nc$All.index[,4], type="o", ylab="CCC", xlab="Number of clusters", col="blue")

dend.bic <- hclust(dist.CCLE.BRCA.bic, method="ward.D")
dend.bic <- as.dendrogram(dend.bic)

patient.sub <- unlist(sapply(colnames(CCLE.BRCA.bic), function(x) BRCA.clin[BRCA.clin$sample_id==x, ]$PAM50.mRNA))
patient.sub[is.na(patient.sub)] <- "Patient(no subtype)"

label <- c(cell.origin[colnames(CCLE.BRCA.bic)[1:675]], 
           "HER2-enriched","Luminal B", "Patient(no subtype)", "Normal-like","Basal-like", "Luminal A")
colors <- rainbow(length(unique(label)))
names(colors) <- unique(label)
edgecolors <- unlist(lapply(labels(dend.bic), function(x) if(grepl("tcga", x)){
  colors[patient.sub[x]]
}else{
  colors[cell.origin[colnames(CCLE.BRCA.bic)[1:675]][x]]
}))
names(edgecolors) <- labels(dend.bic)
samplecolors <- unlist(lapply(labels(dend.bic), function(x){if(grepl("tcga", x)) "black" else "grey"}), use.names = F)
names(samplecolors) <- labels(dend.bic)

dendrapply(dend.bic, function(n){
  if(is.leaf(n)){
    a <- attributes(n)
    attr(n, "edgePar") <- list(col=edgecolors[a$label])
  }
  n;
}) %>% set("labels", "") %>% plot()
# dend.bic %>% rect.dendrogram(k=5, border=8, lty=5, lwd=2, lower_rect = -0.05)
colored_bars(colors=samplecolors, rowLabels = "sample type")
legend(x=570,y=20, legend=c("Breast tumours", "Cell lines"), fill= c("black", "grey"), title="Sample type", cex=0.7)
legend("topright", legend=unique(names(colors)), col= colors[unique(names(colors))], lty=1, lwd=2, cex=0.53)

dendrapply(dend.bic[[1]], function(n){
  if(is.leaf(n)){
    a <- attributes(n)
    attr(n, "edgePar") <- list(col=edgecolors[a$label])
  }
  n;
}) %>% plot()
# dend.bic %>% rect.dendrogram(k=5, border=8, lty=5, lwd=2, lower_rect = -0.05)
legend("topright", legend=names(colors)[colors %in% edgecolors[labels(dend.bic[[1]])]], 
       col=colors[names(colors)[colors %in% edgecolors[labels(dend.bic[[1]])]]], 
       lty=1, lwd=2, cex=0.6)

cut(dend.bic, h = 15)$lower[[1]]


##t-SNE of g.exp
library(Rtsne)
label <- c(cell.origin[colnames(CCLE.sel.dwd)], rep("BREAST_PATIENT", ncol(BRCA.sel.tumor.dwd)))
colors.tsne <- rainbow(length(unique(label)))
names(colors.tsne) <- unique(label)

pdf("/home/users/thomas0915/Research/pc_clustering/Bicluster/result/tSNE.pdf", width=15)
res.tsne <- Rtsne(t(cbind(CCLE.sel.dwd, BRCA.sel.tumor.dwd)), dims=2, perplexity=30)
df.tsne.res <- data.frame(res.tsne$Y)
df.tsne.res <- cbind(df.tsne.res, label)
colnames(df.tsne.res) <- c("tSNE1", "tSNE2", "Label")
df.tsne.res %>% ggplot(aes(x=tSNE1, y=tSNE2, col=Label)) + geom_point(alpha=0.5) +
  scale_color_manual(values=colors.tsne) + ggtitle("t-SNE of CCLE & BRCA patients' gene expression") +
  theme(plot.title=element_text(hjust = 0.5))

##t-SNE of g.exp.norm
res.tsne <- Rtsne(t(cbind(CCLE.norm[, sel.cell], BRCA.exp.tumor.norm.add[, sel.patient])), dims=2, perplexity=30)
df.tsne.res <- data.frame(res.tsne$Y)
df.tsne.res <- cbind(df.tsne.res, label)
colnames(df.tsne.res) <- c("tSNE1", "tSNE2", "Label")
df.tsne.res %>% ggplot(aes(x=tSNE1, y=tSNE2, col=Label)) + geom_point(alpha=0.5) +
  scale_color_manual(values=colors.tsne) + ggtitle("t-SNE of CCLE & BRCA patients' relative gene expression") +
  theme(plot.title=element_text(hjust = 0.5))

##t-SNE of relative modular expression
res.tsne <- Rtsne(t(cbind(CCLE.scva[, sel.cell], BRCA.scva.tumor.add[, sel.patient])), dims=2, perplexity=30)
df.tsne.res <- data.frame(res.tsne$Y)
df.tsne.res <- cbind(df.tsne.res, label)
colnames(df.tsne.res) <- c("tSNE1", "tSNE2", "Label")
df.tsne.res %>% ggplot(aes(x=tSNE1, y=tSNE2, col=Label)) + geom_point(alpha=0.5) +
  scale_color_manual(values=colors.tsne) + ggtitle("t-SNE of CCLE & BRCA patients' relative modular expression") +
  theme(plot.title=element_text(hjust = 0.5))

##t-SNE of Bicluster profile
dist.CCLE.BRCA.bic.exdup <- vegdist(t(CCLE.BRCA.bic[, c(sel.cell, sel.patient)]), method = "jaccard")
res.tsne <- Rtsne(dist.CCLE.BRCA.bic.exdup, is_distance=TRUE, dims=2, perplexity=30)
df.tsne.res <- data.frame(res.tsne$Y)
df.tsne.res <- cbind(df.tsne.res, label)
colnames(df.tsne.res) <- c("tSNE1", "tSNE2", "Label")
df.tsne.res %>% ggplot(aes(x=tSNE1, y=tSNE2, col=Label)) + geom_point(alpha=0.5) +
  scale_color_manual(values=colors.tsne) + ggtitle("t-SNE of CCLE & BRCA patients' Bicluster profile") +
  theme(plot.title=element_text(hjust = 0.5))
dev.off()

##Survival test
library("survival")
BRCA.clin$ev <- 1*(BRCA.clin$vital_status == 1)
BRCA.clin$fut <- as.numeric(BRCA.clin$days_to_last_followup)
index.na <- which(is.na(BRCA.clin$fut))
BRCA.clin[index.na, "fut"] <- as.numeric(BRCA.clin[index.na, "days_to_death"])
BRCA.clin <- BRCA.clin[-which(is.na(BRCA.clin$fut), BRCA.clin$fut),]

#survival among bic cluster index
cluster.idx <- cutree(dend.bic, 5)[grepl("tcga", names(cutree(dend.bic, 5)))]
cluster.idx <- cluster.idx[cluster.idx!=4] #exclude cluster4 because it has only one patient
BRCA.clin.bic <- BRCA.clin
rownames(BRCA.clin.bic) <- BRCA.clin.bic$sample_id
BRCA.clin.bic <- BRCA.clin.bic[names(cluster.idx), ]

su <- Surv(BRCA.clin.bic$fut, BRCA.clin.bic$ev)
surv <- survfit(su~factor(cluster.idx))
plot(surv, col=c("red", "blue"))
survdiff(su~cluster.idx)
res.coxph <- coxph(su~cluster.idx)
exp(coef(res.coxph)) #HR=0.902, p >0.05

#HR of each relative modular expression in BRCA.bic samples
coef <- c()
pval <- c()
for(i in 1:nrow(input.mod.cellpat)){
  mod.exp <- input.mod.cellpat[rownames(input.mod.cellpat)[i], ]
  mod.exp <- mod.exp[names(cluster.idx)]
  mod.exp <- mod.exp*100 # for 0.01 unit change
  res.coxph <- coxph(as.formula("su ~ mod.exp"))
  coef <- c(coef, coef(res.coxph))
  pval <- c(pval, summary(res.coxph)$coefficients[5])
}
length(which(pval < 0.05)) #1221
summary(coef[which(pval < 0.05)])
summary(exp(coef[which(pval < 0.05)]))
df.HR.mod.all <- data.frame(module.idx=which(pval < 0.05),
                            coef=coef[which(pval < 0.05)],
                            HR=exp(coef[which(pval < 0.05)]))
tail(df.HR.mod.all)
#HR of each relative modular expression in each bic cluster
coef <- c()
pval <- c()
for(i in 1:nrow(input.mod.cellpat)){
  mod.exp <- input.mod.cellpat[rownames(input.mod.cellpat)[i], ]
  mod.exp <- mod.exp[names(cluster.idx)]
  mod.exp <- mod.exp*100
  res.coxph <- coxph(as.formula("su ~ mod.exp + strata(cluster.idx)"))
  coef <- c(coef, coef(res.coxph))
  pval <- c(pval, summary(res.coxph)$coefficients[5])
}
length(which(pval < 0.05)) #922
summary(coef[which(pval < 0.05)])
summary(exp(coef[which(pval < 0.05)]))
df.HR.mod.all.strata <- data.frame(module.idx=which(pval < 0.05),
                            coef=coef[which(pval < 0.05)],
                            HR=exp(coef[which(pval < 0.05)]))
tail(df.HR.mod.all.strata)

summary(df.HR.mod.all$HR)
hist(df.HR.mod.all$HR, breaks=seq(0.5,1.5,0.01))
summary(df.HR.mod.all.strata$HR)
hist(df.HR.mod.all.strata$HR, breaks=seq(0.5,1.5,0.01))

#HR of modules in biclusters which are enriched in each cluster
CCLE.BRCA.bic.pat <- CCLE.BRCA.bic[, names(cluster.idx)]
#cluster 1 (modules of max overlapped bicluster)
# mod.idx.c1 <- which.max(apply(CCLE.BRCA.bic.pat[, names(cluster.idx[cluster.idx==1])], 1, mean))
bic.idx.c1 <- which(apply(CCLE.BRCA.bic.pat[, names(cluster.idx[cluster.idx==1])], 1, mean) > 0.05)
mean(unique(unlist(comb.row.QUBIC[comb.idx.sel.pat.bic][bic.idx.c1])) %in% df.HR.mod.all.strata$module.idx)

#cluster 5 (modules of max overlapped bicluster)
# mod.idx.c5 <- which.max(apply(CCLE.BRCA.bic.pat[, names(cluster.idx[cluster.idx==5])], 1, mean))
bic.idx.c5 <- which(apply(CCLE.BRCA.bic.pat[, names(cluster.idx[cluster.idx==5])], 1, mean) > 0.05)
mean(unique(unlist(comb.row.QUBIC[comb.idx.sel.pat.bic][bic.idx.c5])) %in% df.HR.mod.all.strata$module.idx)

#HR of lasso combination of modules in biclusters which are enriched in each cluster
library(penalized)
set.seed(34)
hepato.prof <- profL1(Surv(OS, Death),
                      penalized=hepatoCellularNoMissing[,23:48],
                      standardize=T, fold=10, minlambda1=2, maxlambda1=12)
plot(hepato.prof$cvl ~ hepato.prof$lambda, type="l", log="x",
     xlab="lambda", ylab="Cross-validated log partial likelihood")
hepato.opt <- optL1(Surv(OS, Death),
                    penalized=hepatoCellularNoMissing[,23:48], standardize=T,
                    fold=10)
hepato.opt$lambda
abline(v=hepato.opt$lambda, col="gray")
epato.pen <- penalized(Surv(OS, Death),
                       penalized=hepatoCellularNoMissing[,23:48], standardize=T,
                       lambda1=hepato.opt$lambda)

unique(unlist(comb.row.QUBIC[comb.idx.sel.pat.bic][bic.idx.c1]))



#PAM50 information of bic cluster
table.c1 <- table(BRCA.clin.bic[names(cluster.idx[cluster.idx==1]), ]$PAM50.mRNA)
table.c5 <- table(BRCA.clin.bic[names(cluster.idx[cluster.idx==5]), ]$PAM50.mRNA)
df.pam50 <- data.frame(cluster=c(rep(1, sum(table.c1)), rep(5, sum(table.c5))),
                       PAM50=c(rep("Basal-like", 41), rep("HER2-enriched", 13), rep("Luminal A", 97),
                               rep("Luminal B", 42), rep("Normal-like", 3),
                               rep("Basal-like", 25), rep("HER2-enriched", 7), rep("Luminal A", 10),
                               rep("Luminal B", 15), rep("Normal-like", 2)))
df.pam50 %>% ggplot(aes(x=factor(cluster), fill=PAM50))+geom_bar(position="fill")


# #perform hc & select optimal k that satisfy largest silhouette index in tSNE space
# sil.tSNE <- c()
# for(i in 2:30){
#   dend.tSNE <- hclust(dist(res.tsne$Y), method="ward.D")
#   class.tSNE <- cutree(dend.tSNE, i)
#   sil.tSNE <- c(sil.tSNE, mean(silhouette(class.tSNE, dist=dist(res.tsne$Y))[,3]))
# }
# plot(sil.tSNE, type = "o")
# k <- which.max(sil.tSNE) 
# k <- k+1
# #Calculate silhouette index of tSNE result in real g.exp space
# dend.tSNE <- hclust(dist(res.tsne$Y), method="ward.D")
# class.tSNE <- cutree(dend.tSNE, k)
# sil.tSNE.gexp <- mean(silhouette(class.tSNE, dist=dist(t(cbind(CCLE.sel.dwd, BRCA.sel.tumor.dwd))))[, 3])
# 
# ##hc of g.exp
# #perform hc & select optimal k that satisfy largest silhouette index in g.exp space
# library(foreach)
# library(doParallel)
# registerDoParallel(cores=29)
# sil.gexp <- 
#   foreach(i=2:30, .combine = c) %dopar%{
#     dend.gexp <- hclust(dist(t(cbind(CCLE.sel.dwd, BRCA.sel.tumor.dwd))), method="ward.D")
#     class.gexp <- cutree(dend.gexp, i)
#     mean(silhouette(class.gexp, dist=dist(t(cbind(CCLE.sel.dwd, BRCA.sel.tumor.dwd))))[,3])
#   }
# plot(sil.gexp, type = "o")
# k2 <- which.max(sil.gexp) + 1 
# sil.gexp <- sil.gexp[which.max(sil.gexp)]
# 
# ##hc of mod.act
# #perform hc & select optimal k that satisfy largest silhouette index in module activation space
# sil.mod <- 
#   foreach(i=2:30, .combine = c) %dopar%{
#     dend.mod <- hclust(dist(t(cbind(CCLE.scva, BRCA.scva.tumor.add)[cc.mod, ])), method="ward.D")
#     class.mod <- cutree(dend.mod, i)
#     mean(silhouette(class.mod, dist=dist(t(cbind(CCLE.scva, BRCA.scva.tumor.add)[cc.mod, ])))[,3])
#   }
# stopImplicitCluster()
# plot(sil.mod, type = "o")
# k3 <- which.max(sil.mod)
# k3 <- k3 + 1
# #Calculate silhouette index of tSNE result in real g.exp space
# dend.mod <- hclust(dist(t(cbind(CCLE.scva, BRCA.scva.tumor.add)[cc.mod, ])), method="ward.D")
# class.mod <- cutree(dend.mod, k3)
# sil.mod.gexp <- mean(silhouette(class.mod, dist=dist(t(cbind(CCLE.sel.dwd, BRCA.sel.tumor.dwd))))[, 3])
# 
# sil.gexp
# sil.tSNE.gexp
# sil.mod.gexp


# #stat of selected modules
# length(unique(unlist(row.QUBIC[idx.sel.pat.bic]))) #6301 modules
# length(unique(unlist(col.QUBIC[idx.sel.pat.bic]))) #974 samples
# 
# idx.sel.pat.bic[which(unlist(lapply(col.QUBIC[idx.sel.pat.bic], function(x) sum(x <=526)/length(x) <= 0.05)))]
# col.QUBIC[1183]
# MOD[MOD$ID %in% as.numeric(cc.mod[unlist(row.QUBIC[1183])]), ] %>% as.data.frame

# #test new procedure(QUBIC->combine overlap->DCM, HCM, DAM)
# new.procd.comb <- comb.bicluster(row.QUBIC, col.QUBIC)
# new.comb.row.QUBIC <- new.procd.comb$comb.row.QUBIC
# new.comb.col.QUBIC <- new.procd.comb$comb.col.QUBIC
# 
# new.module.dc.mc <- dc.mc.bicluster(new.comb.row.QUBIC, new.comb.col.QUBIC, input.mod.cellpat)
# new.module_dc_pval <- new.module.dc.mc$module_dc_pval
# new.module_mu_cor <- new.module.dc.mc$module_mu_cor
# 
# new.module_da_pval <- da.bicluster(new.comb.row.QUBIC, new.comb.col.QUBIC, input.mod.cellpat)
# 
# tmp <- intersect(which(new.module_da_pval < 0.05), which(new.module_dc_pval < 0.05)) #379
# new.idx.sel.bic <- intersect(tmp, which(new.module_mu_cor > 0.9)) #238
# new.idx.sel.pat.bic <- new.idx.sel.bic[unlist(lapply(new.comb.col.QUBIC[new.idx.sel.bic], function(x) any(x <=526)))] ##211##
# lapply(new.comb.col.QUBIC, function(x) x[x > 1037]) %>% unlist %>% unique %>% length #451
# lapply(new.comb.col.QUBIC[new.idx.sel.pat.bic], function(x) x[x > 1037]) %>% unlist %>% unique %>% length #234
# 
# hist(unlist(lapply(new.comb.row.QUBIC, length)), breaks=seq(0,2000,10), freq=F,
#      xlab="Number of modules in each bicluster", main="After combining overlapped bicluster in new procedure")  
# hist(unlist(lapply(new.comb.col.QUBIC, length)), breaks=seq(0,1000,10), freq=F,
#      xlab="Number of samples in each bicluster", main="After combining overlapped bicluster in new procedure")  
# ##decide to not use this procedure because of less covered patient and very large bicluster



