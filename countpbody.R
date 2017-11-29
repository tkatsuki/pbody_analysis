library(dipr)
library(EBImage)
library(ggplot2)
dir <- "G:/Yuzuru/Colocalization/171003 P0 DDx6 Dcp1a CDH1/WT 100x Zm1 -2.oif.files"

pbodynum <- function(dir, mol, molthresh=0.1, pbodythresh=0.05, cadch, pbodych, molch){
  filelist <- list.files(dir, pattern=".tif")

  npbody <- list()
  ngerm <- list()
  npbody_per_germ <- list()
  pbodysize <- list()
  intpbody <- list()
  intgerm <- list()
  cells <- list()
  #cadch <- 3
  #pbodych <- 2
  #molch <- 1
  #molthresh=0.1
  #pbodythresh=0.05
  
    cad <- readImage(paste0(dir, "/", filelist[4]))
    pbody <-  readImage(paste0(dir, "/", filelist[3]))
    h <- dim(cad)[1]
    w <- dim(cad)[2]
    
    # Creating germ mask
    cadblur <- gblur(cad, 5) 
    cmask2 <- thresh(cadblur, 15, 15, 0.001)
    cmask3 <- bwlabel(cmask2)
    csize1 <- computeFeatures.shape(cmask3)
    cid1 <- which(csize1[,'s.area'] < 4000)
    cmask4 <- rmObjects(cmask3, cid1)
    cmask4 <- cmask4 > 0
    cmask5 <- thinning(cmask4)
    cmask6 <- fillHull(cmask5)
    cmask7 <- cmask6 - cmask5
    cmask8 <- bwlabel(cmask7)
    csize2 <- computeFeatures.shape(cmask8)
    cid2 <- which(csize2[,'s.area'] < 4000)
    cmask9 <- rmObjects(cmask8, cid2)
    cmask9er <- erode(cmask9, makeBrush(3, shape="diamond"))
    germftr <- computeFeatures.shape(cmask9er)
    rid <- which(germftr[, 's.radius.sd'] > 11)
    cmask10 <- rmObjects(cmask9, rid)
    ngerm <- max(cmask10)
    cmask11 <- cmask10 > 0
    cmask12 <- Image(array(c(cmask10, 0*cmask10, cmask10), dim=c(h, w, 3)))
    colorMode(cmask12) <- "Color"
    cadc <- Image(array(c(cad, cad, cad), dim=c(h, w, 3)))
    colorMode(cadc) <- "Color"
    
    # Creating p-body mask
    pbodyint <- computeFeatures.basic(cmask11, pbody)[1]
    pbodybl <- gblur(pbody, 3, 5)
    pbodybg <- gblur(pbody, 5, 9)
    pbodydn <- pbodybl - pbodybg
    pbodythresh <- 0.5*pbodyint
    pmask <- thresh(pbodydn, 100, 100, pbodythresh)
    pmask3 <- bwlabel(pmask)
    psize <- computeFeatures.shape(pmask3)
    pid <- which(psize[,'s.area'] < 4)
    pmask4 <- rmObjects(pmask3, pid)
    pmask5 <- pmask4 > 0
    pmask6 <- Image(array(c(pbody+pmask5, pbody-pmask5, pbody+pmask5), dim=c(h, w, 3)))
    colorMode(pmask6) <- "Color"
    
    # 
    cellst <- stackObjects(cmask10, cmask10) > 0
    edgest <- dilate(cellst, makeBrush(3, "diamond")) - cellst
    pbodyst <- stackObjects(cmask10, pmask5)
    cellstint <- apply(stackObjects(cmask10, img), 3, mean)
    mergest <- rgbImage(edgest, edgest+(pbodyst>0),edgest)
    pbodystnum <-  apply(bwlabel(pbodyst), 3, max)
    cells <- data.frame(file=rep(paste(dir, "/result", sep="")), 
                             pbodystnum, cellstint)
    
    # Mask pbody with germ cell mask
    germpbody <- pmask5 * cmask11

    # Save processed images
    writeImage(rgbImage(pbody, pbody, pbody), paste(dir,  "/pbody.png", sep=""))
    writeImage(pbodybg, paste(dir, "/pbodybg.png", sep=""))
    writeImage(pbodydn, paste(dir, "/pbodydn.png", sep=""))
    writeImage(pmask, paste(dir, "/pbodymask.png", sep=""))
    writeImage(pmask4, paste(dir,  "/pbodymaskdn.png", sep=""))
    writeImage(pmask6, paste(dir,  "/pbodyaddmask.tif", sep=""))
    writeImage(rgbImage(cad, cad, cad), paste(dir,  "/cad.png", sep=""))
    writeImage(cadblur, paste(dir,  "/cadblur.png", sep=""))
    writeImage(cmask3, paste(dir, "/cadmask3.png", sep=""))
    writeImage(cmask4, paste(dir,  "/cadmask4.png", sep=""))
    writeImage(cmask5, paste(dir,  "/cadmask5.png", sep=""))
    writeImage(normalize(cmask10), paste(dir, "/cadmask10.png", sep=""))
    writeImage(normalize(cad + (cmask9>0)*0.3), paste(dir, "/cadmask9m.png", sep=""))
    writeImage(normalize(cad + (cmask10>0)*0.3), paste(dir,  "/cadmask9m2.png", sep=""))
    writeImage(cadc + (cmask12>0)*0.3, paste(dir,  "/cadmask12m.png", sep=""))
    writeImage(germpbody, paste(dir, "/", "pbodyfinal", ".png", sep=""))
    writeImage(mergest, paste(dir, "/",  "mergest", ".png", sep=""))
    
    npbody <- max(bwlabel(germpbody))
    npbody_per_germ <- max(bwlabel(germpbody))/ngerm
    ifelse(is.null(germpbodyintav), intpbody <- NA, intpbody <- germpbodyintav)
    ifelse(is.null(germintav), intgerm <- NA, intgerm <- germintav[[1]])
    
  }
  result <- data.frame(filename=filelist, 
                       pbody = unlist(npbody),
                       mol = unlist(nmol),
                       nmolpbody = unlist(nmolpbody),
                       germ=unlist(ngerm), 
                       pbody_per_germ= unlist(npbody_per_germ),
                       mol_per_germ= unlist(nmol_per_germ),
                       mol_pbody_per_germ = unlist(nmolpbody_per_germ),
                       intmolpbody= unlist(intmolpbody),
                       intpbody = unlist(intpbody),
                       germintav = unlist(intgerm))
  
  return(list(result, cells))

ptm <- proc.time()
# Control vs Tg +/- series
result1 <- pbodynum(dir1, mol="nanos2", molthresh=0.05, pbodythresh=0.05, cadch=3, pbodych=2, molch=1)
saveRDS(result1, file=paste(dir1, "3-5_048.rds", sep=""))
result2 <- pbodynum(dir2, mol="nanos2", molthresh=0.08, pbodythresh=0.05, cadch=3, pbodych=2, molch=1)
saveRDS(result2, file=paste(dir2, "3-5_048.rds", sep=""))
result3 <- pbodynum(dir3, mol="nanos2", molthresh=0.08, pbodythresh=0.05, cadch=3, pbodych=2, molch=1)
saveRDS(result3, file=paste(dir3, "3-5_048.rds", sep=""))
result4 <- pbodynum(dir4, mol="nanos2", molthresh=0.05, pbodythresh=0.05, cadch=3, pbodych=2, molch=1)
saveRDS(result4, file=paste(dir4, "3-5.rds", sep=""))
result5 <- pbodynum(dir5, mol="nanos2", molthresh=0.08, pbodythresh=0.05, cadch=3, pbodych=2, molch=1)
saveRDS(result5, file=paste(dir5, "3-5_048.rds", sep=""))
result6 <- pbodynum(dir6, mol="nanos2", molthresh=0.08, pbodythresh=0.05, cadch=3, pbodych=2, molch=1)
saveRDS(result6, file=paste(dir6, "3-5_048.rds", sep=""))
result14 <- pbodynum(dir14, mol="nanos2", molthresh=0.08, pbodythresh=0.05, cadch=3, pbodych=2, molch=1)
saveRDS(result14, file=paste(dir14, "3-5_048.rds", sep=""))
proc.time() - ptm
result1 <- readRDS(file=paste(dir1, "3-4.rds", sep=""))
result2 <- readRDS(file=paste(dir2, "3-4.rds", sep=""))
result3 <- readRDS(file=paste(dir3, "3-4.rds", sep=""))
result4 <- readRDS(file=paste(dir4, "3-4.rds", sep=""))
result5 <- readRDS(file=paste(dir5, "3-4.rds", sep=""))
result6 <- readRDS(file=paste(dir6, "3-4.rds", sep=""))
result14 <- readRDS(file=paste(dir14, "3-4.rds", sep=""))

# Image-wise analysis
dfcont <- rbind(result1[[1]][1:12,], result2[[1]][1:15,], result3[[1]][1:11,], result5[[1]][1:15,], result6[[1]][1:15,], result14[[1]][1:15,])
dfnoutr <- rbind(result1[[1]][13:24,], result2[[1]][16:40,], result3[[1]][12:21,], result5[[1]][16:30,], result6[[1]][16:28,], result14[[1]][16:30,])
dfutr <- rbind(result1[[1]][25:36,], result2[[1]][41:55,], result3[[1]][22:31,], result5[[1]][31:46,], result14[[1]][31:45,])
dfall <- rbind(dfcont, dfnoutr, dfutr)

# Cell-wise analysis
dfcellcont <- rbind(do.call(rbind.data.frame, result1[[2]][1:12]), do.call(rbind.data.frame, result2[[2]][1:15]), do.call(rbind.data.frame, result3[[2]][1:11]), do.call(rbind.data.frame, result5[[2]][1:15]), do.call(rbind.data.frame, result6[[2]][1:15]), do.call(rbind.data.frame, result14[[2]][1:15]))
dfcellnoutr <- rbind(do.call(rbind.data.frame, result1[[2]][13:24]), do.call(rbind.data.frame, result2[[2]][16:40]), do.call(rbind.data.frame, result3[[2]][12:21]), do.call(rbind.data.frame, result5[[2]][16:30]), do.call(rbind.data.frame, result6[[2]][16:28]), do.call(rbind.data.frame, result14[[2]][16:30]))
dfcellutr <- rbind(do.call(rbind.data.frame, result1[[2]][25:36]), do.call(rbind.data.frame, result2[[2]][41:55]), do.call(rbind.data.frame, result3[[2]][22:31]), do.call(rbind.data.frame, result5[[2]][31:46]), do.call(rbind.data.frame, result14[[2]][31:45]))
dfcellall <- rbind(dfcellcont, dfcellnoutr, dfcellutr)

dfcellcont[,1] <- rep("Control")
dfcellnoutr[,1] <- rep("Dazl3F;Flp")
dfcellutr[,1] <- rep("Dazl3F")
dfcellcont <- cbind(dfcellcont, molperpbody=dfcellcont$pbodymolstnum/dfcellcont$pbodystnum)
dfcellnoutr <- cbind(dfcellnoutr, molperpbody=dfcellnoutr$pbodymolstnum/dfcellnoutr$pbodystnum)
dfcellutr <- cbind(dfcellutr, molperpbody=dfcellutr$pbodymolstnum/dfcellutr$pbodystnum)
dfcellall <- rbind(dfcellcont, dfcellnoutr, dfcellutr)

# Extract outliers for pbody
pbodycontout <- dfcellcont[which(dfcellcont[,"pbodystnum"]>7),]
pbodynoutrout <- dfcellnoutr[which(dfcellnoutr[,"pbodystnum"]>6),]
pbodyutrout <- dfcellutr[which(dfcellutr[,"pbodystnum"]>7),]
pbodyallout <- rbind(pbodycontout, pbodynoutrout, pbodyutrout)

windowsFonts(calibri = windowsFont("Calibri"))
n_fun <- function(x){
  return(data.frame(y = 0, label = paste0("(", length(x), ")")))
}

m <- ggplot() +
  geom_boxplot(data = dfcellall, aes(factor(file), y=pbodystnum), outlier.colour=NA) +
  geom_dotplot(data = pbodyallout, aes(factor(file), y=pbodystnum),  binaxis = "y", 
               stackdir = "center", position = position_jitter(width=0, height = 0.3), 
               binwidth = 1, dotsize=0.1, stackratio=0.2) +
  scale_y_continuous(breaks=seq(0, 12, by=2), name="Number of granules/cell") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, vjust=0.3, hjust = 1, size=18, color=1), 
        axis.title.x=element_text(size=18),
        strip.text.x = element_text(size = 18),
        axis.title.y=element_text(vjust=0.6, size=18),
        axis.text.y=element_text(size=18, color=1))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1)) 
m
ggsave(filename="pbodynumber0219.png", plot=m, width=2, height=7)

# Extract outliers for Nanos2nanos2contout <- dfcellcont[which(dfcellcont[,"molstnum"]>7),]
nanos2noutrout <- dfcellnoutr[which(dfcellnoutr[,"molstnum"]>5),]
nanos2utrout <- dfcellutr[which(dfcellutr[,"molstnum"]>7),]
nanos2allout <- rbind(nanos2contout, nanos2noutrout, nanos2utrout)

m <- ggplot() +
  geom_boxplot(data = dfcellall, aes(factor(file), y=molstnum), outlier.colour=NA) +
  geom_dotplot(data = nanos2allout, aes(factor(file), y=molstnum),  binaxis = "y", 
               stackdir = "center", position = position_jitter(width=0, height = 0.2), 
               binwidth = 1, dotsize=0.1, stackratio=0.4) +
  scale_y_continuous(breaks=seq(0, 18, by=2), name="Number of granules/cell") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, vjust=0.3, hjust = 1, size=18, color=1), 
        axis.title.x=element_text(size=18),
        strip.text.x = element_text(size = 18),
        axis.title.y=element_text(vjust=0.6, size=18),
        axis.text.y=element_text(size=18, color=1))+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1))
m
ggsave(filename="nanos2number0219.png", plot=m, width=2, height=7)

# Statistical tests
install.packages("RcmdrPlugin.EZR")
library(RcmdrPlugin.EZR)
# Image-wise comparison
group <- rep(1:3, c(nrow(dfcont), nrow(dfnoutr), nrow(dfutr)))
apply(dfall[,-1], 2, function(x) Steel.Dwass(x, group)[,"p"])
# Cell-wise comparison
group2 <- rep(1:3, c(nrow(dfcellcont), nrow(dfcellnoutr), nrow(dfcellutr)))
apply(dfcellall[,-1], 2, function(x) Steel.Dwass(x, group2)[,"p"])

# How many cells were analyzed?
nrow(dfcellcont)
nrow(dfcellutr)
nrow(dfcellnoutr)

# What were the mean values?
library(Hmisc)
round(apply(dfcellcont[,-1], 2, smean.sd), 2)
round(apply(dfcellutr[,-1], 2, smean.sd), 2)
round(apply(dfcellnoutr[,-1], 2, smean.sd), 2)


# Nanos2 vs Dazl series
result7 <- pbodynum(dir7, mol="Dazl", molthresh=0.08, pbodythresh=0.05, cadch=1, pbodych=3, molch=2)
saveRDS(result7, file=paste(dir7, "3-4.rds", sep=""))
result7 <- readRDS(file=paste(dir7, "3-4.rds", sep=""))
result8 <- pbodynum(dir8, mol="Dazl", molthresh=0.08, pbodythresh=0.05, cadch=3, pbodych=2, molch=1)
saveRDS(result8, file=paste(dir8, "3-4.rds", sep=""))
result8 <- readRDS(file=paste(dir8, "3-4.rds", sep=""))
result9 <- pbodynum(dir9, mol="nanos2", molthresh=0.08, pbodythresh=0.05, cadch=3, pbodych=2, molch=1)
saveRDS(result9, file=paste(dir9, "3-4.rds", sep=""))
result9 <- readRDS(file=paste(dir9, "3-4.rds", sep=""))
result10 <- pbodynum(dir10, mol="Dazl", molthresh=0.08, pbodythresh=0.05, cadch=3, pbodych=2, molch=1)
saveRDS(result10, file=paste(dir10, "3-4.rds", sep=""))
result10 <- readRDS(file=paste(dir10, "3-4.rds", sep=""))
result11 <- pbodynum(dir11, mol="nanos2", molthresh=0.08, pbodythresh=0.05, cadch=3, pbodych=2, molch=1)
saveRDS(result11, file=paste(dir11, "3-4.rds", sep=""))
result11 <- readRDS(file=paste(dir11, "3-4.rds", sep=""))
result12 <- pbodynum(dir12, mol="Dazl", molthresh=0.08, pbodythresh=0.05, cadch=3, pbodych=2, molch=1)
saveRDS(result12, file=paste(dir12, "3-4.rds", sep=""))
result12 <- readRDS(file=paste(dir12, "3-4.rds", sep=""))
result13 <- pbodynum(dir13, mol="nanos2", molthresh=0.08, pbodythresh=0.05, cadch=3, pbodych=2, molch=1)
saveRDS(result13, file=paste(dir13, "3-4.rds", sep=""))
result13 <- readRDS(file=paste(dir13, "3-4.rds", sep=""))

result7 <- readRDS(file=paste(dir7, "3-4.rds", sep=""))
result8 <- readRDS(file=paste(dir8, "3-4.rds", sep=""))
result9 <- readRDS(file=paste(dir9, "3-4.rds", sep=""))
result10 <- readRDS(file=paste(dir10, "3-4.rds", sep=""))
result11 <- readRDS(file=paste(dir11, "3-4.rds", sep=""))
result12 <- readRDS(file=paste(dir12, "3-4.rds", sep=""))
result13 <- readRDS(file=paste(dir13, "3-4.rds", sep=""))

dfcellnanos2all <- rbind(do.call(rbind.data.frame, result9[[2]][1:45]), do.call(rbind.data.frame, result11[[2]][1:15]), do.call(rbind.data.frame, result13[[2]][1:30]))
dfcelldazlall <- rbind(do.call(rbind.data.frame, result8[[2]][1:44]), do.call(rbind.data.frame, result10[[2]][1:45]), do.call(rbind.data.frame, result12[[2]][1:29]))
dfcellnanos2all <- cbind(dfcellnanos2all, molperpbody=dfcellnanos2all$pbodymolstnum/dfcellnanos2all$pbodystnum)
dfcelldazlall <- cbind(dfcelldazlall, molperpbody=dfcelldazlall$pbodymolstnum/dfcelldazlall$pbodystnum)

wilcox.test(dfcellnanos2all[,2], dfcelldazlall[,2])
wilcox.test(dfcellnanos2all[,3], dfcelldazlall[,3])
wilcox.test(dfcellnanos2all[,4], dfcelldazlall[,4])
wilcox.test(dfcellnanos2all[,5], dfcelldazlall[,5])

# How many cells were analyzed?
nrow(dfcellnanos2all)
nrow(dfcelldazlall)

library(Hmisc)
#pbodystnum
smean.sd(dfcellnanos2all[,2])
smean.sd(dfcelldazlall[,2])
#molstnum
smean.sd(dfcellnanos2all[,3])
smean.sd(dfcelldazlall[,3])
#pbodymolstnum
smean.sd(dfcellnanos2all[,4])
smean.sd(dfcelldazlall[,4])
#molperpbody
smean.sd(dfcellnanos2all[,5])
smean.sd(dfcelldazlall[,5])

ggplot(dfcellall, aes(x=pbodystnum, fill=file)) + geom_histogram(binwidth=.5, alpha=.8, position="dodge")
ggplot(dfcellall, aes(x=molstnum, fill=file)) + geom_histogram(binwidth=.5, alpha=.8, position="dodge")
ggplot(dfcellall, aes(x=pbodymolstnum, fill=file)) + geom_histogram(binwidth=.5, alpha=.8, position="dodge")

# Extract certain images
extractimg <- function(dir, imgtype){
  filelist <- list.files(dir, pattern=".lsm")
  filename <- unlist(lapply(filelist, function(x) strsplit(x, ".lsm")[[1]]))
  for(i in 1:length(filename))
  {
    imgname <- paste(dir, "/", filename[i], "/", imgtype, sep="")
    tmp <- readImage(imgname)
    writeImage(tmp, paste(dir, "/", filename[i], imgtype, sep=""))
  }
}

extractimg(dir1, "nanos2addmask.png")
extractimg(dir1, "pbodyaddmask.png")
extractimg(dir2, "nanos2addmask.png")
extractimg(dir2, "pbodyaddmask.png")
extractimg(dir3, "nanos2addmask.png")
extractimg(dir3, "pbodyaddmask.png")
extractimg(dir4, "nanos2addmask.png")
extractimg(dir4, "pbodyaddmask.png")
extractimg(dir5, "nanos2addmask.png")
extractimg(dir5, "pbodyaddmask.png")
extractimg(dir6, "nanos2addmask.png")
extractimg(dir6, "pbodyaddmask.png")


# Nanos2-pbody from Cont
round(apply(dfcont[,-1], 2, mean, na.rm=TRUE), 2)
round(apply(dfcont[,-1], 2, sd, na.rm=TRUE), 2)

# Nanos2-pbody from UTR-
round(apply(dfnoutr[,-1], 2, mean, na.rm=TRUE), 2)
round(apply(dfnoutr[,-1], 2, sd, na.rm=TRUE), 2)

# Nanos2-pbody from UTR+
round(apply(dfutr[,-1], 2, mean, na.rm=TRUE), 2)
round(apply(dfutr[,-1], 2, sd, na.rm=TRUE), 2)

pbodystat <- function(df, g1s, g1e, g2s, g2e, g3s, g3e){
  library(DTK)
  # pbody number per image
  c <- 2
  vec <- c(df[g1s:g1e,c], df[g2s:g2e,c], df[g3s:g3e,c])
  vec <- vec[!is.na(vec)]
  f1 <- gl.unequal(n=3,k=c(length(df[g1s:g1e,c][!is.na(df[g1s:g1e,c])]),
                           length(df[g2s:g2e,c][!is.na(df[g2s:g2e,c])]),
                           length(df[g3s:g3e,c][!is.na(df[g3s:g3e,c])])))
  TK.result1<-TK.test(x=vec,f=f1,a=0.05)
  
  # germ cell number per image
  c <- 4
  vec <- c(df[g1s:g1e,c], df[g2s:g2e,c], df[g3s:g3e,c])
  vec <- vec[!is.na(vec)]
  f1 <- gl.unequal(n=3,k=c(length(df[g1s:g1e,c][!is.na(df[g1s:g1e,c])]),
                           length(df[g2s:g2e,c][!is.na(df[g2s:g2e,c])]),
                           length(df[g3s:g3e,c][!is.na(df[g3s:g3e,c])])))
  TK.result2<-TK.test(x=vec,f=f1,a=0.05)
  
  # number of pbody per germ
  c <- 6
  vec <- c(df[g1s:g1e,c], df[g2s:g2e,c], df[g3s:g3e,c])
  vec <- vec[!is.na(vec)]
  f1 <- gl.unequal(n=3,k=c(length(df[g1s:g1e,c][!is.na(df[g1s:g1e,c])]),
                           length(df[g2s:g2e,c][!is.na(df[g2s:g2e,c])]),
                           length(df[g3s:g3e,c][!is.na(df[g3s:g3e,c])])))
  TK.result3<-TK.test(x=vec,f=f1,a=0.05)
  TK.result3
  
  # number of Nanos2 that overlap with pbody per germ
  c <- 8
  vec <- c(df[g1s:g1e,c], df[g2s:g2e,c], df[g3s:g3e,c])
  vec <- vec[!is.na(vec)]
  f1 <- gl.unequal(n=3,k=c(length(df[g1s:g1e,c][!is.na(df[g1s:g1e,c])]),
                           length(df[g2s:g2e,c][!is.na(df[g2s:g2e,c])]),
                           length(df[g3s:g3e,c][!is.na(df[g3s:g3e,c])])))
  TK.result4<-TK.test(x=vec,f=f1,a=0.05)
  TK.result4
  
  # Average intensity of Nanos2 granules that overlap with pbody per germ
  c <- 9
  vec <- c(df[g1s:g1e,c], df[g2s:g2e,c], df[g3s:g3e,c])
  vec <- vec[!is.na(vec)]
  f1 <- gl.unequal(n=3,k=c(length(df[g1s:g1e,c][!is.na(df[g1s:g1e,c])]),
                           length(df[g2s:g2e,c][!is.na(df[g2s:g2e,c])]),
                           length(df[g3s:g3e,c][!is.na(df[g3s:g3e,c])])))
  TK.result5<-TK.test(x=vec,f=f1,a=0.05)
  TK.result5 
  
  # Average intensity of Nanos2 that overlap with pbody per germ
  c <- 10
  vec <- c(df[g1s:g1e,c], df[g2s:g2e,c], df[g3s:g3e,c])
  vec <- vec[!is.na(vec)]
  f1 <- gl.unequal(n=3,k=c(length(df[g1s:g1e,c][!is.na(df[g1s:g1e,c])]),
                           length(df[g2s:g2e,c][!is.na(df[g2s:g2e,c])]),
                           length(df[g3s:g3e,c][!is.na(df[g3s:g3e,c])])))
  TK.result6<-TK.test(x=vec,f=f1,a=0.05)
  TK.result6
  
  # Average intensity of Nanos2 per germ
  c <- 11
  vec <- c(df[g1s:g1e,c], df[g2s:g2e,c], df[g3s:g3e,c])
  vec <- vec[!is.na(vec)]
  f1 <- gl.unequal(n=3,k=c(length(df[g1s:g1e,c][!is.na(df[g1s:g1e,c])]),
                           length(df[g2s:g2e,c][!is.na(df[g2s:g2e,c])]),
                           length(df[g3s:g3e,c][!is.na(df[g3s:g3e,c])])))
  TK.result7<-TK.test(x=vec,f=f1,a=0.05)
  TK.result7$f[,"p adj"]
  
  allres <- data.frame("number of pbody per image" = TK.result1$f[,"p adj"],
                       "number of germ cells per image" = TK.result2$f[,"p adj"],
                       "number of pbody per germ" = TK.result3$f[,"p adj"],
                       "number of Nanos2 in pbody" = TK.result4$f[,"p adj"],
                       "intensity of Nanos2 granules in pbody" = TK.result5$f[,"p adj"],
                       "intensity of Nanos2 in pbody" = TK.result6$f[,"p adj"],
                       "intensity of Nanos2 per germ" = TK.result7$f[,"p adj"])
  return(allres)
}

pbodystat(result2[[1]], g1s=1, g1e=12, g2s=13, g2e=24, g3s=25, g3e=36)
pbodystat(dfall, 
          g1s=1, g1e=nrow(dfcont), 
          g2s=nrow(dfcont)+1, g2e=nrow(dfcont)+nrow(dfnoutr), 
          g3s=nrow(dfcont)+nrow(dfnoutr)+1, g3e=nrow(dfcont)+nrow(dfnoutr)+nrow(dfutr))


df1 <- data.frame(cond=factor(rep("nanos2hetero")), size=unlist(nanos2result[[2]][1:10]))
df2 <- data.frame(cond=factor(rep("nanos2homo")), size=unlist(nanos2result[[2]][11:20]))
df3 <- data.frame(cond=factor(rep("dazlhetero")), size=unlist(dazlresult[[2]][1:10]))
df4 <- data.frame(cond=factor(rep("dazlhomo")), size=unlist(dazlresult[[2]][11:20]))
df <- rbind(df1, df2, df3, df4)

library(ggplot2)
png(file="pbodysize.png", width=1000, res=150)
ggplot(df, aes(x=size, fill=cond)) + geom_histogram(binwidth=10, position="dodge") + xlim(0, 100) + scale_fill_manual(values=c("#3399FF", "#99CCFF", "#FF66CC", "#FFBBCC"))
dev.off()
