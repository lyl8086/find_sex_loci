library('segmented')

# read p values.
pvector <- read.table('infile', header=F)
pvector <- pvector$V1
p.sort  <- sort(pvector,decreasing=F)
# transform p values.
o <- -log10(p.sort)
# generate expected p values.
e <- -log10(ppoints(length(o)))
# qqplot.
png("qqplot.png")
qqplot(e,o,pch=19,cex=1, main="",
xlab=expression(Expected~~-log[10](italic(p))),
ylab=expression(Observed~~-log[10](italic(p))),
xlim=c(0,max(e)), ylim=c(0,max(o)))
dev.off()

# arbitrary retain the top 2000 to reduce compute time.
y<-head(o,2000)[2000:1] # observed
z<-head(e,2000)[2000:1] # expected
png("target.png")
plot(z,y)
dev.off()

# fit the observed and expected p values relationships using GLM based on gaussian distribution
# you may also want to try poisson or quasi, the results are similar for our data.

m<-glm(y~z, gaussian(link = "identity"))

# fits regression models with segmented relationships between observed and expected p values.
# after checking the elbow of qqplot from the above figure, the breakpoint should be around 4.
# we thus give a range of starting values for the breakpoints to be estimated, e.g. 4 to 4.5
# we definitely do not want non-sense estimations.
m.seg<-segmented(m, psi=seq(4,4.5,0.1), control=seg.control(fix.npsi=T, n.boot=0, it.max=1000))

# plot the model fitting.
png("model.png")
plot(m.seg, conf.level=0.95, shade=TRUE)
points(m.seg, link=FALSE, col=2)
dev.off()

# plot the model fitting with breakponts.
png("model_breakpoints.png")
plot(z,y)
plot(m.seg,add=TRUE,link=FALSE,lwd=2, conf.level=0.95, shade=TRUE)
points(m.seg,col=2, link=FALSE)
dev.off()

# show the estimation results.
sink("segmented.out")
print(summary(m.seg))
sink()
print(summary(m.seg))

# find the corresponding observed pvalues.
# we will choose the second estimated breakpoint as the critical values
# please choose breakpoint by your hand.

# find the true expected -log10 p value
elbow_e<-e[length(e[which(e>=m.seg$psi[,2][[2]])])]
# the top number of expected p values
num_e<-length(e[which(e>=elbow_e)])
# the observed -log10 pvalue
cat("========================================================================\n")
cat('observed threshold -log10 p value is:', o[num_e], "\n")
cat('observed threshold p value is:', p.sort[num_e], 'total:', num_e, "\n")
# conventional fdr
ajs_p <- p.adjust(p.sort, method="fdr", n=length(p.sort))
ajs_p.len <- length(ajs_p[which(ajs_p<=0.05)])
cat('fdr corrected threshold p value is:', p.sort[ajs_p.len], 'total:', ajs_p.len, "\n")
cat("========================================================================\n")
