library(reshape2)
library(ggplot2)

data = read.table(file="results/acs_new.txt", sep=",", h=T, stringsAsFactors=F)
acs_ids <- unique(data$id)

# plot individual profiles
plot.prf <- function(id, data) {
    tmp <- data[data$id==id,]
    tmp <- tmp[with(tmp, order(pos)),]
    tmp$sample_num <- 1:600
    tmpm <- melt(tmp, id.vars=c("id", "pos", "sample_num"))
    p<- ggplot(tmpm, aes(sample_num, value))+geom_line(aes(group=variable, color=variable))+ggtitle(id)
    return(p)
}
p <- plot.prf(acs_ids[1], data); p
ggsave(file="results/prf1.pdf")

# aggregate
library(plyr)
plot.aggr <- function(data) {
    tmp <- data
    tmp$sample_num <- 1:600
    tmp <- tmp[,-c(2)]
    tmpm <- melt(tmp, id.vars=c("id", "sample_num"))

    tmpaggr <- ddply(tmpm, .(sample_num, variable), summarize, count=sum(value))
    p<- ggplot(tmpaggr, aes(sample_num, count))+geom_line(aes(group=variable, color=variable))
    return(list(aggr=tmpaggr, plot=p))
}
res<-plot.aggr(data)
res[["plot"]]
ggsave(file="results/prf_aggr.pdf")

prf.aggr <- res[["aggr"]]
skew <- ddply(prf.aggr, .(sample_num), function(df) data.frame(
    ATskew=df[df$variable=="cntA","count"]-df[df$variable=="cntT","count"],
    GCskew=df[df$variable=="cntG","count"]-df[df$variable=="cntC","count"])
    )

skew.m <- melt(skew, id="sample_num")
skew.plot <- ggplot(skew.m, aes(sample_num, value))+geom_line(aes(group=variable, color=variable))
skew.plot
ggsave(file="results/profile_skew.pdf")
