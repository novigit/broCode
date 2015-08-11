# get command line arguments
args<-commandArgs(TRUE)

# jobid
jobid_idx<-grep("^jobid=", args)
jobid<-unlist(strsplit(args[jobid_idx], "="))[2]
if (length(jobid_idx) > 0){
  args <- args[-jobid_idx]
}

# node
node_idx<-grep("^node=", args)
node<-unlist(strsplit(args[node_idx], "="))[2]
if (length(node_idx) > 0){
  args <- args[-node_idx]
}

# jobname
jobname_idx<-grep("^jobname=", args)
jobname<-unlist(strsplit(args[jobname_idx], "="))[2]
if (length(jobname_idx) > 0){
  args <- args[-jobname_idx]
}

# output
output_idx<-grep("^output=", args)
output<-unlist(strsplit(args[output_idx], "="))[2]
if (length(output_idx) > 0){
  args <- args[-output_idx]
}

columns<-c("Time","Seconds","GB_LIMIT","GB_USED","GB_SWAP_USED","CPU1","CPU2","CPU3","CPU4","CPU5","CPU6","CPU7","CPU8","CPU9","CPU10","CPU11","CPU12","CPU13","CPU14","CPU15","CPU16")
log=paste('/sw/share/slurm/milou/uppmax_jobstats/', node, '/', jobid, sep="")
df<-read.table(log,skip=1,header=F,col.names=columns)

# convert time from seconds to hours and so that it starts with 0
time<-(df$Seconds-df$Seconds[1])/3600

# get cpu usage
totalcpu<-rowSums(df[,c(6:21)])

# get memory usage
mem<-df$GB_USED

# plot into output
png(filename=output, width=480, height=480,units="px")

par(mar=c(5, 4, 4, 5))
plot(time,mem,type="l", ylab="", xlab="Time (hours)", lwd=1.5, col="blue", las=1, main=jobname)
par(new=T)
plot(time,totalcpu,type="l", axes=F, ylab="", xlab="", xaxt="n", lty=2, lwd=1, col="darkgreen")
axis(4, las=1, pretty(c(0,max(totalcpu))))
mtext("Memory Usage (GB)", side=2,line=3,col="blue")
mtext("CPU Usage (%)", side=4,line=4,col="darkgreen")

dev.off()
