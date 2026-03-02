#Functions
#=====================================
convert.mod.to.type=function(model){

  if(model == "M4"){
    return("Both")
  }else if(model == "M3" | model == "M0.1" | model == "M0.2"){
    return("No.Med")
  }else if(model == "M1.1" | model == "M2.2"){
    return("Cis")
  }else if(model == "M1.2" | model == "M2.1"){
    return("Trans")
  }else{
    return("Other")
  }

}



#=====================================

convert.truth.to.type=function(model){

  if(model == "M4"){
    return("Both")
  }else if(model == "M3" | model == "M0"){
    return("No.Med")
  }else if(model == "M1"){
    return("Cis")
  }else if(model == "M2"){
    return("Trans")
  }

}



#=====================================

ind.med.edge=function(adj = NULL, model = NULL){

  if(is.null(adj)){
    switch(model, M0 = {
      return(1)
    }, M1.1 = {
      return(1)
    }, M1.2 = {
      return(1)
    }, M2 = {
      return(1)
    }, M3 = {
      return(0)
    }, M4 = {
      return(1)
    }, Other = {
      return(NA)
    })
  }else{
    if(adj[2,3]+adj[3,2]>=1){
      return(1)
    }else{
      return(0)
    }
  }
}

#==================================
ind.gmac = function(cis.trans = NULL){

  if(sum(cis.trans)>=1){
    return(1)
  }else{
    return(0)
  }
}
#==================================
score.med.edge = function(x){
  if(x[1]==0 & x[2]==0 | x[1]==1 & x[2]==1){
    return(1)
  }else{
    return(0)
  }
}
#==================================
precision.recall = function(x){
  tp = diag(x)[2]
  fn = x[1,2]
  fp = x[2,1]
  recall = tp/(tp+fn)
  precision = tp/(tp+fp)
  return(c(recall=recall, precision=precision))
}

make.binary.table = function(inf.vec, true.vec){
  t1 = matrix(0, nrow = 2, ncol = 2)
  colnames(t1) = row.names(t1) = c(0, 1)
  dat.table = cbind.data.frame(inf = inf.vec, truth = true.vec)

  ind = c(0,1)
  for(j in 1:2){
    for(i in 1:2){
      t1[i,j] =  dim(subset(dat.table, inf == (i-1) & truth == (j-1)))[1]
    }
  }
  return(t1)
}


##############################################
get.avg=function(x.mrgn=NULL, x.mrpc=NULL, x.gmac = NULL, params=NULL, which.param=2, use.metric = "prec_rec",
                 include.gmac = T, include.mrpc = T){

  #sort unique simulation parameter values
  x=as.numeric(sort(unique(as.character(params[,which.param]))))

  if(use.metric == "prec_rec"){
    x.mean.mrgn=x.mean.mrpc=x.mean.gmac=as.data.frame(matrix(0, nrow = length(x), ncol = 3))
    colnames(x.mean.mrgn)=colnames(x.mean.mrpc)=colnames(x.mean.gmac)=c("param.value", "recall", "prec")
  }else{
    x.mean.mrgn=x.mean.mrpc=x.mean.gmac=as.data.frame(matrix(0, nrow = length(x), ncol = 5))
    colnames(x.mean.mrgn)=colnames(x.mean.mrpc)=colnames(x.mean.gmac)=c("param.value", "prec", "recall",
                                                                        "se.prec", "se.recall")
  }



  #calculate mean across signal values
  for(i in 1:length(x)){

    #calc means and sdevs for mrgn (remove method and inf columns)
    x.mean.mrgn[i,1]=x[i]
    idx = which(sapply(params[,which.param], function(x,y,z) all.equal(x,y,check.attributes=FALSE), y = x[i])==TRUE)
    if(use.metric == "prec_rec"){
      if(length(idx)>1){
        t1=make.binary.table(inf.vec = x.mrgn$inf[idx], true.vec = x.mrgn$truth[idx])
        x.mean.mrgn[i,2:3] = precision.recall(t1)
      }else{
        x.mean.mrgn[i,2:3] = c(NA,NA)
      }
    }else{
      if(length(which(params[,which.param]==x[i]))==1){
        x.mean.mrgn[i,2:3]=x.mrgn[which(as.factor(params[,which.param])==x[i]),-c(1:3)]
        x.mean.mrgn[i,4:5]=NA

      }else{
        x.mean.mrgn[i,2:3]=colMeans(x.mrgn[which(as.factor(params[,which.param])==x[i]),-1], na.rm=TRUE)
        x.mean.mrgn[i,4:5]=apply(x.mrgn[which(as.factor(params[,which.param])==x[i]),-1], 2, sd, na.rm=TRUE)
      }
    }



    #calc means and sdevs for mrpc (remove method and inf columns)
    x.mean.mrpc[i,1]=x[i]
    if(include.mrpc == T){
      if(use.metric=="prec_rec"){
        if(length(idx)>1){
          t1 = make.binary.table(inf.vec = x.mrpc$inf[idx], true.vec = x.mrpc$truth[idx])
          x.mean.mrpc[i,2:3] = precision.recall(t1)
        }else{
          x.mean.mrpc[i,2:3] = c(NA,NA)
        }
      }else{
        if(length(which(params[,which.param]==x[i]))==1){
          x.mean.mrpc[i,2:3]=x.mrpc[which(as.factor(params[,which.param])==x[i]),-1]
          x.mean.mrpc[i,4:5]=NA

        }else{
          x.mean.mrpc[i,2:3]=colMeans(x.mrpc[which(as.factor(params[,which.param])==x[i]),-1], na.rm=TRUE)
          x.mean.mrpc[i,4:5]=apply(x.mrpc[which(as.factor(params[,which.param])==x[i]),-1],2,sd, na.rm=TRUE)
        }
      }
    }






    #calc means and sdevs for mrpc (remove method and inf columns)
    x.mean.gmac[i,1]=x[i]
    if(include.gmac == T){
      if(use.metric=="prec_rec"){
        if(length(idx)>1){
          t1=make.binary.table(inf.vec = x.gmac$inf[idx], true.vec = x.gmac$truth[idx])
          x.mean.gmac[i,2:3] = precision.recall(t1)
        }else{
          x.mean.gmac[i,2:3]=c(NA,NA)
        }
      }else{
        stop("Cannot Use GMAC for Edge-based metrics!!!")
      }
    }
  }

  #bind for plotting
  all.means.x=rbind.data.frame(cbind.data.frame(Method=rep("MRGN", length(x)), x.mean.mrgn),
                               cbind.data.frame(Method=rep("MRPC", length(x)), x.mean.mrpc),
                               cbind.data.frame(Method=rep("GMAC", length(x)), x.mean.gmac))
  return(all.means.x)

}
##############################################
plot.spline=function(Data=NULL, pd=NULL, spline.int=NULL, params=NULL, which.param=NULL, rplot=4, nplot=1,
                     brks1=seq(0,1,0.2), brks2=seq(0,1,0.2), xbrks=NULL, use.metric = "prec_rec"){

  if(is.null(xbrks)){
    xbrks=unique(Data$param.value)
  }
  #precision
  if(use.metric == "prec_rec"){
    A=ggplot(data=Data, aes(x=param.value, y=prec, color=Method))+
      geom_line(position = pd)+
      geom_point(position = pd)+
      stat_smooth(n=spline.int)+
      theme_classic()+
      #scale_y_continuous(limits = lims, breaks = brks)+
      scale_y_continuous(breaks = brks1)+
      scale_x_continuous(breaks = xbrks)+
      xlab(paste0(colnames(params)[which.param]))+
      ylab("Precision T1 - T2 edge")
    #Recall
    B=ggplot(data=Data, aes(x=param.value, y=recall, color=Method))+
      geom_line(position = pd)+
      geom_point(position = pd)+
      stat_smooth(n=spline.int)+
      theme_classic()+
      #scale_y_continuous(limits = lims, breaks = brks)+
      scale_y_continuous(breaks = brks2)+
      scale_x_continuous(breaks = xbrks)+
      xlab(paste0(colnames(params)[which.param]))+
      ylab("Recall T1 - T2 edge")
    Z=ggarrange(A,B, labels=c("A","B"), nrow=rplot, ncol=nplot)
    plot(Z)

  }else{
    A=ggplot(data=Data, aes(x=param.value, y=prec, color=Method))+
      geom_line(position = pd)+
      geom_point(position = pd)+
      stat_smooth(n=spline.int)+
      theme_classic()+
      #scale_y_continuous(limits = lims, breaks = brks)+
      scale_y_continuous(breaks = brks1)+
      scale_x_continuous(breaks = xbrks)+
      xlab(paste0(colnames(params)[which.param]))+
      ylab("Mean Edge-Based Precision")
    #Recall
    B=ggplot(data=Data, aes(x=param.value, y=recall, color=Method))+
      geom_line(position = pd)+
      geom_point(position = pd)+
      stat_smooth(n=spline.int)+
      theme_classic()+
      #scale_y_continuous(limits = lims, breaks = brks)+
      scale_y_continuous(breaks = brks2)+
      scale_x_continuous(breaks = xbrks)+
      xlab(paste0(colnames(params)[which.param]))+
      ylab("Mean Edge-Based Recall")
    Z=ggarrange(A,B, labels=c("A","B"), nrow=rplot, ncol=nplot)
    plot(Z)

  }


  # #Specificity
  # C=ggplot(data=Data, aes(x=param.value, y=specificity, color=Method))+
  #   geom_line(position = pd)+
  #   geom_point(position = pd)+
  #   stat_smooth(n=spline.int)+
  #   xlab(paste0(colnames(params)[which.param]))+
  #   ylab("Mean Specificity")
  # #FOR
  # D=ggplot(data=Data, aes(x=param.value, y=FOR, color=Method))+
  #   geom_line(position = pd)+
  #   geom_point(position = pd)+
  #   stat_smooth(n=spline.int)+
  #   xlab(paste0(colnames(params)[which.param]))+
  #   ylab("Mean FOR")
  #Combine
}
##############################################
plot.regular=function(Data=NULL, pd=NULL, params=NULL, which.param=NULL, rplot=4, nplot=1,
                      brks1=seq(0,1,0.2), brks2=seq(0,1,0.2), xbrks=NULL){

  if(is.null(xbrks)){
    xbrks=unique(Data$param.value)
  }
  #Precision
  A=ggplot(data=Data, aes(x=param.value, y=acc, color=Method))+
    geom_line(position = pd)+
    geom_point(position = pd)+
    #scale_y_continuous(limits = lims, breaks = brks)+
    scale_y_continuous(breaks = brks1)+
    scale_x_continuous(breaks = xbrks)+
    theme_classic()+
    xlab(paste0(colnames(params)[which.param]))+
    ylab("Mean Acc of T1 - T2 edge")
  # #Recall
  # B=ggplot(data=Data, aes(x=param.value, y=recall, color=Method))+
  #   geom_line(position = pd)+
  #   geom_point(position = pd)+
  #   #scale_y_continuous(limits = lims, breaks = brks)+
  #   scale_y_continuous(breaks = brks2)+
  #   scale_x_continuous(breaks = xbrks)+
  #   xlab(paste0(colnames(params)[which.param]))+
  #   ylab("Mean Recall")
  # #Specificity
  # C=ggplot(data=Data, aes(x=param.value, y=specificity, color=Method))+
  #   geom_line(position = pd)+
  #   geom_point(position = pd)+
  #   xlab(paste0(colnames(params)[which.param]))+
  #   ylab("Mean Specificity")
  # #FOR
  # D=ggplot(data=Data, aes(x=param.value, y=FOR, color=Method))+
  #   geom_line(position = pd)+
  #   geom_point(position = pd)+
  #   xlab(paste0(colnames(params)[which.param]))+
  #   ylab("Mean FOR")
  #Combine
  #Z=ggarrange(A,B, labels=c("A","B"), nrow=rplot, ncol=nplot)
  #plot(Z)
  plot(A)
}


##############################################
plot.sim.metrics=function(metrics.mrgn=NULL, metrics.mrpc=NULL, metrics.gmac = NULL, params=NULL, which.param=2,
                          by.class=FALSE, which.class=NULL, plot.it=TRUE, spline.it=TRUE, spline.int=10, dodge=0.1,
                          return.means=TRUE, rplot=4, nplot=1, brks1=seq(0,1,0.2), brks2=seq(0,1,0.2),
                          xbrks=NULL, save.plot=FALSE, plot.path="", use.metric = "prec_rec", include.gmac = T,
                          include.mrpc = T){

  library(ggpubr)
  pd=position_dodge(dodge)

  if(by.class==TRUE){
    #calc averages within a certain model class
    metrics.mrgn.m=subset(metrics.mrgn, truth==which.class)
    metrics.mrpc.m=subset(metrics.mrpc, truth==which.class)
    metrics.gmac.m=subset(metrics.gmac, truth==which.class)
    #get mean scores across parameter within class
    avg.scores.all=get.avg(x.mrgn = metrics.mrgn.m,
                           x.mrpc = metrics.mrpc.m,
                           x.gmac = metrics.gmac.m,
                           params = params,
                           use.metrics = use.metrics,
                           which.param = which.param,
                           include.gmac = include.gmac,
                           include.mrpc = include.mrpc)
    if(include.gmac==F){
      avg.scores.all = subset(avg.scores.all, Method == "MRGN" | Method == "MRPC")
    }
    if(include.mrpc==F){
      avg.scores.all = subset(avg.scores.all, Method == "MRGN" | Method == "GMAC")
    }

  }else{
    #get mean scores across parameter
    avg.scores.all=get.avg(x.mrgn = metrics.mrgn,
                           x.mrpc = metrics.mrpc,
                           x.gmac = metrics.gmac,
                           params = params,
                           use.metric = use.metric,
                           which.param = which.param,
                           include.gmac = include.gmac)
    if(include.gmac==F){
      avg.scores.all = subset(avg.scores.all, Method == "MRGN" | Method == "MRPC")
    }
    if(include.mrpc==F){
      avg.scores.all = subset(avg.scores.all, Method == "MRGN" | Method == "GMAC")
    }

  }


  #plotting - splined plots
  if(plot.it==TRUE & spline.it==TRUE){

    #for plot saving
    if(save.plot==TRUE){
      pdf(plot.path, height = 8, width = 12)
      plot.spline(avg.scores.all,
                  pd = pd,
                  spline.int = spline.int,
                  params = params,
                  which.param = which.param,
                  rplot = rplot,
                  nplot = nplot,
                  brks1 = brks1,
                  brks2 = brks2,
                  xbrks = xbrks,
                  use.metric = use.metric)
      dev.off()
    }else{
      plot.spline(avg.scores.all,
                  pd = pd,
                  spline.int = spline.int,
                  params = params,
                  which.param = which.param,
                  rplot = rplot,
                  nplot = nplot,
                  brks1 = brks1,
                  brks2 = brks2,
                  xbrks = xbrks,
                  use.metric = use.metric)
    }

    #plotting - regular plots
  }else if(plot.it==TRUE & spline.it==FALSE){
    if(save.plot==TRUE){
      pdf(plot.path, height = 8, width = 12)
      plot.regular(Data = avg.scores.all,
                   pd = pd,
                   params = params,
                   which.param = which.param,
                   rplot = rplot,
                   nplot = nplot,
                   brks1 = brks1,
                   brks2 = brks2)
      dev.off()
    }else{
      plot.regular(Data = avg.scores.all,
                   pd = pd,
                   params = params,
                   which.param = which.param,
                   rplot = rplot,
                   nplot = nplot,
                   brks1 = brks1,
                   brks2 = brks2)
    }
  }

  #return mean score and sdev data
  if(return.means==TRUE){
    return(avg.scores.all)
  }

}






#function to calculate and plot the Confounder selection performance

calc.css=function(true.datasets, gmac.confs.list, mrgn.confs.list, conf.names,plot.graphs = TRUE){
  recall.gmac = NULL
  recall.mrgn = NULL
  prec.gmac = NULL
  prec.mrgn = NULL
  total.true.confs = NULL
  selected.gmac = NULL
  selected.mrgn = NULL
  gmac.correct = NULL
  mrgn.correct = NULL
  #find the selected and true number of confs and construct CSS
  for(i in 1:length(true.datasets)){
    true.confs = colnames(true.datasets[[i]])[-c(1:6, (dim(true.datasets[[i]])[2]-1):(dim(true.datasets[[i]])[[2]]))]
    total.true.confs[i] = length(true.confs)
    selected.gmac[i] = length(which(gmac.confs.list[[i]][[1]] == 1))
    selected.mrgn[i] = length(colnames(mrgn.confs.list[[i]]))
    gmac.correct[i] = length(na.omit(match(conf.names[which(gmac.confs.list[[i]][[1]] == 1)], true.confs)))
    mrgn.correct[i] = length(na.omit(match(true.confs, colnames(mrgn.confs.list[[i]]))))
    recall.gmac[i] = gmac.correct[i]/total.true.confs[i]
    recall.mrgn[i] = mrgn.correct[i]/total.true.confs[i]
    prec.gmac[i] = gmac.correct[i]/selected.gmac[i]
    prec.mrgn[i] = mrgn.correct[i]/selected.mrgn[i]
  }
  conf.data = cbind.data.frame(GMAC.Recall = recall.gmac, MRGN.MRPC.Recall = recall.mrgn, GMAC.Prec = prec.gmac,
                               MRGN.MRPC.Prec = prec.mrgn, Num.of.True = total.true.confs, Selected.GMAC = selected.gmac,
                               Selected.MRGN.MRPC = selected.mrgn, Correct.GMAC = gmac.correct,
                               Correct.MRGN.MRPC = mrgn.correct)

  #constuct summary table of means
  unq.conf = sort(unique(conf.data$Num.of.True))
  tab = as.data.frame(matrix(0, nrow = length(unq.conf), ncol = ncol(conf.data)))
  colnames(tab) = c("Number of True Confounders", paste0("Mean.",colnames(conf.data)[-5]))
  tab[,1] = unq.conf
  for(i in 1:length(unq.conf)){
    tab[i, -1] = colMeans(conf.data[which(conf.data$Num.of.True==unq.conf[i]),-5], na.rm = T)
  }
  #longdata for last plot A5
  LD = cbind.data.frame(Method = c(rep("GMAC", length(unq.conf)), rep("MRGN/MRPC",length(unq.conf))),
                        rbind(as.matrix(tab[,c(1,2,4),]), as.matrix(tab[,c(1,3,5)])))
  colnames(LD) = c("Method", "Num.of.True", "Mean.Recall", "Mean.Precision")

  #plotting
  if(plot.graphs == TRUE){

    library(ggpubr)
    A1 = ggplot(data = conf.data, aes(x = as.factor(Num.of.True), y = GMAC.Recall))+
      geom_point()+geom_boxplot()+
      scale_y_continuous(breaks=seq(0,1,0.1))+
      theme_classic()+
      xlab("Number of True Confounders")+
      ylab("Recall (GMAC)")
    A2 = ggplot(data = conf.data, aes(x = as.factor(Num.of.True), y = GMAC.Prec))+
      geom_point()+geom_boxplot()+#geom_abline(intercept = 0, slope = 1, colour = "Red")+
      scale_y_continuous(breaks = seq(0,1,.1))+
      theme_classic()+
      xlab("Number of True Confounders")+
      ylab("Precision (GMAC)")
    A3 = ggplot(data = conf.data, aes(x = as.factor(Num.of.True), y = MRGN.MRPC.Recall))+
      geom_point()+geom_boxplot()+
      scale_y_continuous(breaks = seq(0,1,0.1))+
      theme_classic()+
      xlab("Number of True Confounders")+
      ylab("Recall (MRGN/MRPC)")
    A4 = ggplot(data = conf.data, aes(x = as.factor(Num.of.True), y = MRGN.MRPC.Prec))+
      geom_point()+geom_boxplot()+#geom_abline(intercept = 0, slope = 1, colour = "Red")+
      scale_y_continuous(breaks = seq(0,1,.1))+
      theme_classic()+
      xlab("Number of True Confounders")+
      ylab("Precision (MRGN/MRPC)")
    
    if(max(unq.conf > 15)){
      xscale = seq(1,max(unq.conf),5)
    }else{
      xscale = seq(1,max(unq.conf),1)
    }
    A5 = ggplot(data = LD, aes(x = Num.of.True, y = Mean.Recall, color = Method))+
      geom_point(size = 1.5)+geom_line(size = 1)+
      scale_y_continuous(breaks = seq(0.2,1,0.05))+
      scale_x_continuous(breaks = xscale, name = "Number of True Confounders")
    theme_classic()+
      ylab("Mean Recall")
    A6 = ggplot(data = LD, aes(x = Num.of.True, y = Mean.Precision, color = Method))+
      geom_point(size = 1.5)+geom_line(size = 1)+
      scale_y_continuous(breaks = seq(0.2,1,0.1))+
      scale_x_continuous(breaks = xscale, name = "Number of True Confounders")
    theme_classic()+
      ylab("Mean Precision")

    B1 = ggarrange(A1, A2, A3, A4, nrow = 2, ncol = 2)
    #plot(B1)
    B2 = ggarrange(A5, A6, nrow = 1 , ncol = 2, common.legend = T)
    #plot(B2)
    B3 = ggarrange(B1, B2, labels = c("A","B"), nrow = 2 , ncol = 1, heights = c(1.5, 1))
    plot(B3)

  }

  return(list(data = conf.data, table = tab))
}









#function to calculate and plot the Confounder selection performance

calc.css2=function(true.datasets, mrgn.confs.list, conf.names, plot.graphs = TRUE){
  recall.mrgn = NULL
  prec.mrgn = NULL
  total.true.confs = NULL
  selected.mrgn = NULL
  mrgn.correct = NULL
  #find the selected and true number of confs and construct CSS
  for(i in 1:length(true.datasets)){
    true.confs = colnames(true.datasets[[i]])[-c(1:6, (dim(true.datasets[[i]])[2]-1):(dim(true.datasets[[i]])[2]))]
    total.true.confs[i] = length(true.confs)
    selected.mrgn[i] = length(mrgn.confs.list[[i]])
    mrgn.correct[i] = length(na.omit(match(true.confs, mrgn.confs.list[[i]])))
    recall.mrgn[i] = mrgn.correct[i]/total.true.confs[i]
    prec.mrgn[i] = mrgn.correct[i]/selected.mrgn[i]
  }
  conf.data = cbind.data.frame(Recall = recall.mrgn, Prec = prec.mrgn, Num.of.True = total.true.confs, 
                               Selected = selected.mrgn, Correct = mrgn.correct)
  
  #constuct summary table of means
  unq.conf = sort(unique(conf.data$Num.of.True))
  tab = as.data.frame(matrix(0, nrow = length(unq.conf), ncol = ncol(conf.data)))
  colnames(tab) = c("Number of True Confounders", paste0("Mean.",colnames(conf.data)[-3]))
  tab[,1] = unq.conf
  for(i in 1:length(unq.conf)){
    tab[i, -1] = colMeans(conf.data[which(conf.data$Num.of.True==unq.conf[i]),-3], na.rm = T)
  }
  #longdata for last plot A5
  LD = as.data.frame(tab[,c(1,2,3),])
  colnames(LD) = c("Num.of.True", "Mean.Recall", "Mean.Precision")
  
  #plotting
  if(plot.graphs == TRUE){
    
    library(ggpubr)
    library(ggthemes)
    # A1 = ggplot(data = conf.data, aes(x = as.factor(Num.of.True), y = GMAC.Recall))+
    #   geom_point()+geom_boxplot()+
    #   scale_y_continuous(breaks=seq(0,1,0.1))+
    #   theme_classic()+
    #   xlab("Number of True Confounders")+
    #   ylab("Recall (GMAC)")
    # A2 = ggplot(data = conf.data, aes(x = as.factor(Num.of.True), y = GMAC.Prec))+
    #   geom_point()+geom_boxplot()+#geom_abline(intercept = 0, slope = 1, colour = "Red")+
    #   scale_y_continuous(breaks = seq(0,1,.1))+
    #   theme_classic()+
    #   xlab("Number of True Confounders")+
    #   ylab("Precision (GMAC)")
    # A3 = ggplot(data = conf.data, aes(x = as.factor(Num.of.True), y = MRGN.MRPC.Recall))+
    #   geom_point()+geom_boxplot()+
    #   scale_y_continuous(breaks = seq(0,1,0.1))+
    #   theme_classic()+
    #   xlab("Number of True Confounders")+
    #   ylab("Recall (MRGN/MRPC)")
    # A4 = ggplot(data = conf.data, aes(x = as.factor(Num.of.True), y = MRGN.MRPC.Prec))+
    #   geom_point()+geom_boxplot()+#geom_abline(intercept = 0, slope = 1, colour = "Red")+
    #   scale_y_continuous(breaks = seq(0,1,.1))+
    #   theme_classic()+
    #   xlab("Number of True Confounders")+
    #   ylab("Precision (MRGN/MRPC)")
    
    rps = dim(LD)[1]
    LD.resized = cbind.data.frame(Metric = c(rep("Precision", rps),rep("Recall", rps)), 
                                  Num.of.True = rep(LD$Num.of.True, 2),
                                  Performance = c(LD$Mean.Precision, LD$Mean.Recall))
    
    color.codes = c("#0073C2FF", "#EFC000FF", "#868686FF")
    A5 = ggplot(data = LD.resized, aes(x = Num.of.True, y = Performance, color = Metric))+
      geom_point(size = 1.5)+geom_line(size = 1)+
      theme_hc()+
      geom_smooth()+
      scale_fill_manual(values = color.codes)+
      scale_color_manual(values = color.codes)+
      theme(legend.position = "top")+
      scale_y_continuous(breaks = seq(0,1,0.1))+
      scale_x_continuous(breaks = seq(0,50,5), name = "Number of True Confounders")
      ylab("Mean Precision/Recall")
    
    plot(A5)
    # rps2 = dim(tab)[1]
    # LD2 = cbind.data.frame(Metric = c(rep("Avg Confs Selected", rps2), rep("Avg Confs Correct"))
    # A6 = ggplot(data = tab, aes(x = Num.of.True, y = Performance, color = Metric))+
    #   geom_point(size = 1.5)+geom_line(size = 1)+
    #   scale_y_continuous(breaks = seq(0.2,1,0.1))+
    #   scale_x_continuous(breaks = seq(1,15,1), name = "Number of True Confounders")
    # theme_classic()+
    #   ylab("Mean Precision/Recall")
    
    # B1 = ggarrange(A1, A2, A3, A4, nrow = 2, ncol = 2)
    #plot(B1)
    #B2 = ggarrange(A5, A6, nrow = 1 , ncol = 2, common.legend = T)
    #plot(B2)
    # B3 = ggarrange(B1, B2, labels = c("A","B"), nrow = 2 , ncol = 1, heights = c(1.5, 1))
    #plot(B2)
    
  }
  
  return(list(data = conf.data, table = tab))
}




# a function to easily calculate precision and recall for all class-based and edge-based metrics
#---------------function---------------------------
tub.calc.prec.rec = function(reg.res, params, Name){
  
  models = c("M0","M1","M2","M3","M4","Other")
  preds = convert.cats(unlist(reg.res[14,]))
  truth = convert.truth(params$model)
  
  pred.adj = lapply(c(1:300), function(x,y,z) get.adj.from.class(y[x], reg.vec = unlist(z[1:6,x])), 
                    y = unlist(reg.res[14,]), z = reg.res)
  
  true.adj = lapply(c(1:300), function(x,y) get.adj.from.class(y[x]),
                    y = truth)
  
  true.edge.ind = unlist(lapply(true.adj, ind.med.edge))
  pred.edge.ind = unlist(lapply(pred.adj, ind.med.edge))

  all.res = cbind.data.frame(truth = truth, preds = preds)
  cm = matrix(0, nrow = 6, ncol = 5)
  colnames(cm) = models[1:5]
  row.names(cm) = models
  
  for(i in 1:length(models)){
    for(j in 1:length(models[1:5])){
      cm[i, j] = dim(subset(all.res, truth == models[j] & preds == models[i]))[1]
    }
  }
  
  rec = diag(cm)/colSums(cm) 
  prec = diag(cm)/rowSums(cm)[1:5]
  
  edge.prec.rec = colMeans(do.call('rbind', lapply(c(1:300), 
                           function(x,y,z) get.metrics(Truth = convert.truth(many.conf.params$model)[x],
                                                      Inferred = pred.adj[[x]],
                                                      get.adj.truth = TRUE))), na.rm = T)
  
  t1.t2.tab = table(MRGN = pred.edge.ind, Truth = true.edge.ind)
  t1.t2.rec = (diag(t1.t2.tab)/colSums(t1.t2.tab))[2]
  t1.t2.prec = (diag(t1.t2.tab)/rowSums(t1.t2.tab))[2]
  
  
  recall.final = c(rec, t1.t2.rec, edge.prec.rec[2])
  precision.final = c(prec, t1.t2.prec, edge.prec.rec[1])
  
  df = cbind.data.frame(Case = rep(Name, 7), recall = recall.final, precision = precision.final)
  row.names(df) = c(models[1:5], "T1.T2.Edge", "All.Edges")
  
  return(df)
  
}

