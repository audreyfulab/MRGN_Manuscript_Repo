
library(MRGN)
library(ggthemes)
#needed functions

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
  if(sum(x[2,])>0){
    precision = tp/(tp+fp)
  }else{
    precision = 0
  }
  return(c(recall=recall, precision=precision))
}


#==================================
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
get.avg=function(x.mrgn=NULL, x.mrgn.libconf=NULL, x.mrpc=NULL, x.gmac = NULL, params=NULL, which.param=2, 
                 use.metric = "t1t2", include.gmac = T, include.mrpc = T){
  
  #sort unique simulation parameter values
  x=as.numeric(sort(unique(as.character(params[,which.param]))))
  
  if(use.metric == "t1t2"){
    x.mean.mrgn=x.mean.mrgn.libconf=x.mean.mrpc=x.mean.gmac=as.data.frame(matrix(0, nrow = length(x), ncol = 3))
    colnames(x.mean.mrgn)=colnames(x.mean.mrgn.libconf)=colnames(x.mean.mrpc)=colnames(x.mean.gmac)=c("param.value", 
                                                                                                      "recall", 
                                                                                                      "prec")
  }else{
    x.mean.mrgn=x.mean.mrgn.libconf=x.mean.mrpc=x.mean.gmac=as.data.frame(matrix(0, nrow = length(x), ncol = 5))
    colnames(x.mean.mrgn)=colnames(x.mean.mrgn.libconf)=colnames(x.mean.mrpc)=colnames(x.mean.gmac)=c("param.value", 
                                                                                                      "prec", 
                                                                                                      "recall",
                                                                                                      "se.prec", 
                                                                                                      "se.recall")
  }
  
  
  
  #calculate mean across signal values
  for(i in 1:length(x)){
    
    #calc means and sdevs for mrgn (remove method and inf columns)
    x.mean.mrgn[i,1]=x[i]
    idx = which(sapply(params[,which.param], function(x,y,z) all.equal(x,y,check.attributes=FALSE), y = x[i])==TRUE)
    if(use.metric == "t1t2"){
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
        x.mean.mrgn[i,2:3]=colMeans(x.mrgn[which(as.factor(params[,which.param])==x[i]),-c(1:3)], na.rm=TRUE)
        x.mean.mrgn[i,4:5]=apply(x.mrgn[which(as.factor(params[,which.param])==x[i]),-c(1:3)], 2, sd, na.rm=TRUE)
      }
    }
    
    
    #calc means and sdevs for mrgn liberal.conf (remove method and inf columns)
    x.mean.mrgn.libconf[i,1]=x[i]
    idx = which(sapply(params[,which.param], function(x,y,z) all.equal(x,y,check.attributes=FALSE), y = x[i])==TRUE)
    if(use.metric == "t1t2"){
      if(length(idx)>1){
        t1=make.binary.table(inf.vec = x.mrgn.libconf$inf[idx], true.vec = x.mrgn.libconf$truth[idx])
        x.mean.mrgn.libconf[i,2:3] = precision.recall(t1)
      }else{
        x.mean.mrgn.libconf[i,2:3] = c(NA,NA)
      }
    }else{
      if(length(which(params[,which.param]==x[i]))==1){
        x.mean.mrgn.libconf[i,2:3]=x.mrgn.libconf[which(as.factor(params[,which.param])==x[i]),-c(1:3)]
        x.mean.mrgn.libconf[i,4:5]=NA
        
      }else{
        x.mean.mrgn.libconf[i,2:3]=colMeans(x.mrgn.libconf[which(as.factor(params[,which.param])==x[i]),-c(1:3)], na.rm=TRUE)
        x.mean.mrgn.libconf[i,4:5]=apply(x.mrgn.libconf[which(as.factor(params[,which.param])==x[i]),-c(1:3)], 2, sd, na.rm=TRUE)
      }
    }
    
    
    
    #calc means and sdevs for mrpc (remove method and inf columns)
    x.mean.mrpc[i,1]=x[i]
    if(include.mrpc == T){
      if(use.metric=="t1t2"){
        if(length(idx)>1){
          t1 = make.binary.table(inf.vec = x.mrpc$inf[idx], true.vec = x.mrpc$truth[idx])
          x.mean.mrpc[i,2:3] = precision.recall(t1)
        }else{
          x.mean.mrpc[i,2:3] = c(NA,NA)
        }
      }else{
        if(length(which(params[,which.param]==x[i]))==1){
          x.mean.mrpc[i,2:3]=x.mrpc[which(as.factor(params[,which.param])==x[i]),-c(1:3)]
          x.mean.mrpc[i,4:5]=NA
          
        }else{
          x.mean.mrpc[i,2:3]=colMeans(x.mrpc[which(as.factor(params[,which.param])==x[i]),-c(1:3)], na.rm=TRUE)
          x.mean.mrpc[i,4:5]=apply(x.mrpc[which(as.factor(params[,which.param])==x[i]),-c(1:3)],2,sd, na.rm=TRUE)
        }
      }
    }
    
    
    
    
    
    
    #calc means and sdevs for mrpc (remove method and inf columns)
    x.mean.gmac[i,1]=x[i]
    if(include.gmac == T){
      if(use.metric=="t1t2"){
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
  all.means.x=rbind.data.frame(cbind.data.frame(Method=rep("MRGN + CSFDR", length(x)), x.mean.mrgn),
                               cbind.data.frame(Method=rep("MRGN + CSnoFDR", length(x)), x.mean.mrgn.libconf),
                               cbind.data.frame(Method=rep("MRPC-ADDIS + CSFDR", length(x)), x.mean.mrpc),
                               cbind.data.frame(Method=rep("GMAC", length(x)), x.mean.gmac))
  return(all.means.x)
  
}


##############################################
plot.spline=function(Data=NULL, pd=NULL, spline.int=NULL, params=NULL, which.param=NULL, rplot=4, nplot=1,
                     brks1=seq(0,1,0.2), brks2=seq(0,1,0.2), xbrks=NULL, use.metric = "t1t2", sp.method = "loess",
                     color.codes = NULL, lims1 = c(0,1), lims2 = c(0,1), remove.y = FALSE){
  
  if(is.null(xbrks)){
    xbrks=unique(Data$param.value)
  }
  #precision
  if(use.metric == "t1t2"){
    A=ggplot(data=Data, aes(x=param.value, y=prec, color=Method))+
      geom_line(position = pd)+
      geom_point(position = pd)+
      stat_smooth(aes(linetype = Method), method = sp.method, n=spline.int)+
      theme_hc()+
      theme(legend.position="bottom", legend.title = element_blank(),
            legend.text = element_text(size = 16),
            axis.text = element_text(size = 16),
            axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
            axis.title.x = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 18),
            legend.spacing.y = unit(1.0, 'cm'))+
      #color_palette(palette = "jco")+
      scale_fill_manual(values = color.codes)+
      scale_color_manual(values = color.codes)+
      #scale_y_continuous(limits = lims, breaks = brks)+
      scale_y_continuous(limits = lims1, breaks = brks1)+
      scale_x_continuous(breaks = xbrks)+
      xlab(paste0(colnames(params)[which.param]))+
      ylab("Precision")
    #Recall
    B=ggplot(data=Data, aes(x=param.value, y=recall, color=Method))+
      geom_line(position = pd)+
      geom_point(position = pd)+
      stat_smooth(aes(linetype = Method), method = sp.method, n=spline.int)+
      theme_hc()+
      theme(legend.position="bottom", legend.title = element_blank(),
            legend.text = element_text(size = 16),
            axis.text = element_text(size = 16),
            axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
            axis.title.x = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 18),
            legend.spacing.y = unit(1.0, 'cm'))+
      #color_palette(palette = "jco")+
      scale_fill_manual(values = color.codes)+
      scale_color_manual(values = color.codes)+
      scale_y_continuous(limits = lims2, breaks = brks2)+
      scale_x_continuous(breaks = xbrks)+
      xlab(paste0(colnames(params)[which.param]))+
      ylab("Recall")
    lgd.grob = ggpubr::get_legend(A, position = "bottom")
    if(remove.y == TRUE){
      Z=ggarrange(A+ rremove("xlab")+rremove("ylab"),B+ rremove("xlab")+rremove("ylab"), nrow=rplot, ncol=nplot, legend = "none",
                  font.label = list(size = 16, face = "bold", color = "black"))
      Z = annotate_figure(Z, bottom = textGrob(paste0(colnames(params)[which.param]), gp = gpar(cex = 1, fontsize = 18)))
    }else{
      Z=ggarrange(A+ rremove("xlab"),B+ rremove("xlab"), nrow=rplot, ncol=nplot, legend = "none", 
                  font.label = list(size = 16, face = "bold", color = "black"))
      Z = annotate_figure(Z, bottom = textGrob(paste0(colnames(params)[which.param]), gp = gpar(cex = 1, fontsize = 18)))
    }
    
    plot(Z)
    return(list(plot = Z, lgd.grob = lgd.grob))
    
  }else{
    A=ggplot(data=Data, aes(x=param.value, y=prec, color=Method))+
      geom_line(position = pd)+
      geom_point(position = pd)+
      stat_smooth(aes(linetype = Method), method = sp.method, n=spline.int)+
      theme_hc()+
      theme(legend.position="bottom", legend.title = element_blank(),
            legend.text = element_text(size = 16),
            axis.text = element_text(size = 16),
            axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
            axis.title.x = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 18),
            legend.spacing.y = unit(1.0, 'cm'))+
      scale_fill_manual(values = color.codes)+
      scale_color_manual(values = color.codes)+
      scale_y_continuous(limits = lims1, breaks = brks1)+
      scale_x_continuous(breaks = xbrks)+
      xlab(paste0(colnames(params)[which.param]))+
      ylab("Precision")
    #Recall
    B=ggplot(data=Data, aes(x=param.value, y=recall, color=Method))+
      geom_line(position = pd)+
      geom_point(position = pd)+
      stat_smooth(aes(linetype = Method), method = sp.method, n=spline.int)+
      theme_hc()+
      theme(legend.position="bottom", legend.title = element_blank(),
            legend.text = element_text(size = 16),
            axis.text = element_text(size = 16),
            axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
            axis.title.x = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 18),
            legend.spacing.y = unit(1.0, 'cm'))+
      scale_fill_manual(values = color.codes)+
      scale_color_manual(values = color.codes)+
      scale_y_continuous(limits = lims2, breaks = brks2)+
      scale_x_continuous(breaks = xbrks)+
      xlab(paste0(colnames(params)[which.param]))+
      ylab("Recall")
    lgd.grob = ggpubr::get_legend(A, position = "bottom")
    if(remove.y == TRUE){
      Z=ggarrange(A+ rremove("xlab")+rremove("ylab"),B+ rremove("xlab")+rremove("ylab"), nrow=rplot, ncol=nplot, legend = "none",
                  font.label = list(size = 16, face = "bold", color = "black"))
      Z = annotate_figure(Z, bottom = textGrob(paste0(colnames(params)[which.param]), gp = gpar(cex = 1, fontsize = 18)))
    }else{
      Z=ggarrange(A+ rremove("xlab"),B+ rremove("xlab"), nrow=rplot, ncol=nplot, legend = "none",
                  font.label = list(size = 16, face = "bold", color = "black"))
      Z = annotate_figure(Z, bottom = textGrob(paste0(colnames(params)[which.param]), gp = gpar(cex = 1, fontsize = 18)))
    }
    plot(Z)
    return(list(plot = Z, lgd.grob = lgd.grob))
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
                      brks1=seq(0,1,0.2), brks2=seq(0,1,0.2), xbrks=NULL, lims = c(0,1)){
  
  if(is.null(xbrks)){
    xbrks=unique(Data$param.value)
  }
  #Precision
  A=ggplot(data=Data, aes(x=param.value, y=acc, color=Method))+
    geom_line(position = pd)+
    geom_point(position = pd)+
    #scale_y_continuous(limits = lims, breaks = brks)+
    scale_y_continuous(limits = lims, breaks = brks1)+
    scale_x_continuous(breaks = xbrks)+
    theme_hc()+
    theme(legend.position="bottom", legend.title = element_blank(),
          legend.text = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
          axis.title.x = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 18),
          legend.spacing.y = unit(1.0, 'cm'))+
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
plot.sim.metrics=function(metrics.mrgn=NULL, metrics.mrgn.libconf = NULL, metrics.mrpc=NULL, metrics.gmac = NULL, 
                          params=NULL, which.param=2, by.class=FALSE, which.class=NULL, plot.it=TRUE, spline.it=TRUE,
                          spline.int=10, dodge=0.1, return.means=TRUE, rplot=4, nplot=1, brks1=seq(0,1,0.2), 
                          brks2=seq(0,1,0.2), lmts1=c(0,1), lmts2 = c(0,1), xbrks=NULL, save.plot=FALSE, plot.path="", 
                          use.metric = "t1t2", include.gmac = T, include.mrpc = T, sp.method = "loess", 
                          remove.ylab = FALSE){
  
  library(ggpubr)
  pd=position_dodge(dodge)
  
  if(by.class==TRUE){
    #calc averages within a certain model class
    metrics.mrgn.m=subset(metrics.mrgn, truth==which.class)
    metrics.mrgn.libconf.m=subset(metrics.mrgn.libconf, truth==which.class)
    metrics.mrpc.m=subset(metrics.mrpc, truth==which.class)
    metrics.gmac.m=subset(metrics.gmac, truth==which.class)
    #get mean scores across parameter within class
    avg.scores.all=get.avg(x.mrgn = metrics.mrgn.m,
                           x.mrgn.libconf = metrics.mrgn.libconf.m,
                           x.mrpc = metrics.mrpc.m,
                           x.gmac = metrics.gmac.m,
                           params = params,
                           use.metrics = use.metrics,
                           which.param = which.param,
                           include.gmac = include.gmac,
                           include.mrpc = include.mrpc)
    avg.scores.all$Method = factor(avg.scores.all$Method, levels = c("MRGN + CSnoFDR", "MRGN + CSFDR", 
                                                                     "GMAC", "MRPC-ADDIS + CSFDR"))
    
    if(include.gmac==F){
      avg.scores.all = subset(avg.scores.all, Method != 'GMAC')
    }
    if(include.mrpc==F){
      avg.scores.all = subset(avg.scores.all, Method != 'MRPC-ADDIS + CSFDR')
    }
    
  }else{
    
    #get mean scores across parameter
    avg.scores.all=get.avg(x.mrgn = metrics.mrgn,
                           x.mrgn.libconf = metrics.mrgn.libconf,
                           x.mrpc = metrics.mrpc,
                           x.gmac = metrics.gmac,
                           params = params,
                           use.metric = use.metric,
                           which.param = which.param,
                           include.gmac = include.gmac)
    avg.scores.all$Method = factor(avg.scores.all$Method, levels = c("MRGN + CSnoFDR", "MRGN + CSFDR", 
                                                                     "GMAC", "MRPC-ADDIS + CSFDR"))
    pal_col = c("#0073C2FF","#7AA6DCFF","#EFC000FF","#868686FF")
    
    if(include.gmac==F){
      avg.scores.all = subset(avg.scores.all, Method != 'GMAC')
      pal_col = pal_col[c(1:2, 4)]
    }
    if(include.mrpc==F){
      avg.scores.all = subset(avg.scores.all, Method != 'MRPC-ADDIS + CSFDR')
      pal_col = pal_col[c(1:3)]
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
                  lims1 = lmts1,
                  lims2 = lmts2,
                  use.metric = use.metric,
                  color.codes = pal_col,
                  remove.y = remove.ylab)
      dev.off()
    }else{
      duel_plots = plot.spline(avg.scores.all,
                               pd = pd,
                               spline.int = spline.int,
                               params = params,
                               which.param = which.param,
                               rplot = rplot,
                               nplot = nplot,
                               brks1 = brks1,
                               brks2 = brks2,
                               xbrks = xbrks,
                               lims1 = lmts1,
                               lims2 = lmts2,
                               use.metric = use.metric,
                               sp.method = sp.method,
                               color.codes = pal_col,
                               remove.y = remove.ylab)
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
                   brks2 = brks2,
                   lims = lmts)
      dev.off()
    }else{
      plot.regular(Data = avg.scores.all,
                   pd = pd,
                   params = params,
                   which.param = which.param,
                   rplot = rplot,
                   nplot = nplot,
                   brks1 = brks1,
                   brks2 = brks2,
                   lims = lmts)
    }
  }
  
  #return mean score and sdev data
  if(return.means==TRUE){
    return(list(avg.data = avg.scores.all, plot = duel_plots$plot, legend = duel_plots$lgd.grob))
  }
  
}

