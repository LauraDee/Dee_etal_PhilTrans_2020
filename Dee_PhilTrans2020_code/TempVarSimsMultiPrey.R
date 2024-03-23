library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(future)
library(future.apply)

# Load simulation functions
source("TempVarSimFunctions.R")
source("TempVarSimFunctionsMultiPrey.R")
source("TempVarPlotFunctions.R")


# Load simulation configurations
sim.configs = fread("./simconfigs/multiprey_configs.csv")


# Run through each config. We'll essentially do a 3x2 set of experiments, for each one changing
# the temperature regime. The 3x2 design is:
# - Type of temperature change: variability only, warming only, or both.
# - TPC: same tpc; prey shifted right and pred shifted left

# We'll examine a range of variability to show the dependency

trange.seq = seq(from=0, to=8, by=0.5)
n.exper = length(trange.seq)

expers.to.run = 1#1:4
plan(multiprocess, workers=1)#4)
all.exp.res = future_lapply(expers.to.run, FUN=function(config.row) {
  source("TempVarSimFunctionsMultiPrey.R")
  source("TempVarPlotFunctions.R")
  sim.specs = do.call('rbind', replicate(n=n.exper, expr=sim.configs[config.row,], simplify=F))[,tamp:=trange.seq]
  exp.res = RunExperimentMultiPrey(sim.specs, return.tseries = T)
  #exp.res[,trange:=trange.seq]
  plots = PlotExperimentPair(results = exp.res[c(2,17)],
                             configs = sim.specs,
                             pred.color = "black",
                             prey.color = "red")
  #return(plots)
  return(exp.res)
})
future:::ClusterRegistry("stop")

sim.data = melt(as.data.table(all.exp.res[[1]][[17]]$sim.res), id.vars=c("time"))
ggplot(sim.data, aes(x=time, y=value, group=variable, colour=variable)) + geom_line()



# Reproduce Figure 2 in the main text but where there is a second prey species with TPC at a fixed position
sim.configs = fread("./simconfigs/heatmap_configs.csv")




RunHeatmapSimsMultiPrey = function(sim.configs, row, n.procs = 6) {
  
  heatmap.configs = GenerateHeatmapConfigsForRow(sim.configs, row)
  
  # Add a second prey species with growth TPC peaking at environmental mean temp
  heatmap.configs[,PreyGrowthAtTemp:=paste0(PreyGrowthAtTemp,";MakeTPCAG(",mean.temp, ",",mean.temp+1,",1.184321,20)")]
  heatmap.configs[,PreyCapAtTemp:=paste0(PreyCapAtTemp,";MakeTPCConst(1)")]
  
  # initialize predator and first prey to values they are in main paper fig 2
  heatmap.configs$prey.pop.init = heatmap.configs$pred.mort.at.topt/(heatmap.configs$pred.attack.at.topt*(heatmap.configs$pred.eff.at.topt-heatmap.configs$pred.mort.at.topt*heatmap.configs$pred.handle.at.topt))
  heatmap.configs$pred.pop.init = heatmap.configs$prey.growth.at.topt/heatmap.configs$pred.attack.at.topt * (1-heatmap.configs$prey.pop.init/heatmap.configs$prey.cap.at.topt) * (heatmap.configs$pred.attack.at.topt*heatmap.configs$pred.handle.at.topt*heatmap.configs$prey.pop.init + 1)

  # add remainder of second prey parameters
  heatmap.configs[,prey.growth.at.topt:=paste0(prey.growth.at.topt,";0.6")]
  heatmap.configs[,prey.cap.at.topt:=paste0(prey.cap.at.topt,";20")]
  
    
  # initialize second prey to same level as first
  heatmap.configs[,prey.pop.init:=paste0(prey.pop.init,";",prey.pop.init)]
  
  indices.by.proc = split(1:nrow(heatmap.configs), 1:n.procs)
  plan(multiprocess, workers=n.procs)
  st=Sys.time()
  heatmap.exp.res = future_lapply(indices.by.proc, FUN=function(indices) {
    source("TempVarSimFunctions.R")
    source("TempVarPlotFunctions.R")
    exp.res = RunExperimentMultiPrey(heatmap.configs[indices,])
    return(list(config=heatmap.configs[indices,],
                results=exp.res))
  })
  en=Sys.time()
  future:::ClusterRegistry("stop")
  all.heatmap.res = rbindlist(lapply(1:n.procs, FUN=function(proc) { 
    # the call to rbindlist will create multiple rows per config: one per each
    # prey species. Code after casts to wide format, combines, so that we get
    # different columns for summary stats for each simulation configuration.
    tmp = rbindlist(heatmap.exp.res[[proc]]$results)
    prey.cols = colnames(tmp)[grep("prey", colnames(tmp))]
    pred.cols = colnames(tmp)[grep("pred", colnames(tmp))]
    tmp[,prey.num:=rep(1:2,times=length(heatmap.exp.res[[proc]]$results))]
    tmp[,config.num:=rep(1:length(heatmap.exp.res[[proc]]$results),each=2)]
    prey.data.wide = dcast(tmp, config.num ~ prey.num, value.var=prey.cols)
    pred.data = unique(tmp[,pred.cols, with=F])
    result.data = cbind(pred.data, prey.data.wide[,config.num:=NULL])
  }))
  all.heatmap.config = rbindlist(lapply(1:n.procs, FUN=function(proc) { heatmap.exp.res[[proc]]$config}))
  all.heatmap.res = cbind(all.heatmap.config, all.heatmap.res)
  fwrite(x=all.heatmap.res, file=paste0("../results/multiprey_heatmap_res",row,".csv"))
}

RunHeatmapSimsMultiPrey(sim.configs, 1, n.procs = 6)
RunHeatmapSimsMultiPrey(sim.configs, 2, n.procs = 6)
RunHeatmapSimsMultiPrey(sim.configs, 3, n.procs = 6)

#note with multiple prey, there will be multiple result rows per config
row = 1
config.dat = fread(paste0("../results/multiprey_heatmap_res1.csv"))


MakeThreePanelHeatmapPlotMultiPrey = function(fname1, fname2, fname3, out.file) {
  
  heatmap.res.1 = fread(fname1)
  heatmap.res.2 = fread(fname2)
  heatmap.res.3 = fread(fname3)
  pred.range = range(c(heatmap.res.1$pred.end.mean, 
                       heatmap.res.2$pred.end.mean, 
                       heatmap.res.3$pred.end.mean))
  col.range = c(0, pred.range[2])
  
  
  heatmap.1 = MakeHeatmapPlotMultiPrey(heatmap.res.1, col.lim = col.range) +
    geom_text(x=0, y=1.5, label="Unstable Coexistence", colour="white", angle=90) +
    geom_text(x=1.2, y=4.0, label="Stable Coexistence", colour="white") +
    geom_text(x=2.4, y=0.8, label="Collapse", colour="white", angle=90)
  
  heatmap.2 = MakeHeatmapPlotMultiPrey(heatmap.res.2, col.lim = col.range) +
    geom_text(x=-0.5, y=1.2, label="Stable Coexistence", colour="white") +
    geom_text(x=2.1, y=1.2, label="Collapse", colour="white")
  
  heatmap.3 = MakeHeatmapPlotMultiPrey(heatmap.res.3, col.lim = col.range) + 
    geom_text(x=0, y=2.5, label="Unstable Coexistence", colour="white", angle=90) +
    geom_text(x=1.9, y=2.1, label="Stable Coexistence", angle=75, colour="white") +
    geom_text(x=2.4, y=0.6, label="Collapse", angle=90, colour="white")
  
  
  pspace.plots = arrangeGrob(LabelPanel(heatmap.1,"a"), 
                             LabelPanel(heatmap.2,"b"), 
                             LabelPanel(heatmap.3,"c"),
                             get_legend(heatmap.1),
                             textGrob("Base Scenario"),
                             textGrob("Efficiency=0.1"),
                             textGrob("Max attack=0.5"),
                             layout_matrix = rbind(c(5,6,7),
                                                   c(1,2,3),
                                                   c(4,4,4)),
                             heights=c(0.1, 0.75, 0.15)
  )
  
  png(out.file, width=900, height=300)
  grid.draw(pspace.plots)
  dev.off()
  grid.draw(pspace.plots)
}


MakeHeatmapPlotMultiPrey = function(heatmap.res, col.lim=c(0,5)) {
  heatmap.res[,min.cv:=pmin(pred.end.cv, prey.end.cv_1, prey.end.cv_2)]
  ggplot(heatmap.res, aes(x=topt.shift.prey, y=tamp, fill=pred.end.mean)) + 
    geom_tile() + 
    stat_contour(aes(z=pred.end.mean), breaks=c(1E-3), colour="white") +
    stat_contour(aes(#z=pred.end.cv/pred.start.cv * pred.end.cv),
      z=min.cv), 
      breaks=c(0.05), 
      colour="white")+#,
    #fill="gray",
    #alpha=0.1,
    #geom="polygon") + 
    scale_fill_viridis_c(name="Long-run predator population",
                         limits=col.lim
    ) +
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) +
    xlab("Prey TPC offset from mean temperature (°C)") +
    ylab("Temperature amplitude (°C)") +
    theme_bw() + 
    theme(legend.position="bottom")
  
}


MakeThreePanelHeatmapPlotMultiPrey("../results/multiprey_heatmap_res1.csv",
                                   "../results/multiprey_heatmap_res2.csv",
                                   "../results/multiprey_heatmap_res3.csv",
                                   "../results/multiprey_heatmaps_eqtype.png")
