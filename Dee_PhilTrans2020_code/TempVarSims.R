library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(future)
library(future.apply)

# Load simulation functions
source("TempVarSimFunctions.R")
source("TempVarPlotFunctions.R")


# Load simulation configurations
sim.configs = fread("./simconfigs/final_configs.csv")


# Run through each config. We'll essentially do a 3x2 set of experiments, for each one changing
# the temperature regime. The 3x2 design is:
# - Type of temperature change: variability only, warming only, or both.
# - TPC: same tpc; prey shifted right and pred shifted left

# We'll examine a range of variability to show the dependency

trange.seq = seq(from=0, to=8, by=0.5)
n.exper = length(trange.seq)

expers.to.run = 1:4
plan(multiprocess, workers=4)
all.exp.res = future_lapply(expers.to.run, FUN=function(config.row) {
  source("TempVarSimFunctions.R")
  source("TempVarPlotFunctions.R")
  sim.specs = do.call('rbind', replicate(n=n.exper, expr=sim.configs[config.row,], simplify=F))[,tamp:=trange.seq]
  exp.res = RunExperiment(sim.specs, return.tseries = T)
  #exp.res[,trange:=trange.seq]
  plots = PlotExperimentPair(results = exp.res[c(2,17)],
                             configs = sim.specs,
                             pred.color = "black",
                             prey.color = "red")
  return(plots)
  #return(exp.res)
})
future:::ClusterRegistry("stop")





same.tpc.plot = MakeTPCPlot(sim.configs[1,])
offset.tpc.plot = MakeTPCPlot(sim.configs[3,])



main.plot.grid = arrangeGrob(LabelPanel(same.tpc.plot,"a"),
                             LabelPanel(all.exp.res[[1]]$rate.plot, "b"),
                             LabelPanel(all.exp.res[[2]]$pop.plot, "c"),
                             LabelPanel(all.exp.res[[1]]$pop.plot, "d"),
                             LabelPanel(offset.tpc.plot, "e"),
                             LabelPanel(all.exp.res[[3]]$rate.plot, "f"),
                             LabelPanel(all.exp.res[[4]]$pop.plot, "g"),
                             LabelPanel(all.exp.res[[3]]$pop.plot, "h"),
                             get_legend(all.exp.res[[1]]$pop.plot + theme(legend.position="bottom")),
                             textGrob("TPCs", hjust= 0),
                             textGrob("Rates", hjust = 0),
                             textGrob("Population Trajectories"),
                             layout_matrix = rbind(c(10,11,12,12),
                                                   c(1,2,3,4),
                                                   c(5,6,7,8),
                                                   c(9,9,9,9)),
                             heights=c(0.1,0.4,0.4,0.1))

png("../results/mainplotgrid.png", width=800, height=400)
grid.draw(main.plot.grid)
dev.off()



# Now run experiments for indirect effects, where we set one species' 
# TPC to a constant (it is unaffected by temperature). We'll run 4 sims:
# 1. Neither affected, 2. Predator unaffected, 3. Prey unaffected, 4. Both affected
expers.to.run = 5:8
plan(multiprocess, workers=4)
all.exp.res = future_lapply(expers.to.run, FUN=function(config.row) {
  source("TempVarSimFunctions.R")
  source("TempVarPlotFunctions.R")
  sim.specs = do.call('rbind', replicate(n=n.exper, expr=sim.configs[config.row,], simplify=F))[,tamp:=trange.seq]
  exp.res = RunExperiment(sim.specs, return.tseries = T)
  #exp.res[,trange:=trange.seq]
  # plots = PlotExperimentPair(results = exp.res[c(2,17)],
  #                            configs = sim.specs,
  #                            pred.color = "black",
  #                            prey.color = "red")
  return(exp.res)
})
future:::ClusterRegistry("stop")

# make combined plots showing predator, prey trajectories for 4 cases.  
var.index = 10
indirect.plot.colors = c("No temperature effects"="black", 
                         "Direct effect only"="blue", 
                         "Indirect effect only"="red", 
                         "Both effects"="purple")
pop.trajs = data.table(t=all.exp.res[[1]][[var.index]]$sim.res[,1],
                       prey.preyfixed.predfixed=all.exp.res[[1]][[var.index]]$sim.res[,2],
                       prey.preyfixed.predvar=all.exp.res[[2]][[var.index]]$sim.res[,2],
                       prey.preyvar.predfixed=all.exp.res[[3]][[var.index]]$sim.res[,2],
                       prey.preyvar.predvar=all.exp.res[[4]][[var.index]]$sim.res[,2],
                       pred.preyfixed.predfixed=all.exp.res[[1]][[var.index]]$sim.res[,3],
                       pred.preyfixed.predvar=all.exp.res[[2]][[var.index]]$sim.res[,3],
                       pred.preyvar.predfixed=all.exp.res[[3]][[var.index]]$sim.res[,3],
                       pred.preyvar.predvar=all.exp.res[[4]][[var.index]]$sim.res[,3])

# Prey
prey.indirect.plot = ggplot(data=pop.trajs[t<=600,], aes(x=t)) + 
  geom_line(aes(y=prey.preyfixed.predfixed, colour="No temperature effects")) + 
  geom_line(aes(y=prey.preyvar.predfixed, colour="Direct effect only")) + 
  geom_line(aes(y=prey.preyfixed.predvar, colour="Indirect effect only")) + 
  geom_line(aes(y=prey.preyvar.predvar, colour="Both effects")) + 
  scale_y_continuous(expand=c(0,0), limits = c(0,NA)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_colour_manual(name="Experiment",
                      breaks=names(indirect.plot.colors),
                      values=indirect.plot.colors) + 
  xlab("Time") +
  ylab("Prey population") + 
  theme_bw() + 
  theme(legend.position="bottom")

# Predator
pred.indirect.plot = ggplot(data=pop.trajs[t<=600,], aes(x=t)) + 
  geom_line(aes(y=pred.preyfixed.predfixed, colour="No temperature effects")) + 
  geom_line(aes(y=pred.preyfixed.predvar, colour="Direct effect only")) + 
  geom_line(aes(y=pred.preyvar.predfixed, colour="Indirect effect only")) + 
  geom_line(aes(y=pred.preyvar.predvar, colour="Both effects")) + 
  scale_y_continuous(expand=c(0,0), limits = c(0,NA)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_colour_manual(name="Experiment",
                      breaks=names(indirect.plot.colors),
                      values=indirect.plot.colors) + 
  xlab("Time") +
  ylab("Predator population") + 
  theme_bw() + 
  theme(legend.position="bottom")

indirect.grid = grid.arrange(prey.indirect.plot + nolegend,
                             pred.indirect.plot + nolegend,
                             get_legend(prey.indirect.plot),
                             layout_matrix=rbind(c(1,2), c(3,3)),
                             heights=c(0.9, 0.1))

png("../results/indirectplotgrid.png", width=800, height=300)
plot(indirect.grid)
dev.off()



#####################
# Warming scenarios

expers.to.run = 9:20
plan(multiprocess, workers=6)
all.exp.res = future_lapply(expers.to.run, FUN=function(config.row) {
  source("TempVarSimFunctions.R")
  source("TempVarPlotFunctions.R")
  sim.specs = do.call('rbind', replicate(n=n.exper, expr=sim.configs[config.row,], simplify=F))[,tamp:=trange.seq]
  exp.res = RunExperiment(sim.specs, return.tseries = T)
  #exp.res[,trange:=trange.seq]
  # plots = PlotExperimentPair(results = exp.res[c(2,17)],
  #                            configs = sim.specs,
  #                            pred.color = "black",
  #                            prey.color = "red")
  return(exp.res)
})
future:::ClusterRegistry("stop")

# Make a single plot with a) variability only, b) warming only, c) both
# for a high temperature range.
var.only.pop = all.exp.res[[1]][[17]]$sim.res
warm.only.pop = all.exp.res[[2]][[17]]$sim.res
warm.var.pop = all.exp.res[[3]][[17]]$sim.res
warm.inc.var.pop = all.exp.res[[4]][[17]]$sim.res

plot.dt = data.table(cbind(var.only.pop, 
                           warm.only.pop[,2:3], 
                           warm.var.pop[,2:3],
                           warm.inc.var.pop[,2:3]))
setnames(plot.dt, c("time", 
                    "var.prey.pop", "var.pred.pop",
                    "warm.prey.pop", "warm.pred.pop",
                    "warm.var.prey.pop", "warm.var.pred.pop",
                    "warm.inc.var.prey.pop", "warm.inc.var.pred.pop"))

prey.warming.plot = MakeTempScenarioPlot(plot.dt, "prey")
pred.warming.plot = MakeTempScenarioPlot(plot.dt, "pred")


var.only.pop = all.exp.res[[5]][[17]]$sim.res
warm.only.pop = all.exp.res[[6]][[17]]$sim.res
warm.var.pop = all.exp.res[[7]][[17]]$sim.res
warm.inc.var.pop = all.exp.res[[8]][[17]]$sim.res

offset.plot.dt = data.table(cbind(var.only.pop, 
                                  warm.only.pop[,2:3], 
                                  warm.var.pop[,2:3],
                                  warm.inc.var.pop[,2:3]))
setnames(offset.plot.dt, c("time", 
                           "var.prey.pop", "var.pred.pop",
                           "warm.prey.pop", "warm.pred.pop",
                           "warm.var.prey.pop", "warm.var.pred.pop",
                           "warm.inc.var.prey.pop", "warm.inc.var.pred.pop"))

offset.prey.warming.plot = MakeTempScenarioPlot(offset.plot.dt, "prey")
offset.pred.warming.plot = MakeTempScenarioPlot(offset.plot.dt, "pred")


var.only.pop = all.exp.res[[9]][[5]]$sim.res
warm.only.pop = all.exp.res[[10]][[5]]$sim.res
warm.var.pop = all.exp.res[[11]][[5]]$sim.res
warm.inc.var.pop = all.exp.res[[12]][[5]]$sim.res

low.eff.plot.dt = data.table(cbind(var.only.pop, 
                                 warm.only.pop[,2:3], 
                                 warm.var.pop[,2:3],
                                 warm.inc.var.pop[,2:3]))
setnames(low.eff.plot.dt, c("time", 
                    "var.prey.pop", "var.pred.pop",
                    "warm.prey.pop", "warm.pred.pop",
                    "warm.var.prey.pop", "warm.var.pred.pop",
                    "warm.inc.var.prey.pop", "warm.inc.var.pred.pop"))

low.eff.prey.warming.plot = MakeTempScenarioPlot(low.eff.plot.dt, "prey")
low.eff.pred.warming.plot = MakeTempScenarioPlot(low.eff.plot.dt, "pred")



warming.plot.grid = arrangeGrob(#LabelPanel(same.tpc.plot + scale_colour_manual(values=c("Predator"="grey", "Prey"="black")) + rightmargin,"a"),
                                #LabelPanel(offset.tpc.plot + scale_colour_manual(values=c("Predator"="grey", "Prey"="black")) + rightmargin,"d"),
                                LabelPanel(prey.warming.plot + rightmargin,"a"),
                                LabelPanel(pred.warming.plot + rightmargin,"b"),
                                LabelPanel(low.eff.prey.warming.plot + rightmargin,"c"),
                                LabelPanel(low.eff.pred.warming.plot + rightmargin,"d"), 
                                #get_legend(offset.tpc.plot + theme(legend.position="bottom") + scale_colour_manual(name="Species",values=c("Predator"="grey", "Prey"="black"))),
                                get_legend(prey.warming.plot + guides(colour=guide_legend(nrow=2))),
                                #textGrob("TPCs",hjust=-0.05),
                                textGrob("Prey",hjust=0),
                                textGrob("Predator",hjust=0.3),
                                #layout_matrix=rbind(c(9,10,11), c(1,3,4), c(2,5,6), c(7,8,8)),
                                layout_matrix=rbind(c(6,7), c(1,2), c(3,4), c(5,5)),
                                heights=c(0.06, 0.37, 0.37, 0.2),
                                widths=c(0.5, 0.5))

png("../results/warmingplotgrid_new.png", width=800, height=400)
grid.draw(warming.plot.grid)
dev.off()

pdf("../results/warmingplotgrid_new_hires.pdf", width=8, height=4)
grid.draw(warming.plot.grid)
dev.off()


################################################
# simulations for heatmaps; vary both 
# 1) variance of temperature
# 2) offset of TPCs
# Later: shift of temperature mean
sim.configs = fread("./simconfigs/heatmap_configs.csv")



RunHeatmapSims(sim.configs, 1, n.procs = 6)
RunHeatmapSims(sim.configs, 2, n.procs = 6)
RunHeatmapSims(sim.configs, 3, n.procs = 6)


MakeThreePanelHeatmapPlot("../results/heatmap_res1.csv",
                          "../results/heatmap_res2.csv",
                          "../results/heatmap_res3.csv",
                          "../results/heatmaps_eqtype.png")

MakeThreePanelHeatmapPlot("../results/heatmap_res1.csv",
                          "../results/heatmap_res2.csv",
                          "../results/heatmap_res3.csv",
                          "../results/heatmaps_eqtype_hires.pdf",
                          fig.width=12,
                          fig.height=4)

# For understanding, we could produce a similar 3-panel heatmap, but for each location
# in each panel, we'll use a set of fixed parameter values equal to the mean value of the
# parameter in the variable temperature simulations for that location. This is to see if
# the effects we see on populations and stability are driven primarily by the effects of variability
# on average rates or if there are residual effects that depend on the fact that rates are
# time-varying.


RunHeatmapAvgRateSims = function(sim.configs, row, n.procs = 6) {
  
  # generate the configs under variable temperature regimes
  heatmap.configs = GenerateHeatmapConfigsForRow(sim.configs, row)
  # Find the average rates that would arise under the specified temperature regimes.
  # Generate new TPCs with fixed values at those levels
  for(i in 1:nrow(heatmap.configs)) {
    sim.spec = ParseFunctions(as.list(heatmap.configs[i,]))
    temp.fun = MakeSampleTempFun(type=sim.spec$TempAtTime,
                                 mean.temp=sim.spec$mean.temp,
                                 temp.amp=sim.spec$tamp,
                                 warm.rate=sim.spec$tamp/sim.spec$time.end)
    temps = temp.fun(seq(from=sim.spec$time.start, to=sim.spec$time.end, by=sim.spec$time.step))
    heatmap.configs[i,"PredAttackAtTemp"]=paste0("MakeTPCConst(",mean(sim.spec$PredAttackAtTemp(temps)),")")
    heatmap.configs[i,"PreyGrowthAtTemp"]=paste0("MakeTPCConst(",mean(sim.spec$PreyGrowthAtTemp(temps)),")")
  }
  
  # Run simulations with the TPCs fixed at average rates
  indices.by.proc = split(1:nrow(heatmap.configs), 1:n.procs)
  plan(multiprocess, workers=n.procs)
  st=Sys.time()
  heatmap.exp.res = future_lapply(indices.by.proc, FUN=function(indices) {
    source("TempVarSimFunctions.R")
    source("TempVarPlotFunctions.R")
    exp.res = RunExperiment(heatmap.configs[indices,])
    return(list(config=heatmap.configs[indices,],
                results=exp.res))
  })
  en=Sys.time()
  future:::ClusterRegistry("stop")
  all.heatmap.res = rbindlist(lapply(1:n.procs, FUN=function(proc) { rbindlist(heatmap.exp.res[[proc]]$results)}))
  all.heatmap.config = rbindlist(lapply(1:n.procs, FUN=function(proc) { heatmap.exp.res[[proc]]$config}))
  all.heatmap.res = cbind(all.heatmap.config, all.heatmap.res)
  fwrite(x=all.heatmap.res, file=paste0("../results/heatmap_avgrates_res",row,".csv"))
}

RunHeatmapAvgRateSims(sim.configs, 1, n.procs = 6)
RunHeatmapAvgRateSims(sim.configs, 2, n.procs = 6)
RunHeatmapAvgRateSims(sim.configs, 3, n.procs = 6)


MakeThreePanelHeatmapPlot("../results/heatmap_avgrates_res1.csv",
                          "../results/heatmap_avgrates_res2.csv",
                          "../results/heatmap_avgrates_res3.csv",
                          "../results/heatmaps_avgrates_eqtype.png")




##########################
# Direct/indirect effect heatmaps.
# 1. Base case: no temp effects
##########################
RunHeatmapSims(sim.configs = sim.configs, row = 4, n.procs = 6)

# 2. Direct effects on predator only
RunHeatmapSims(sim.configs = sim.configs, row = 5, n.procs = 6)

MakeNetDirectPlot = function(net.eff.file, dir.eff.file, no.eff.file) {
  # load, merge, compute effects
  net.eff.res = fread(net.eff.file)
  dir.eff.res = fread(dir.eff.file)
  no.eff.res = fread(no.eff.file)
  net.eff.res[,min.cv:=pmin(pred.end.cv, prey.end.cv)]
  all.eff.res = merge(net.eff.res[,.(topt.shift.prey=topt.shift.prey, 
                                     tamp=tamp, 
                                     pred.end.net=pred.end.ext.mean,
                                     attack.net=pred.attack.mean,
                                     growth.net=prey.growth.mean)], 
                      dir.eff.res[,.(topt.shift.prey=topt.shift.prey, 
                                     tamp=tamp, 
                                     pred.end.direct=pred.end.ext.mean,
                                     attack.direct=pred.attack.mean)],
                      by=c("topt.shift.prey", "tamp"))
  
  all.eff.res = merge(all.eff.res, 
                      no.eff.res[,.(topt.shift.prey=topt.shift.prey, 
                                    tamp=tamp, 
                                    pred.end.none=pred.end.ext.mean,
                                    attack.none = pred.attack.mean,
                                    growth.none = prey.growth.mean)],
                      by=c("topt.shift.prey", "tamp"))
  
  # compute direct and net effects
  all.eff.res[,`:=`(dir.eff=pred.end.direct-pred.end.none,
                    net.eff=pred.end.net-pred.end.none,
                    attack.rate.eff=attack.net-attack.none,
                    growth.rate.eff=growth.net-growth.none)]
  all.eff.res[,eff.signs:=paste0(net.eff>=0,"_", dir.eff>=0)]
  
  # Net effects with outlines for regions where direct and net effects conflict
  net.eff.plot = ggplot(all.eff.res, aes(x=topt.shift.prey, y=tamp)) + 
    geom_tile(aes(fill=net.eff<0, colour=(dir.eff*net.eff<0), width=0.09, height=0.09), 
              size=0.75, alpha=0.3) + #, size=dir.eff*net.eff<0)) +
    scale_fill_manual(name = "Net effect",
                      labels=c("Positive", "Negative"),
                      values=c("blue", "red"))+
    scale_colour_manual(name="Direct and net effects",
                        labels=c("Same sign", "Opposite sign"),
                        values=c("grey80", "black")) +
    xlab(expression("Prey TPC offset from mean temperature ("*degree*"C)")) + 
    ylab("Temperature amplitude (°C)") +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() + 
    guides(colour=guide_legend(override.aes=list(fill=NA)))# +     
    #stat_contour(data=net.eff.res, aes(z=pred.end.mean), breaks=c(1E-3), colour="white") +
    #stat_contour(data = net.eff.res, aes(#z=pred.end.cv/pred.start.cv * pred.end.cv),
    #  z=min.cv), 
    #  breaks=c(0.05), 
    #  colour="white")
  
  return(net.eff.plot)
}

net.dir.plot.a = MakeNetDirectPlot(net.eff.file = "../results/heatmap_res1.csv",
                                   dir.eff.file = "../results/heatmap_res5.csv",
                                   no.eff.file = "../results/heatmap_res4.csv") + 
                 ggtitle("Net vs. direct effects on predator population")

ggsave(filename="../results/heatmap_netdirect.png", 
       plot = net.dir.plot.a, 
       width=6, height=4, units = "in")

# hi res version
ggsave(filename="../results/heatmap_netdirect_hires.pdf",
       net.dir.plot.a, 
       width=6, height=4, units = "in")

# Redo this, but for parameters matching panels b and c of Fig 2

# Panel b
RunHeatmapSims(sim.configs = sim.configs, row = 6, n.procs = 6) # no temp effects
RunHeatmapSims(sim.configs = sim.configs, row = 7, n.procs = 6) # direct effects on pred only

# Panel c
RunHeatmapSims(sim.configs = sim.configs, row = 8, n.procs = 6) # no temp effects
RunHeatmapSims(sim.configs = sim.configs, row = 9, n.procs = 6) # direct effects on pred only


net.dir.plot.b = MakeNetDirectPlot(net.eff.file = "../results/heatmap_res2.csv",
                                   dir.eff.file = "../results/heatmap_res7.csv",
                                   no.eff.file = "../results/heatmap_res6.csv")

net.dir.plot.c = MakeNetDirectPlot(net.eff.file = "../results/heatmap_res3.csv",
                                   dir.eff.file = "../results/heatmap_res9.csv",
                                   no.eff.file = "../results/heatmap_res8.csv")


net.dir.plot.grid = arrangeGrob(LabelPanel(net.dir.plot.a, "a"), 
                                LabelPanel(net.dir.plot.b, "b"), 
                                LabelPanel(net.dir.plot.c, "c"),
                                get_legend(net.dir.plot.a + theme(legend.position = "bottom")),
                                layout_matrix = rbind(c(1,2,3),
                                                      c(4,4,4)),
                                heights=c(0.9, 0.1))


png("../results/heatmap_netdirect_3panel.png", width=900, height=300, type="cairo")
grid.draw(net.dir.plot.grid)
dev.off()

# Plot rates to help with interpretation

# base rates with no temp variability
ggplot(no.eff.res, aes(x=topt.shift.prey, y=tamp, fill=pred.attack.mean)) + 
  geom_tile() + 
  stat_contour(aes(z=pred.attack.mean), breaks=c(0.0833333333333333), colour="red")

# effects of temp variability 
ggplot(all.eff.res, aes(x=topt.shift.prey, y=tamp, fill=attack.rate.eff)) + 
  geom_tile() + 
  scale_fill_gradient2() +
  stat_contour(aes(z=attack.rate.eff), breaks=c(0), colour="black") + 
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(legend.position = "bottom")

ggplot(all.eff.res, aes(x=topt.shift.prey, y=tamp, fill=growth.rate.eff)) + 
  geom_tile() + 
  scale_fill_gradient2() +
  stat_contour(aes(z=growth.rate.eff), breaks=c(0), colour="black")





##############################
# Redo Fig 2 but with somewhat more realistic temperature series with
# a) non-constant variance,
# b) autoregressive structure, and
# c) seasonality
# d) warming trends
# The baseline scenarios in the main text do c) and d) for expositional purposes
# since it's easy to contrast those extreme cases with identical variance. 
sim.configs = fread("./simconfigs/stoch_configs.csv")

# First, stochastic events

trange.seq = c(0.5, 8)
n.exper = length(trange.seq)

expers.to.run = 1:2
plan(multiprocess, workers=2)
all.exp.res = future_lapply(expers.to.run, FUN=function(config.row) {
  source("TempVarSimFunctions.R")
  source("TempVarPlotFunctions.R")
  sim.specs = do.call('rbind', replicate(n=n.exper, expr=sim.configs[config.row,], simplify=F))[,tamp:=trange.seq]
  exp.res = RunExperiment(sim.specs, return.tseries = T)
  #exp.res[,trange:=trange.seq]
  plots = PlotExperimentPair(results = exp.res[c(1,2)],
                             configs = sim.specs,
                             pred.color = "black",
                             prey.color = "red", 
                             tmax = 2000)
  return(plots)
  #return(exp.res)
})
future:::ClusterRegistry("stop")

png("../results/popplot_stochevents.png", width=800, height=800)
grid.arrange(all.exp.res[[1]]$pop.plot,
             all.exp.res[[2]]$pop.plot)
dev.off()

# redo stochastic events, but for different transition probabilities and effect of event.
expers.to.run = 1
prob.mat = rbind(c(0.93, 0.07), c(0, 1)) # if a shock happens, it will persist
event.shift = 2
alt.exp.res = lapply(expers.to.run, FUN=function(config.row) {
  sim.specs = do.call('rbind', replicate(n=n.exper, expr=sim.configs[config.row,], simplify=F))[,tamp:=trange.seq]
  exp.res = RunExperiment(sim.specs, return.tseries = T)
  #exp.res[,trange:=trange.seq]
  plots = PlotExperimentPair(results = exp.res[c(1,2)],
                             configs = sim.specs,
                             pred.color = "black",
                             prey.color = "red", 
                             tmax = 2000)
  return(plots)
  #return(exp.res)
})

png("../results/popplot_stochevents_highac.png", width=800, height=400)
alt.exp.res[[1]]$pop.plot
dev.off()

