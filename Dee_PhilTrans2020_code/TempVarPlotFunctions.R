library(ggplot2)

nolegend = theme(legend.position = "none")
rightmargin = theme(plot.margin=margin(0,10,0,0))

LabelPanel = function(panel, l) {
  return(panel + nolegend + rightmargin + ggtitle(l) + theme(plot.title=element_text(vjust=0)))
}



PlotExperimentPair = function(results, configs, 
                              pred.color="green", prey.color="black", 
                              y.limits=c(-1, 21),
                              tmax=600) {
  
  pop.plot.data = data.table(t=results[[1]]$sim.res[,1],
                             prey.pop.1 = results[[1]]$sim.res[,2],
                             prey.pop.2 = results[[2]]$sim.res[,2],
                             pred.pop.1 = results[[1]]$sim.res[,3],
                             pred.pop.2 = results[[2]]$sim.res[,3])
  pop.plot = ggplot(pop.plot.data[t<=tmax,], aes(x=t)) + 
                    geom_line(aes(y=prey.pop.1, colour="Prey", linetype="Low")) + 
                    geom_line(aes(y=prey.pop.2, colour="Prey", linetype="High")) + 
                    geom_line(aes(y=pred.pop.1, colour="Predator", linetype="Low")) + 
                    geom_line(aes(y=pred.pop.2, colour="Predator", linetype="High")) + 
                    scale_colour_manual(name="Species",
                                        breaks=c("Prey", "Predator"),
                                        values=c(prey.color, pred.color)) +
                    scale_linetype_manual(name="Temp. Variability", 
                                          breaks=c("Low", "High"),
                                          values=c("solid", "dashed")) + 
                    ylab("Population") +
                    xlab("Time") +
                    scale_x_continuous(expand=c(0,0)) +
                    scale_y_continuous(expand=c(0,0), limits=y.limits) +
                    theme_bw()
  
  rate.plot.data = data.table(t=results[[1]]$sim.res[,1],
                              attack.rate.1 = results[[1]]$attack.seq,
                              attack.rate.2 = results[[2]]$attack.seq,
                              growth.rate.1 = results[[1]]$growth.seq,
                              growth.rate.2 = results[[2]]$growth.seq)
  
  rate.plot =  ggplot(rate.plot.data[t<=1,], aes(x=t)) + 
                      geom_line(aes(y=growth.rate.1, colour="Prey", linetype="Low")) + 
                      geom_line(aes(y=growth.rate.2, colour="Prey", linetype="High")) + 
                      geom_line(aes(y=attack.rate.1, colour="Predator", linetype="Low")) + 
                      geom_line(aes(y=attack.rate.2, colour="Predator", linetype="High")) + 
                      scale_colour_manual(name="Species",
                                          breaks=c("Prey", "Predator"),
                                          values=c(prey.color, pred.color)) +
                      scale_linetype_manual(name="Temp. Variability", 
                                            breaks=c("Low", "High"),
                                            values=c("solid", "dashed")) + 
                      ylab("Rate") +
                      xlab("Time") +
                      scale_x_continuous(expand=c(0,0)) +
                      scale_y_continuous(expand=c(0,0), limits=c(0,1.1)) +
                      theme_bw()
  
  return(list(pop.plot=pop.plot,
              rate.plot=rate.plot))
}






PlotResults = function(results, setup, pred.color="black", prey.color="grey", min.temp=0, max.temp=25) {
  
  setup = ParseFunctions(setup)
  pop.max = max(results[,2:3])
  results.dt = melt(data.table(results), id.vars='time')
  pop.plot = ggplot(results.dt, 
                    aes(x=time, 
                        y=value, 
                        group=variable, 
                        colour=variable)) + 
    geom_line() + 
    xlab("Period") + 
    ylab("Population (kg)") + 
    ylim(0,pop.max*1.05) +
    scale_colour_manual(values=c(prey.color, pred.color), 
                        labels=c("Prey", "Predator"),
                        name="Species") + 
    theme_minimal() + 
    theme(legend.position = "bottom")
  
  # Visualize the TPCs and how they relate to the temperatures
  PredTPC = function(x) { setup$PredAttackAtTemp(x) }
  PreyTPC = function(x) { setup$PreyGrowthAtTemp(x) }
  tseq = seq(from=0, to=max(results[,"time"]), by=0.01)
  temps = data.frame(temp=setup$TempAtTime(tseq))
  tpc.plot = ggplot(data.frame(temp=seq(from=min.temp, to=max.temp, by=0.001)), aes(x=temp)) + 
    stat_function(fun=PredTPC, mapping=aes(colour="Predator")) + 
    stat_function(fun=PreyTPC, mapping=aes(colour="Prey")) + 
    stat_density(data=temps, geom="line", mapping=aes(colour="Temperature")) +
    scale_colour_manual("", values = c(pred.color, prey.color, "red")) + 
    xlab("Temperature") + 
    ylab("TPC, Temperature Density") + 
    ylim(0, 1) +
    theme_minimal() + 
    theme(legend.position = "bottom")
  
  # Visualize rates (pred attack, prey growth) through time and averages.
  # Given speed of oscillation, only do this over a shorter time interval
  tseq = seq(from=0, to=10, by=0.01)
  temps = data.frame(temp=setup$TempAtTime(tseq))
  rate.data = data.frame(t = tseq,
                         temp = temps$temp, 
                         a = PredTPC(temps$temp),
                         r = PreyTPC(temps$temp))
  rate.plot = ggplot(rate.data, aes(x=t)) + 
    geom_line(aes(y=a, colour="Attack", linetype="Value at time")) + 
    geom_line(aes(y=r, colour="Prey Growth", linetype="Value at time")) + 
    geom_hline(aes(yintercept = mean(a), colour="Attack", linetype="Average")) + 
    geom_hline(aes(yintercept = mean(r), colour="Prey Growth", linetype="Average")) + 
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_colour_manual(name="", 
                        breaks = c("Attack", "Prey Growth"),
                        values=c("black", "grey")) + 
    scale_linetype_manual(name="", 
                          breaks = c("Average", "Value at time"),
                          values=c("dashed", "solid")) + 
    ylab("Parameter Value (relative to max)") + 
    xlab("Time") + 
    theme_minimal()
  
  return(list(pop.plot=pop.plot,
              tpc.plot=tpc.plot,
              rate.plot = rate.plot))
}

PlotExperiment = function(setup) {
  kTempAmp <<- setup$tamp
  kWarmRate <<- setup$tamp/setup$time.end
  topt.shift.pred = setup$topt.shift.pred
  topt.shift.prey = setup$topt.shift.prey
  
  setup$PredAttackAtTemp=MakeTPC(t.opt = kMeanTemp+topt.shift.pred, t.max = kMeanTemp+topt.shift.pred+2, t.sigma = 6)
  setup$PreyGrowthAtTemp=MakeTPC(t.opt = kMeanTemp+topt.shift.prey, t.max = kMeanTemp+topt.shift.prey+2, t.sigma = 6)
  
  setup = ParseFunctions(setup)
  
  # Plot 1: visualize temperature variability through time,
  # putting TPCs on y axis as secondary plots
  tseq = seq(from=0, to=4, by=0.01)
  temps = data.frame(temp=setup$TempAtTime(tseq))
  rate.data = data.frame(t = tseq,
                         temp = temps$temp, 
                         a = setup$PredAttackAtTemp(temps$temp),
                         r = setup$PreyGrowthAtTemp(temps$temp))
  
  
  setup.plot = ggplot(rate.data, aes(x=t)) + 
    geom_line(aes(y=temp), colour="red") +
    geom_hline(aes(yintercept=min(temp)), linetype="dashed", colour="red") +
    geom_hline(aes(yintercept=max(temp)), linetype="dashed", colour="red") +
    geom_segment(aes(y = min(temp),
                     yend = min(temp)+2*setup$topt.shift.pred,
                     x = 3.5,
                     xend = 3.5),
                 arrow=arrow(length = unit(0.03, "npc")),
                 colour="red") +
    geom_segment(aes(y = max(temp),
                     yend = max(temp)+2*setup$topt.shift.prey,
                     x = 3.5,
                     xend = 3.5),
                 arrow=arrow(length = unit(0.03, "npc")),
                 colour="red") +
    annotate(geom = "point", x = 3.5, y=23, size=8, colour="red") +
    annotate(geom = "text", x = 3.5, y=23, label="1", size=6, colour="white") +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), limits=c(0, 30)) +
    scale_colour_manual(name="", 
                        breaks = c("Attack", "Prey Growth"),
                        values=c("black", "grey")) + 
    scale_linetype_manual(name="", 
                          breaks = c("Average", "Value at time"),
                          values=c("dashed", "solid")) + 
    ylab("Temperature") + 
    xlab("Time") + 
    theme_classic()
  
  tpc.ymax = 1.3
  tpc.arrow.y = 1+(tpc.ymax-1)*0.3
  tpc.plot = ggplot(rate.data, aes(x=temp)) + 
    stat_function(fun=setup$PredAttackAtTemp, colour="black") +
    stat_function(fun=setup$PreyGrowthAtTemp, colour="gray") + 
    #geom_vline(aes(xintercept=min(temp)), linetype="dashed", colour="red") +
    #geom_vline(aes(xintercept=max(temp)), linetype="dashed", colour="red") +
    geom_segment(aes(x = kMeanTemp+setup$topt.shift.pred,
                     xend = kMeanTemp+setup$topt.shift.pred,
                     y = 1,
                     yend=tpc.ymax),
                 linetype="dashed",
                 colour="black") +
    geom_segment(aes(x = kMeanTemp+setup$topt.shift.pred,
                     xend = kMeanTemp+3*setup$topt.shift.pred,
                     y = tpc.arrow.y,
                     yend=tpc.arrow.y),
                 arrow=arrow(length = unit(0.03, "npc")),
                 colour="black") +
    geom_segment(aes(x = kMeanTemp+setup$topt.shift.prey,
                     xend = kMeanTemp+setup$topt.shift.prey,
                     y = 1,
                     yend=tpc.ymax),
                 linetype="dashed",
                 colour="gray") +
    geom_segment(aes(x = kMeanTemp+setup$topt.shift.prey,
                     xend = kMeanTemp+3*setup$topt.shift.prey,
                     y = tpc.arrow.y,
                     yend=tpc.arrow.y),
                 arrow=arrow(length = unit(0.03, "npc")),
                 colour="gray") +
    annotate(geom = "point", x=kMeanTemp+4*setup$topt.shift.prey, y=tpc.arrow.y, size=8, colour="black") +
    annotate(geom = "text", x=kMeanTemp+4*setup$topt.shift.prey, y=tpc.arrow.y, label="2", size=6, colour="white") +
    scale_x_continuous(expand=c(0,0), limits=c(0,30)) + 
    scale_y_continuous(expand=c(0,0), limits=c(0, tpc.ymax)) + 
    ylab("Rate scaling") +
    coord_flip() +
    theme_classic() + 
    theme(axis.text.y = element_text(colour="white"),
          axis.title.y=element_text(colour="white"))

  #tpc.plot
  return(grid.arrange(setup.plot,
               tpc.plot,
               nrow=1))
}


MakeTPCPlot = function(sim.config, prey.color="black", pred.color="red") {
  topt.pred = sim.config[,mean.temp+topt.shift.pred]
  td.pred = topt.pred + 1
  ed.pred = sim.config[,pred.ed]
  tpc.pred =MakeTPCAG(t.ref = topt.pred,
                      t.d = topt.pred+1,
                      e.a = ed.pred/(1+exp(ed.pred*(td.pred/topt.pred-1)/(kBoltz*td.pred))), 
                      e.d = ed.pred)
  
  topt.prey = sim.config[,mean.temp+topt.shift.prey]
  td.prey = topt.prey + 1
  ed.prey = sim.config[,prey.ed]
  tpc.prey =MakeTPCAG(t.ref = topt.prey,
                      t.d = topt.prey+1,
                      e.a = ed.prey/(1+exp(ed.prey*(td.prey/topt.prey-1)/(kBoltz*td.prey))), 
                      e.d = ed.prey)
  temps = seq(from=273, to=300, by=0.1)
  temp.dt = data.table(temps=temps,
                       pred.tpc = tpc.pred(temps),
                       prey.tpc = tpc.prey(temps))
  tpc.plot = ggplot(temp.dt, aes(x=temps-273.15)) + 
    geom_line(aes(y=pred.tpc, colour="Predator"), size=2) + 
    geom_line(aes(y=prey.tpc, colour="Prey"), size=1.2) + 
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), limits=c(0, 1.1)) +
    scale_colour_manual(name="Species",
                        breaks=c("Predator", "Prey"),
                        values=c(pred.color, prey.color)) +
    xlab("Temperature (C)") + 
    ylab("Parameter (relative to max)") + 
    theme_bw() 
  return(tpc.plot)
}



MakeTempScenarioPlot = function(plot.dt, sp) {
  sp.plot = ggplot(plot.dt, aes(x=time)) + 
    geom_line(aes(y=get(paste0("var.",sp, ".pop")), colour="Var")) + 
    geom_line(aes(y=get(paste0("warm.",sp, ".pop")), colour="Warm")) + 
    geom_line(aes(y=get(paste0("warm.var.",sp, ".pop")), colour="WarmVar")) +
    geom_line(aes(y=get(paste0("warm.inc.var.",sp, ".pop")), colour="WarmIncVar")) + 
    xlab("Time") + 
    ylab("Population") + 
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_colour_manual(name="Temperature Scenario",
                        breaks=c("Var", "Warm", "WarmVar", "WarmIncVar"),
                        values=c("Var"="blue",
                                 "Warm"="red",
                                 "WarmVar" = "purple",
                                 "WarmIncVar" = "darkorange1"),
                        labels=c("Fixed Variability, No Warming", "Warming + No Variability", "Warming + Fixed Variability", "Warming + Increasing Variability")) +
    theme_bw() + 
    theme(legend.position="bottom")
  return(sp.plot)
}

MakeHeatmapPlot = function(heatmap.res, col.lim=c(0,5)) {
  heatmap.res[,min.cv:=pmin(pred.end.cv, prey.end.cv)]
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
    xlab("Prey TPC offset from mean temperature (C)") +
    ylab("Temperature amplitude (C)") +
    theme_bw() + 
    theme(legend.position="bottom") 
}


MakeThreePanelHeatmapPlot = function(fname1, fname2, fname3, out.file,
                                     fig.width=900,
                                     fig.height=300) {
  
  heatmap.res.1 = fread(fname1)
  heatmap.res.2 = fread(fname2)
  heatmap.res.3 = fread(fname3)
  pred.range = range(c(heatmap.res.1$pred.end.mean, 
                       heatmap.res.2$pred.end.mean, 
                       heatmap.res.3$pred.end.mean))
  col.range = c(0, pred.range[2])
  
  
  heatmap.1 = MakeHeatmapPlot(heatmap.res.1, col.lim = col.range) +
    geom_text(x=-0.5, y=1.5, label="Unstable Coexistence", colour="white") +
    geom_text(x=1.0, y=4.0, label="Stable Coexistence", colour="white") +
    geom_text(x=2.3, y=1.2, label="Collapse", colour="white", angle=90)
  
  heatmap.2 = MakeHeatmapPlot(heatmap.res.2, col.lim = col.range) +
    geom_text(x=-0.5, y=1.2, label="Stable Coexistence", colour="white") +
    geom_text(x=1.5, y=2.5, label="Collapse", colour="white")
  
  heatmap.3 = MakeHeatmapPlot(heatmap.res.3, col.lim = col.range) + 
    geom_text(x=-0.5, y=2.5, label="Unstable Coexistence", colour="white") +
    geom_text(x=1.8, y=2.5, label="Stable Coexistence", angle=75, colour="white") +
    geom_text(x=2.3, y=0.7, label="Collapse", angle=90, colour="white")
  
  
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
                             heights=c(0.1, 0.7, 0.2)
  )
  pdf(file=out.file, width = fig.width, height=fig.height)
  grid.draw(pspace.plots)
  dev.off()
}

