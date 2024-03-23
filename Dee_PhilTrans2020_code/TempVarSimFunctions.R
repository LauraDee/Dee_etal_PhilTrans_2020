library(deSolve)
CtoK = 273.15
kBoltz = 8.62E-5

# Stochastic events (a la ENSO)
# For now, probabilities fixed through time
prob.mat = rbind(c(0.93, 0.07), # prob if not currently El Niño
                 c(0.5, 0.5))   # prob if currently El Niño
event.shift = 0.5 # temperature shift (degrees K) if event is active

# if functions are actually text (e.g., read in from config file) 
# convert to functions

ParseFunctions = function(setup) {
  for(fn in c("PredConsumptionAtTemp")) {
    if(class(setup[[fn]]) == "character") {
      if(exists(setup[[fn]])) {
        setup[[fn]] = eval(parse(text=setup[[fn]]))
      }
    }
  }
  
  # eventually switch functions over to this model, where we don't need to check for existence
  # because the config file just contains a call to a generator as a column entry
  for(fn in c("PredHandleAtTemp", "PredMortAtTemp", "PredEffAtTemp", "PreyCapAtTemp","PredAttackAtTemp", "PreyGrowthAtTemp")) {
    if(class(setup[[fn]]) == "character") {
      setup[[fn]] = eval(parse(text=setup[[fn]]))
    }
  }
  
  return(setup)
}

##########
# Functions to return parameters which are 
# temperature dependent, i.e. TPCs
##########

# Basic Gaussian and inverted Gaussian
GaussianTPC = function(temp, t.opt, t.sigma, inverted=F) {
  exp.mult = 2*(inverted-0.5)
  return(exp(exp.mult*(temp-t.opt)^2/(2*t.sigma)))
}

# Basic exponential
ExponentialTPC = function(temp, t.ref, e.a) {
  return(exp(-e.a/kBoltz*(1/temp-1/t.ref)))
}

# From Deutsch et al. 2008 PNAS
AsymmTPC = function(temp, t.opt, t.max, t.sigma) {
  tpc = exp(-(temp-t.opt)^2/(2*t.sigma))
  tpc[temp>t.opt] = pmax(0, 1-((temp[temp>t.opt]-t.opt)/(t.opt-t.max))^2)
  return(tpc)
}


# From Gårdmark et al. in prep. Note it deals with temps in Kelvin, not C.
AsymmTPCAG = function(temp, t.ref, t.d, e.a, e.d) {
  exp(-e.a/kBoltz * (1/temp-1/t.ref)) / (1 + exp(-e.d/kBoltz*(1/temp-1/t.d))) * (1 + exp(-e.d/kBoltz*(1/t.ref-1/t.d)))
}

# quadratic TPC
SymmTPC = function(temp, t.opt, t.width) {
  return(pmax(1-t.width^2*(temp-t.opt)^2, 0))
}

# Fixed performance (constant TPC)
ConstTPC = function(temp) {
  return(rep(1, length(temp)))
}

# Generators

MakeTPCConst = function(value) {
  return(function(temp) {
    return(value)
  })
}

MakeTPCAsymm = function(t.opt, t.max, t.sigma) {
  return(function(temp) {
    return(AsymmTPC(temp, t.opt, t.max, t.sigma))
  })
}

MakeTPCAG = function(t.ref, t.d, e.a, e.d) {
  return(function(temp) {
    return(AsymmTPCAG(temp, t.ref, t.d, e.a, e.d))
  })
}

MakeTPCQuad = function(t.opt, t.max) {
  t.width = sqrt(1/(t.max-t.opt)^2)
  return(function(temp) {
    return(SymmTPC(temp, t.opt, t.width))
  })
}

MakeTPCAGAtMean = function(mean.temp, t.ref, t.d, e.a, e.d) {
  tpc.val = AsymmTPCAG(mean.temp, t.ref, t.d, e.a, e.d)
  return(function(temp) {
    return(tpc.val)
  })
}

MakeTPCGaussian = function(t.opt, t.sigma) {
  return(function(temp) {
    return(GaussianTPC(temp, t.opt, t.sigma, inverted=F))
  })
}

MakeTPCInvertedGaussian = function(t.opt, t.sigma) {
  return(function(temp) {
    return(GaussianTPC(temp, t.opt, t.sigma, inverted=T))
  })
}

MakeTPCExponential = function(t.ref, e.a) {
  return(function(temp) {
    return(ExponentialTPC(temp, t.ref, e.a))
  })
}


##################
# Functions to get temperature as a function of t.
# will provide primitives and generators for 
# functions that are only a function of t
##################

# Primitives
SampleTempCos = function(t, mean.temp, temp.amp) {
  return(mean.temp-temp.amp*cos(2*pi*t))
}

SampleTempSaw = function(t, mean.temp, temp.amp) {
  return(mean.temp + temp.amp*(0.5 - 2*abs((t-floor(t))-0.5)))
}

ConstantTemp = function(t, mean.temp, temp.amp) {
  return(mean.temp)
}

SampleTempWarming = function(t, mean.temp, temp.amp, warm.rate) {
  return(mean.temp - temp.amp/2 + t*warm.rate)
}

SampleTempWarmingAndSaw = function(t, mean.temp, temp.amp, warm.rate) {
  return(mean.temp - temp.amp/2 + t*warm.rate + temp.amp*(0.5 - 2*abs((t-floor(t))-0.5)))
}

SampleTempWarmingAndGrowingSaw = function(t, mean.temp, temp.amp, warm.rate) {
  return(mean.temp - temp.amp/2 + t*warm.rate + (temp.amp + t*warm.rate)*(0.5 - 2*abs((t-floor(t))-0.5)))
}



# DKO simulator function which takes
# x: time steps
# l_x: length-scale transfer function to modulate smoothness/autocorrelation as a function of (temporal) distance
# n: number of time series to draw
# phi: modulate length-scale overall
# var_offset: increment to variance (presumably to preserve/ensure positive-definiteness)
# sigma: vector (or scalar) of standard deviations
exp_ND_fun <- function(x, l_x, n = 1, phi = 0.5, var_offset = 0.0001, sigma = 1) {
  l_x_2 <- l_x / mean(l_x) / phi
  distance <- as.matrix(dist(x, upper = T, diag = T))
  C_mat <- matrix(NA, ncol = length(x), nrow = length(x))
  for (i in 1 : (length(x) - 1)) {
    for (j in (i + 1) : (length(x))) {
      C_mat[i, j] = ((l_x_2[i] ^ 2) ^ 0.25 * (l_x_2[j] ^ 2) ^ 0.25) * (0.5 * l_x_2[i] ^ 2 + 0.5 * l_x_2[j] ^ 2) ^ (- 0.5) * exp( - distance[i, j] ^ 2 / (0.5 * l_x_2[i] ^ 2 + 0.5 * l_x_2[j] ^ 2));
      C_mat[j, i] = C_mat[i, j];
    }
  }
  for (k in 1 : (length(x))) C_mat[k, k] = 1 + var_offset;
  rand <- matrix(rnorm(length(x) * n, 0, 1), ncol = n)
  temp <- t(chol(C_mat)) %*% (sigma * rand)
  return(temp)
}



# Generator

MakeSampleTempFun = function(type, mean.temp, temp.amp, warm.rate, t.max) {
  
  if(type=="Cos") {
    return(function(t) {
      SampleTempCos(t, mean.temp, temp.amp)
    }) 
  } else if(type=="Saw") {
    return(function(t) {
      SampleTempSaw(t, mean.temp, temp.amp)
    }) 
  } else if(type=="Warming") {
    return(function(t) {
      SampleTempWarming(t, mean.temp, temp.amp, warm.rate)
    }) 
  } else if(type=="WarmingAndSaw") {
    return(function(t) {
      SampleTempWarmingAndSaw(t, mean.temp, temp.amp, warm.rate)
    }) 
  } else if(type=="WarmingAndGrowingSaw") {
    return(function(t) {
      SampleTempWarmingAndGrowingSaw(t, mean.temp, temp.amp, warm.rate)
    }) 
  } else if(type=="Stoch") {
    # to allow for autocorrelation, we either need to track temperature as a state variable
    # or we can pre-sample a series w/autocorrelation, and generate a function that 
    # interpolates between those sampled points. The latter lets us keep the basic structure
    # of the simulations with only populations as state variables.
    # we'll do that using DKO's function which samples from a gaussian process with covariance 
    # structure that depends upon temporal distance (exp_ND_fun above)
    
    # If we want to keep the variance comparable to the base linear oscillation,
    # we could halve the amplitude of the oscillation and set sigma to that same value
    # for the added gaussian process. We'll do that for the start of the simulation (allowing sigma
    # to grow over time)
    
    #temp.deviations = filter(rnorm(t.max+1), filter=rep(1,3), method="convolution", circular=T)
    steps.per.year = 4
    time.steps = 1:(t.max*steps.per.year)
    # logistic transfer function: length-scale increasing through time so that we get more prolonged
    # "events" in the future
    scaled.tsteps <- time.steps / (length(time.steps) / 10) - 5
    l_x <- exp(scaled.tsteps)/(1 + exp(scaled.tsteps)) 
    
    temp.deviations = exp_ND_fun(time.steps, l_x, n = 1, phi = 1, sigma = (temp.amp + (1:t.max)*warm.rate)/16) 
    
    return(function(t) {
      #base.temp = SampleTempWarmingAndSaw(t, mean.temp, temp.amp/2, warm.rate) #start with some noiseless trend/seasonality
      base.temp = SampleTempSaw(t, mean.temp, temp.amp*15/16)
      return(base.temp + temp.deviations[floor(t*steps.per.year+1)])  # add autocorrelated noise to get at shocks, extreme events
    }) 
  } else if(type=="StochEvents") {
    # Introduce stochasticity/climate variability in the form of "events" a la El Niño.
    # We'll use a simple Markov process where a shock could come in the form of a hot year or years a la El Niño.
    # To facilitate comparisons, we'll keep the underlying seasonal oscillation but just change the range.
    time.steps = 1:t.max
    
    # Sample a stochastic series of events
    # approximate ENSO transitions between El Niño and not (ignoring La Niña) using
    # a simple markov process.
    regimes = rep(0, length(time.steps))
    
    for(t in time.steps[2:length(time.steps)]) {
      regimes[t] = sample(x=c(0,1), size=1, prob=prob.mat[regimes[t-1]+1,])  
    }

    return(function(t) { 
      cur.regime = regimes[floor(t)+1]
      mean.temp - temp.amp/2 + temp.amp*(0.5 - 2*abs((t-floor(t))-0.5)) + cur.regime*event.shift
    })
    
  } else if(type=="Constant") {
    return(function(t) {
      ConstantTemp(t, mean.temp, temp.amp)
    }) 
  }
}



#############
# Function(s) defining interaction between predator and prey
#############
# could pull pred.pop multiplier outside, but this allows for more general
# nonlinear role of predator population.
TypeIIResponse = function(temp, pred.pop, prey.pop, params) {
  with(as.list(params), {
    pred.attack = pred.attack.at.topt*PredAttackAtTemp(temp)
    pred.handle = pred.handle.at.topt*PredHandleAtTemp(temp)
    return(pred.pop * pred.attack*prey.pop/(pred.attack*pred.handle*prey.pop + 1))
  })
}




NoInteractionFn = function(pred.pop, prey.pop, params) {
  return(0)
}


# Function to describe continuous dynamics
SystemDynamics <- function(t, state, params) {
  with(as.list(c(state,params)), {
    # get temp based on time
    temp = TempAtTime(t)

    # Get parameters affected by temp.
    # Prey:
    # - growth
    prey.growth = prey.growth.at.topt*PreyGrowthAtTemp(temp)
    # 
    # - carrying capacity (for now, fixed)
    prey.cap = prey.cap.at.topt * PreyCapAtTemp(temp)   

    # Predator:
    # - consumption - total, not per capita
    pred.cons = PredConsumptionAtTemp(temp, pred.pop, prey.pop, params)
    
    # - conversion efficiency (for now, fixed)
    pred.eff = pred.eff.at.topt * PredEffAtTemp(temp)
    
    # - mortality (for now, fixed)
    pred.mort = pred.mort.at.topt * PredMortAtTemp(temp)
    
    
    # Implement differential equations using parts constructed above
    dprey = prey.growth * prey.pop * (1 - prey.pop/prey.cap) - pred.cons - prey.harvest.f*prey.pop# prey dynamics
    dpred = pred.eff*pred.cons - pred.pop*(pred.mort + pred.harvest.f)
    return(list(c(dprey, dpred)))
  })
}

# Solve the system over a defined time horizon
# given initial conditions, dynamics, and parameters
SimulateModel = function(sim.config) {
  initial.values = c(prey.pop=sim.config$prey.pop.init,  
                     pred.pop=sim.config$pred.pop.init)
  times = seq(from=sim.config$time.start, to=sim.config$time.end, by=sim.config$time.step)
  sim.config = ParseFunctions(sim.config)
  population.paths = ode(y = initial.values, 
                         times = times, 
                         func = SystemDynamics, 
                         parms = sim.config)
  return(population.paths)
}

GetCV = function(sim.results, start.frac=0.9, end.frac=1) {
  sim.length=nrow(sim.results)
  cv.start.index = floor(sim.length*(start.frac))
  cv.end.index = floor(sim.length*(end.frac))
  pred.traj = sim.results[cv.start.index:cv.end.index, "pred.pop"]
  prey.traj = sim.results[cv.start.index:cv.end.index, "prey.pop"]
  return(list(prey.cv = sd(prey.traj)/mean(prey.traj),
              pred.cv = sd(pred.traj)/mean(pred.traj)))
  
}


RunExperiment = function(sim.specs, return.tseries=F) {
  n.exper = nrow(sim.specs)
  exp.res = lapply(1:n.exper, function(i) {
    sim.spec = as.list(sim.specs[i,])
    sim.spec = ParseFunctions(sim.spec)   
    
    # Construct temp sampling function
    sim.spec$TempAtTime = MakeSampleTempFun(type=sim.spec$TempAtTime,
                                            mean.temp=sim.spec$mean.temp,
                                            temp.amp=sim.spec$tamp,
                                            warm.rate=sim.spec$tamp/sim.spec$time.end,
                                            t.max=sim.spec$time.end)
    
    
    # Initialize populations to fixed rate coexistence equilibrium if not given
    if(is.na(sim.spec$prey.pop.init)) {
      sim.spec$prey.pop.init = sim.spec$pred.mort.at.topt/(sim.spec$pred.attack.at.topt*(sim.spec$pred.eff.at.topt-sim.spec$pred.mort.at.topt*sim.spec$pred.handle.at.topt))
    }
    if(is.na(sim.spec$pred.pop.init)) {
      sim.spec$pred.pop.init = sim.spec$prey.growth.at.topt/sim.spec$pred.attack.at.topt * (1-sim.spec$prey.pop.init/sim.spec$prey.cap.at.topt) * (sim.spec$pred.attack.at.topt*sim.spec$pred.handle.at.topt*sim.spec$prey.pop.init + 1)
    }
    
    
    sim.res = SimulateModel(as.list(sim.spec))
    # average pops over last 10% of sample
    start = floor(0.9*nrow(sim.res))
    prey.end.mean = mean(sim.res[start:nrow(sim.res),"prey.pop"])
    pred.end.mean = mean(sim.res[start:nrow(sim.res),"pred.pop"])
    
    # Because of cylcing and possibility of picking up partial cycles,
    # as an alternate summary of the long-run pop, we'll also look at
    # the mean of the extrema in that same final portion
    prey.end.ext.mean = mean(range(sim.res[start:nrow(sim.res),"prey.pop"]))
    pred.end.ext.mean = mean(range(sim.res[start:nrow(sim.res),"pred.pop"]))
    
    
    prey.range = range(sim.res[,"prey.pop"])
    pred.range = range(sim.res[,"pred.pop"])
    end.cv = GetCV(sim.res, start.frac = 0.9, end.frac=1)
    prey.end.cv = end.cv$prey.cv
    pred.end.cv = end.cv$pred.cv
    start.cv = GetCV(sim.res, start.frac = 0, end.frac=0.1)
    prey.start.cv = start.cv$prey.cv
    pred.start.cv = start.cv$pred.cv
    
    
    # Also compute and return average rates over the course of the 
    # simulation
    if(class(sim.spec$TempAtTime) == "character") {
      sim.spec$TempAtTime = eval(parse(text=sim.spec$TempAtTime))
    }
    temp.seq = sim.spec$TempAtTime(seq(from=sim.spec$time.start, to=sim.spec$time.end, by=sim.spec$time.step))
    attack.seq = sim.spec$PredAttackAtTemp(temp.seq)
    pred.attack.mean = mean(attack.seq)
    growth.seq = sim.spec$PreyGrowthAtTemp(temp.seq)
    prey.growth.mean = mean(growth.seq)
    res = list(pred.end.mean=pred.end.mean,
               prey.end.mean=prey.end.mean,
               pred.end.ext.mean=pred.end.ext.mean,
               prey.end.ext.mean=prey.end.ext.mean,
               pred.range.min = pred.range[1],
               pred.range.max = pred.range[2],
               prey.range.min = prey.range[1],
               prey.range.max = prey.range[2],
               pred.end.cv = pred.end.cv,
               prey.end.cv = prey.end.cv,
               pred.start.cv = pred.start.cv,
               prey.start.cv = prey.start.cv,
               pred.attack.mean = pred.attack.mean,
               prey.growth.mean = prey.growth.mean)
    if(return.tseries) {
      res = c(res, list(sim.res=sim.res, attack.seq=attack.seq, growth.seq=growth.seq))
    }
    return(res)
  })
  return(exp.res)
}




ParseTPCSpec = function(tpc.fn.text) {
  call.parts = strsplit(substr(tpc.fn.text, start=1, stop=nchar(tpc.fn.text)-1), split="\\(")[[1]]
  return(list(fn=call.parts[1],
              params=strsplit(call.parts[2], ",")[[1]]))
}

GenerateHeatmapConfigsForRow = function(sim.configs, row) {
  topt.shift.pred = seq(from=-2.5, to=2.5, by=0.1)
  tamp = seq(from=0, to=5, by=0.1)
  param.grid = data.table(expand.grid(topt.shift.pred=topt.shift.pred, tamp=tamp))
  param.grid[,topt.shift.prey:=-topt.shift.pred]
  
  prey.growth.tpc.info = ParseTPCSpec(sim.configs$PreyGrowthAtTemp[row])
  pred.attack.tpc.info = ParseTPCSpec(sim.configs$PredAttackAtTemp[row])
  
  heatmap.configs = sim.configs[rep(row, nrow(param.grid)),]
  heatmap.configs[,`:=`(topt.shift.pred=param.grid$topt.shift.pred,
                        topt.shift.prey=param.grid$topt.shift.prey,
                        tamp=param.grid$tamp)]
  # now actually update the functions for predator attack and prey growth TPCs
  # based on these shift parameters
  
  prey.topt = heatmap.configs[,mean.temp+topt.shift.prey]
  prey.td = prey.topt+1 # should infer this from base function call in config file
  prey.ed = 20 # should parse this from base function call in config file
  
  # these are values such that the t.ref provided is actually t.opt, permitting a cleaner
  # interpretation of the reference temperature.
  prey.ea = prey.ed/(1+exp(prey.ed*(prey.td/prey.topt-1)/(kBoltz*prey.td)))
  prey.mean.param = ifelse(prey.growth.tpc.info$fn=="MakeTPCAGAtMean", 
                           paste0(sim.configs$mean.temp[row], ","),
                           "")
  
  pred.topt = heatmap.configs[,mean.temp+topt.shift.pred]
  pred.td = pred.topt + 1 # should infer this from base function call in config file
  pred.ed = 20 # should parse this from base function call in config file
  
  # these are values such that the t.ref provided is actually t.opt, permitting a cleaner
  # interpretation of the reference temperature.
  pred.ea = pred.ed/(1+exp(pred.ed*(pred.td/pred.topt-1)/(kBoltz*pred.td)))
  pred.mean.param = ifelse(pred.attack.tpc.info$fn=="MakeTPCAGAtMean", 
                           paste0(sim.configs$mean.temp[row], ","),
                           "")
  
  heatmap.configs[,`:=`(PreyGrowthAtTemp=paste0(prey.growth.tpc.info$fn,"(", 
                                                prey.mean.param,
                                                prey.topt,",", 
                                                prey.td,",",
                                                prey.ea,",",
                                                prey.ed,
                                                ")"))]
  
  heatmap.configs[,`:=`(PredAttackAtTemp=paste0(pred.attack.tpc.info$fn,"(", 
                                                pred.mean.param,
                                                pred.topt,",", 
                                                pred.td,",",
                                                pred.ea,",",
                                                pred.ed,
                                                ")"))]
  return(heatmap.configs)
}

RunHeatmapSims = function(sim.configs, row, n.procs = 6) {
  
  heatmap.configs = GenerateHeatmapConfigsForRow(sim.configs, row)
  
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
  fwrite(x=all.heatmap.res, file=paste0("../results/heatmap_res",row,".csv"))
}
