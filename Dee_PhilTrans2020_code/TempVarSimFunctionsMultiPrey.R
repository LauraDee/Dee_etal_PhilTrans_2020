library(deSolve)
CtoK = 273.15
kBoltz = 8.62E-5

# if functions are actually text (e.g., read in from config file) 
# convert to functions

ParseFunctionsMultiPrey = function(setup) {
  for(fn in c("PredConsumptionAtTemp")) {
    if(class(setup[[fn]]) == "character") {
      if(exists(setup[[fn]])) {
        setup[[fn]] = eval(parse(text=setup[[fn]]))
      }
    }
  }
  
  # eventually switch functions over to this model, where we don't need to check for existence
  # because the config file just contains a call to a generator as a column entry
  for(fn in c("PredHandleAtTemp", "PredMortAtTemp", "PredEffAtTemp", "PredAttackAtTemp")) {
    if(class(setup[[fn]]) == "character") {
      setup[[fn]] = eval(parse(text=setup[[fn]]))
    }
  }

  # For now we simply have the thermal responses of prey growth differ across potential prey.
  # These will be in config files in the same column buy separated by semicolons. We'll store
  # the resulting functions as a list.
  for(fn in c("PreyCapAtTemp", "PreyGrowthAtTemp")) {
    if(class(setup[[fn]]) == "character") {
      prey.fns = strsplit(setup[[fn]], split=";", fixed=T)[[1]]
      setup[[fn]] = lapply(prey.fns, FUN=function(prey.fn) {
        eval(parse(text=prey.fn))
      })
    }
  }
  
    
  return(setup)
}



# Allow for multiple prey populations. For now, assume attack rate is shared,
# handling time is identical, so that differences in functional response are driven
# entirely by differences in the prey populations, which are in turn affected by how
# temperature modulates prey growth across species.
# The prey.pops argument should be a vector.
TypeIIResponseMultiPrey = function(temp, pred.pops, prey.pop, params) {
  with(as.list(params), {
    pred.attack = pred.attack.at.topt*PredAttackAtTemp(temp)
    pred.handle = pred.handle.at.topt*PredHandleAtTemp(temp)
    return(pred.pop * pred.attack*prey.pops/(sum(pred.attack*pred.handle*prey.pops) + 1))
  })
}

# Function to describe continuous dynamics
SystemDynamicsMultiPrey <- function(t, state, params) {
  with(as.list(c(state,params)), {
    # get prey populations as a single vector to simplify later calcs
    prey.pop = state[grep("prey.pop", names(state))]
    # get temp based on time
    temp = TempAtTime(t)
    
    # Get parameters affected by temp.
    # Prey:
    # - growth
    prey.growth = prey.growth.at.topt*sapply(PreyGrowthAtTemp, FUN=function(prey.growth.fn) { 
      prey.growth.fn(temp)
    })
    # 
    # - carrying capacity (for now, fixed)
    prey.cap = prey.cap.at.topt * sapply(PreyCapAtTemp, FUN=function(prey.cap.fn) { 
      prey.cap.fn(temp)
    })
    
    # Predator:
    # - consumption - total, not per capita
    pred.cons = PredConsumptionAtTemp(temp, pred.pop, prey.pop, params)
    
    # - conversion efficiency (for now, fixed)
    pred.eff = pred.eff.at.topt * PredEffAtTemp(temp)
    
    # - mortality (for now, fixed)
    pred.mort = pred.mort.at.topt * PredMortAtTemp(temp)
    
    
    # Implement differential equations using parts constructed above
    dprey = prey.growth * prey.pop * (1 - prey.pop/prey.cap) - pred.cons - prey.harvest.f*prey.pop# prey dynamics
    dpred = sum(pred.eff*pred.cons) - pred.pop*(pred.mort + pred.harvest.f)
    return(list(c(dprey, dpred)))
  })
}

# Solve the system over a defined time horizon
# given initial conditions, dynamics, and parameters
SimulateModelMultiPrey = function(sim.config) {
  sim.config = ParseFunctionsMultiPrey(sim.config)
  n.prey = length(sim.config$prey.pop.init)
  initial.values = c(sim.config$prey.pop.init, sim.config$pred.pop.init)
  names(initial.values) = c(paste0("prey.pop", 1:n.prey), "pred.pop")

  times = seq(from=sim.config$time.start, to=sim.config$time.end, by=sim.config$time.step)

  population.paths = ode(y = initial.values, 
                         times = times, 
                         func = SystemDynamicsMultiPrey, 
                         parms = sim.config)
  return(population.paths)
}



GetCVMultiPrey = function(sim.results, start.frac=0.9, end.frac=1) {
  prey.colnums = grep("prey", colnames(sim.results))
  sim.length=nrow(sim.results)
  cv.start.index = floor(sim.length*(start.frac))
  cv.end.index = floor(sim.length*(end.frac))
  pred.traj = sim.results[cv.start.index:cv.end.index, "pred.pop"]
  pred.cv = sd(pred.traj)/mean(pred.traj)
  prey.cvs = sapply(prey.colnums, FUN=function(cnum) {
    prey.traj = sim.results[cv.start.index:cv.end.index, cnum]  
    prey.cv = sd(prey.traj)/mean(prey.traj)
  })
  cvs = c(pred.cv, prey.cvs)
  names(cvs) = c("pred.cv", paste0("prey.cv", 1:length(prey.cvs)))
  return(cvs)
}


RunExperimentMultiPrey = function(sim.specs, return.tseries=F) {
  n.exper = nrow(sim.specs)
  exp.res = lapply(1:n.exper, function(i) {
    sim.spec = as.list(sim.specs[i,])
    sim.spec = ParseFunctionsMultiPrey(sim.spec)
    
    # Parse parameters for different prey:
    if(is.character(sim.spec$prey.growth.at.topt)) {
      sim.spec$prey.growth.at.topt = as.numeric(strsplit(sim.spec$prey.growth.at.topt, ";")[[1]])
    }
    if(is.character(sim.spec$prey.cap.at.topt)) {
      sim.spec$prey.cap.at.topt = as.numeric(strsplit(sim.spec$prey.cap.at.topt, ";")[[1]])
    }
    if(is.character(sim.spec$prey.pop.init)) {
      sim.spec$prey.pop.init = as.numeric(strsplit(sim.spec$prey.pop.init, ";")[[1]])
    }
    
    # Construct temp sampling function
    sim.spec$TempAtTime = MakeSampleTempFun(type=sim.spec$TempAtTime,
                                            mean.temp=sim.spec$mean.temp,
                                            temp.amp=sim.spec$tamp,
                                            warm.rate=sim.spec$tamp/sim.spec$time.end,
                                            t.max=sim.spec$time.end)
    
    

    
    
    sim.res = SimulateModelMultiPrey(as.list(sim.spec))
    # average pops over last 10% of sample
    start = floor(0.9*nrow(sim.res))
    prey.colnums = grep("prey", colnames(sim.res))
    
    prey.end.mean = sapply(prey.colnums, FUN=function(cnum) {
      mean(sim.res[start:nrow(sim.res),cnum])
    })
    pred.end.mean = mean(sim.res[start:nrow(sim.res),"pred.pop"])
    
    # Because of cylcing and possibility of picking up partial cycles,
    # as an alternate summary of the long-run pop, we'll also look at
    # the mean of the extrema in that same final portion
    prey.end.ext.mean = sapply(prey.colnums, FUN=function(cnum) {
      mean(range(sim.res[start:nrow(sim.res),cnum]))
    })
    pred.end.ext.mean = mean(range(sim.res[start:nrow(sim.res),"pred.pop"]))
    
    
    prey.range = sapply(prey.colnums, FUN=function(cnum) {
      range(sim.res[,cnum])
    })
    
    pred.range = range(sim.res[,"pred.pop"])
    
    end.cv = GetCVMultiPrey(sim.res, start.frac = 0.9, end.frac=1)
    prey.end.cv = end.cv[prey.colnums]
    pred.end.cv = end.cv[1]
    start.cv = GetCVMultiPrey(sim.res, start.frac = 0, end.frac=0.1)
    prey.start.cv = start.cv[prey.colnums]
    pred.start.cv = start.cv[1]
    
    
    # Also compute and return average rates over the course of the 
    # simulation
    if(class(sim.spec$TempAtTime) == "character") {
      sim.spec$TempAtTime = eval(parse(text=sim.spec$TempAtTime))
    }
    temp.seq = sim.spec$TempAtTime(seq(from=sim.spec$time.start, to=sim.spec$time.end, by=sim.spec$time.step))
    attack.seq = sim.spec$PredAttackAtTemp(temp.seq)
    pred.attack.mean = mean(attack.seq)
    growth.seq = sapply(sim.spec$PreyGrowthAtTemp, FUN=function(prey.growth.fun) {
      prey.growth.fun(temp.seq)
    })
    prey.growth.mean = apply(growth.seq, MARGIN=2, FUN=mean)
    res = list(pred.end.mean=pred.end.mean,
               prey.end.mean=prey.end.mean,
               pred.end.ext.mean=pred.end.ext.mean,
               prey.end.ext.mean=prey.end.ext.mean,
               pred.range.min = pred.range[1],
               pred.range.max = pred.range[2],
               prey.range.min = prey.range[1,],
               prey.range.max = prey.range[2,],
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

