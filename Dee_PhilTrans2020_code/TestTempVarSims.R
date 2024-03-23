# Test the simulation code for a deterministic case with known equilibria
library(data.table) #for loading stocahstic temps (irrelevant, but gives something to use)

# Load simulation functions
source("TempVarSimFunctions.R")
source("TempVarPlotFunctions.R")

# Set up TPCs that do not vary with temp
PreyTPC = ConstTPC
PredTPC = ConstTPC

# Set up temperature sampling based on historical monthly temperatures
lme.month.data = fread("../data/derived/sstByLMEandMonth.csv")
setorder(lme.month.data, year, month)
lme.month.temps = lme.month.data[LME==7, temp-273.15] # pick out NE US CS



# linearly interpolate between observed values
SampleTempNEUSCS = function(t) {
  index = floor(t)+1
  return(lme.month.temps[index] + (t-index+1)*(lme.month.temps[index+1]-lme.month.temps[index]))
}

testsim.constants = c(
  # Pred params
  
  # Attack rate of 654.9 from Moustahfid et al. supplement Table 2.
  # Note those params were estimated based on prey densities measured
  # in thousands of fish per km^2.
  # can load change in units onto attack rate. 
  # /1000 for change to actual numbers of fish/nautical mile sq
  # /0.2 to change to kg/nmisq instead of individuals/nmisq
  # /35000 to account for density vs raw numbers 
  # (approx size of Scotian Shelf; our carrying cap # comes from canadian herring pop
  pred.attack.max = 654.9/1000/0.2/35000, 
  pred.ind.size = 2000,         # predator individual size in g
  pred.handle = 0.001,          # from Moustahfid et al. supplement Table 2
  pred.mort = 0.2,              # Can we find a better number?
  # 0.2 commonly used in management.
  # See https://pdfs.semanticscholar.org/d21c/d49086ccbd42a9185341a743a1ab197913d3.pdf -- experiment group H is ~
  # 8 months, suggesting (3.1+22.6+3.1+15.6+25.8+29.0+87.9)/7*12/8= approximately 0.4 annual M.
  pred.eff = 0.42,   # 0.31 from control group in Lemieux et al. Table 1: (https://link.springer.com/article/10.1023/A:1007791019523)
  # Also broadly in line with peak efficiency in Fig 3 of (https://www.nrcresearchpress.com/doi/full/10.1139/F08-159#.XadcJ5NKi35)
  
  # Prey params
  # Capacity and growth pulled from grid search estimates for Canadian Herring
  # in https://www.sciencedirect.com/science/article/pii/S0304380012005236.
  # Specifically, text just after Table 1: "For Canadian herring the Grid Search estimates are obtained as r = 0.52, a = 800, θ = 0.0006, K = 1210"
  # change units for capacity to be kg: First *1000 is for thousands of tonnes, second is for kg instead of tonnes)
  prey.cap = 1210*1000*1000,    
  prey.max.growth = 0.52, 
  
  # Select functions of interest
  
  # Pred/prey interaction
  PredConsumptionAtTemp = function(temp, pred.pop, prey.pop, params) {
    # Test Type II response given in weight of prey consumed per predator,
    # converted to consumption per g of predator.
    with(params, {
      return(TypeIIResponse(temp, pred.pop, prey.pop, params)/pred.ind.size)
    })
  },
  PredAttackAtTemp = PredTPC,
  PreyGrowthAtTemp = PreyTPC,
  #TempAtTime = SampleTempCos,
  TempAtTime = SampleTempNEUSCS,
  
  # Harvest rates
  pred.harvest.f = 0.00,
  prey.harvest.f = 0.00
)


# Initialize to coexistence equilibrium, ensure it remains
# there
prey.pop.init = with(testsim.constants, {
  pred.mort*pred.ind.size/(pred.attack.max*(pred.eff-pred.mort*pred.handle*pred.ind.size))
})
pred.pop.init = with(testsim.constants, {
  pred.ind.size*prey.max.growth*(1-prey.pop.init/prey.cap)*(pred.attack.max*pred.handle*prey.pop.init+1)/pred.attack.max
})

test.initial.values = c(prey.pop=prey.pop.init,  # Initial prey biomass (kg)
                        pred.pop=pred.pop.init)   # Initial predator biomass (kg)

# Test system dynamics in isolation
SystemDynamics(0, test.initial.values, testsim.constants)


# Test full simulation
times = seq(0,length(lme.month.temps)-1,by=0.01)
times = times[-length(times)]


# Run simulations
test.results = SimulateModel(testsim.constants, test.initial.values, times)
test.plots = PlotResults(results=test.results,
                         setup=testsim.constants,
                         max.temp = 30)

test.plots$pop.plot
