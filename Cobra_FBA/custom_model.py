# Toy central carbon metabolism model for E. coli using COBRApy
# Objective: biomass drain reaction
#
# pip install cobra

from cobra import Model, Reaction, Metabolite

def add_exchange(model, met, rid, lb=-1000.0, ub=1000.0):
    r = Reaction(rid)
    r.lower_bound = lb
    r.upper_bound = ub
    r.add_metabolites({met: -1.0})
    model.add_reactions([r])
    return r

def add_sink(model, met, rid, lb=0.0, ub=1000.0):
    r = Reaction(rid)
    r.lower_bound = lb
    r.upper_bound = ub
    r.add_metabolites({met: -1.0})
    model.add_reactions([r])
    return r

model = Model("ecoli_central_carbon_toy")

c, e = "c", "e"

# Metabolites
glc__D_e = Metabolite("glc__D_e", compartment=e)
glc__D_c = Metabolite("glc__D_c", compartment=c)

pyr_c    = Metabolite("pyr_c", compartment=c)
lac__D_c = Metabolite("lac__D_c", compartment=c)
ac_c     = Metabolite("ac_c", compartment=c)
co2_c    = Metabolite("co2_c", compartment=c)

atp_c  = Metabolite("atp_c", compartment=c)
adp_c  = Metabolite("adp_c", compartment=c)
pi_c   = Metabolite("pi_c", compartment=c)
h2o_c  = Metabolite("h2o_c", compartment=c)
h_c    = Metabolite("h_c", compartment=c)
nad_c  = Metabolite("nad_c", compartment=c)
nadh_c = Metabolite("nadh_c", compartment=c)

model.add_metabolites([
    glc__D_e, glc__D_c,
    pyr_c, lac__D_c, ac_c, co2_c,
    atp_c, adp_c, pi_c, h2o_c, h_c, nad_c, nadh_c
])

# Exchanges / supplies
add_exchange(model, glc__D_e, "EX_glc__D_e", lb=-10.0, ub=1000.0)

# Toy “base pool” supplies
add_exchange(model, h2o_c, "EX_h2o_c")
add_exchange(model, pi_c,  "EX_pi_c")
add_exchange(model, adp_c, "EX_adp_c")
add_exchange(model, nad_c, "EX_nad_c")
add_exchange(model, h_c,   "EX_h_c")

# Transport
GLCtex = Reaction("GLCtex")
GLCtex.lower_bound = 0.0
GLCtex.upper_bound = 1000.0
GLCtex.add_metabolites({glc__D_e: -1.0, glc__D_c: 1.0})
model.add_reactions([GLCtex])

# Glycolysis (lumped)
GLY = Reaction("GLYCOLYSIS")
GLY.lower_bound = 0.0
GLY.upper_bound = 1000.0
GLY.add_metabolites({
    glc__D_c: -1.0,
    adp_c: -2.0,
    pi_c: -2.0,
    nad_c: -2.0,
    pyr_c: 2.0,
    atp_c: 2.0,
    nadh_c: 2.0,
    h2o_c: 2.0,
    h_c: 2.0
})
model.add_reactions([GLY])

# Fermentation / byproducts (optional routes)
LDH = Reaction("LDH_D")
LDH.lower_bound = 0.0
LDH.upper_bound = 1000.0
LDH.add_metabolites({pyr_c: -1.0, nadh_c: -1.0, h_c: -1.0, lac__D_c: 1.0, nad_c: 1.0})
model.add_reactions([LDH])

ACET = Reaction("ACETATE_PROD")
ACET.lower_bound = 0.0
ACET.upper_bound = 1000.0
ACET.add_metabolites({pyr_c: -1.0, nad_c: -1.0, ac_c: 1.0, co2_c: 1.0, nadh_c: 1.0})
model.add_reactions([ACET])

# Maintenance
ATPM = Reaction("ATPM")
ATPM.lower_bound = 0.0
ATPM.upper_bound = 1000.0
ATPM.add_metabolites({atp_c: -1.0, h2o_c: -1.0, adp_c: 1.0, pi_c: 1.0, h_c: 1.0})
model.add_reactions([ATPM])

# Byproduct sinks
add_sink(model, lac__D_c, "DM_lac__D_c")
add_sink(model, ac_c,     "DM_ac_c")
add_sink(model, co2_c,    "DM_co2_c")

# Biomass as a DRAIN reaction (no biomass metabolite)
# pyr + ATP + NADH -> ADP + Pi + NAD + H
BIOMASS = Reaction("BIOMASS_Ecoli_core")
BIOMASS.lower_bound = 0.0
BIOMASS.upper_bound = 1000.0
BIOMASS.add_metabolites({
    pyr_c: -1.0,
    atp_c: -1.0,
    nadh_c: -1.0,
    adp_c: 1.0,
    pi_c: 1.0,
    nad_c: 1.0,
    h_c: 1.0
})
model.add_reactions([BIOMASS])

model.objective = "BIOMASS_Ecoli_core"

sol = model.optimize()
print("Status:", sol.status)
print("Objective (biomass flux):", sol.objective_value)
for rid in ["EX_glc__D_e","GLCtex","GLYCOLYSIS","BIOMASS_Ecoli_core","LDH_D","ACETATE_PROD","ATPM"]:
    print(f"{rid:18s} {sol.fluxes[rid]: .6f}")
