
"""
function get_response_type(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage})

Returns temperature response type for appropriate dispatch of temperature-sensitive functions.
"""
function get_response_type(node::Node, species::Type{<:Species}, life_stage::Type{<:LifeStage})
return typeof(node.organisms[species].all_stages[life_stage].μ_temperature_response)
end

"""
function temperature_responsive_inheritence(response_type::Type{NoResponse}, ctemp::Float64, genetics::Genetics, gene_index::Int64)

Dispatches on `NoResponse` type to return offspring likelihoods for each genotype that are not responsive to temperature.
"""
function temperature_responsive_inheritence(response_type::Type{NoResponse}, ctemp::Float64, genetics::Genetics, gene_index::Int64)
return genetics.likelihood[:,:,gene_index]
end

"""
function temperature_responsive_inheritence(response_type::Type{T}, ctemp::Float64, genetics::Genetics, gene_index::Int64) where {T <:TemperatureResponse}

Dispatches on `TemperatureResponse` type(s) to return offspring likelihoods for each genotype that are responsive to temperature.
"""
function temperature_responsive_inheritence(response_type::Type{T}, ctemp::Float64, genetics::Genetics, gene_index::Int64) where {T <:TemperatureResponse}
# TODO citation: Fig 7B + associated data for Ross et al (2019)
threshold_options = [31.5, 33.0, 35.0] #sensitivity analysis 
threshold = threshold_options[3]
if ctemp >= threshold
genetics.likelihood .= [0.0;0.0;;0.0;0.0 ;;; 1.0;1.0;;1.0;1.0]
else
genetics.likelihood .= [1.0;1.0;;1.0;0.0 ;;; 0.0;0.0;;0.0;1.0]
end
return genetics.likelihood
end

"""
function temperature_responsive_hatching(response_type::Type{NoResponse}, ctemp::Float64, genetics::Genetics, gene_index::Int64)

Dispatches on `NoResponse` type to return hatch proportions for each genotype that are not responsive to temperature.
"""
function temperature_responsive_hatching(response_type::Type{NoResponse}, ctemp::Float64, genetics::Genetics, gene_index::Int64)
return genetics.Τ[:,:,gene_index]
end

"""
function temperature_responsive_hatching(response_type::Type{NoResponse}, ctemp::Float64, genetics::Genetics, gene_index::Int64)

Dispatches on `TemperatureResponse` type(s) to return hatch proportions for each genotype that are responsive to temperature.
"""
function temperature_responsive_hatching(response_type::Type{T},
ctemp::Float64, genetics::Genetics, gene_index::Int64) where {T <:TemperatureResponse}
# TODO: citation: Fig 7A + associated data for Ross et al (2019)
if gene_index == 2
return genetics.Τ[:,:,gene_index]
else
hatch = 1.0
if ctemp <= 26.0
    hatch = 0.91833
elseif ctemp > 36.0
    hatch = 0.0
else
    hatch = 0.9157166084543669 - 0.01950769717539987*(ctemp - 26.0) - 0.006212541487762658*exp(0.4741354738481292*(ctemp - 26.0))
end
updated_tau = genetics.Τ[:,:,gene_index]
updated_tau = [hatch hatch; 0.0 hatch]
end

# TODO citation: Fig 7B + 7C of Ross et al (2019) 
threshold_options = [31.5, 33.0, 35.0] #sensitivity analysis 
threshold = threshold_options[3]
if ctemp >= threshold && gene_index == 1
updated_tau[2,1] = 1.0
return updated_tau
elseif ctemp >= threshold && gene_index == 2
updated_tau[2,1] = 1.0
return updated_tau 
else
return updated_tau
end

end

"""
function oviposit(F, data::WolbachiaExperiment, key_species, genetics::Genetics{Wolbachia},
    gene_index::Int64, inputs::ExogenousInputs, t)

Dispatches on `WolbachiaExperiment` type(s) to return oviposited eggs, accounting for temperature-responsiveness in
hatching, cytoplasmic incompatibility, and inheritence.
"""
function oviposit(F, data::WolbachiaExperiment, key_species,
genetics::Genetics{Wolbachia},
gene_index::Int64, inputs::ExogenousInputs, t)

node = data.node
node_name = GeneDrive.get_name(node)
response_type = get_response_type(node, key_species, Egg)

ctemp = get_smoothed_temperature_value(data.smoothed_temperature_data, t,
node.temperature, inputs.temperature[node_name])

likelihood = temperature_responsive_inheritence(response_type, ctemp, genetics, gene_index)
Τ = temperature_responsive_hatching(response_type, ctemp, genetics, gene_index)
S = genetics.S
Β = genetics.Β

ΒS = Matrix{Float64}(undef, length(Β), length(S))
for i in 1:length(Β)
ΒS[i,:] = Β[i].*S'
end

O = likelihood[:,:,gene_index].*Τ.*ΒS.*F

return sum(O)
end

"""
function population_model_wolbachia(du, u, (data, inputs), t)

Builds dynamic population model that incorporates new Wolbachia-specific
temperature-responsive functions into the GeneDrive.jl ODE model.
"""
function population_model_wolbachia(du, u, (data, inputs), t)

node = data.node

##################
#   Parameters   #
##################

for (index_organism, key_species) in enumerate(keys(node.organisms))

genetics = get_genetics(node, key_species)
gN = count_genotypes(genetics)
likelihood = genetics.likelihood
S = genetics.S
Τ = genetics.Τ
Φ = genetics.Φ
Ξ_m = genetics.Ξ_m
Ξ_f = genetics.Ξ_f
Ω = genetics.Ω
Β = genetics.Β
Η = genetics.Η

n = count_substages(node, key_species)
nE = n[1]
nL = n[2]
nP = n[3]
nJuv = nE+nL+nP
nM = n[4]
nF = n[5]*gN + nM + nJuv

##################
#   State space  #
##################

E = u.x[index_organism][1:nE, :]
dE = @view du.x[index_organism][1:nE, :]

L = u.x[index_organism][1+nE : nE+nL, :]
dL = @view du.x[index_organism][1+nE : nE+nL, :]

P = u.x[index_organism][1+nE+nL : nE+nL+nP, :]
dP = @view du.x[index_organism][1+nE+nL : nE+nL+nP, :]

M = u.x[index_organism][1+nJuv, :]
dM = @view du.x[index_organism][1+nJuv, :]

F = u.x[index_organism][1+nM+nJuv : nF, :]
dF = @view du.x[index_organism][1+nM+nJuv : nF, :]

##################
#   Life stages  #
##################

for gene_index in 1:gN

    eggsnew = oviposit(F, data, key_species, genetics, gene_index, inputs, t)

    GeneDrive.create_egg!(dE, E, node, key_species, eggsnew, gene_index, inputs, t)

    GeneDrive.create_larva!(dL, L, E, node, key_species, gene_index, inputs, t)

    GeneDrive.create_pupa!(dP, P, L, node, key_species, gene_index, inputs, t)

    GeneDrive.create_male!(dM, M, P, Φ, Ξ_m, Ω, node, key_species, gene_index, inputs, t)

    matematrix = GeneDrive.mate(P, M, Φ, Ξ_f, Η, node, key_species, gene_index, inputs, t)

    GeneDrive.create_female!(dF, F, Ω, Ξ_f, node, key_species, matematrix, gene_index, inputs, t)

end

end

end

"""
function solve_dynamic_model_wolbachia(data::WolbachiaExperiment,
    releases::Vector, tempseries::Vector, algorithm, tspan)

Solves dynamic population model that incorporates new Wolbachia-specific
temperature-responsive functions into the GeneDrive.jl ODE model.
"""
function solve_dynamic_model_wolbachia(data::WolbachiaExperiment,
releases::Vector, tempseries::Vector, algorithm, tspan)

node = data.node
tstops = Vector()
callbacks = Vector()
for release in releases
    push!(tstops, release.times...)
    push!(callbacks, release.callbacks...)
end
for temp in tempseries
    push!(tstops, temp.times...)
    push!(callbacks, temp.set.discrete_callbacks...)
end
collected_callback_set = diffeq.CallbackSet((), tuple(callbacks...))
inputs = ExogenousInputs(node)
u0, dens = init_node!(node)
@info "Simulation initialized successfully."

prob = diffeq.ODEProblem(population_model_wolbachia, u0, tspan, (data, inputs))
@info "Problem built successfully."

@info "Beginning model run with **Wolbachia experimental struct**, releases, and time series temperature."
sol = diffeq.solve(prob, algorithm, callback = collected_callback_set,
tstops = unique!(tstops), reltol = 1e-9)

@info "Model run complete."
return sol
end
