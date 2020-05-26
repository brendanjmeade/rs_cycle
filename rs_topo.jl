using DifferentialEquations
using PyCall
using PyPlot

function plottimeseries(sol, siay, titlelabel)
    t1 = [x / siay for x in sol.t]
    θ1 = [x[1] for x in sol.u]
    v1 = [x[2] for x in sol.u]
    σ1 = [x[3] for x in sol.u]
    τ1 = [x[4] for x in sol.u]
    h1 = [x[5] for x in sol.u]

    figure(figsize = (12, 18))
    subplot(5, 2, 1)
    plot(t1, θ1, "-b", linewidth = 0.5)
    yscale("log")
    ylabel(L"\theta")
    subplot(5, 2, 2)
    plot(1:1:length(t1), θ1, "-b", linewidth = 0.5)
    yscale("log")
    ylabel(L"\theta")

    subplot(5, 2, 3)
    plot(t1, v1, "-b", linewidth = 0.5)
    yscale("log")
    xlabel("time (years)")
    ylabel("v (m/s)")
    subplot(5, 2, 4)
    plot(1:1:length(t1), v1, "-b", linewidth = 0.5)
    yscale("log")
    xlabel("time (step)")
    ylabel("v (m/s)")
    suptitle(titlelabel)

    subplot(5, 2, 5)
    plot(t1, σ1, "-b", linewidth = 0.5)
    yscale("log")
    xlabel("time (years)")
    ylabel("normal stress (Pa)")
    subplot(5, 2, 6)
    plot(1:1:length(t1), σ1, "-b", linewidth = 0.5)
    yscale("log")
    xlabel("time (step)")
    ylabel("normal stress (Pa)")

    subplot(5, 2, 7)
    plot(t1, τ1, "-b", linewidth = 0.5)
    yscale("log")
    xlabel("time (years)")
    ylabel("shear stress (Pa)")
    subplot(5, 2, 8)
    plot(1:1:length(t1), τ1, "-b", linewidth = 0.5)
    yscale("log")
    xlabel("time (step)")
    ylabel("shear stress (Pa)")

    subplot(5, 2, 9)
    plot(t1, h1, "-b", linewidth = 0.5)
    # yscale("log")
    xlabel("time (years)")
    ylabel("height (m)")
    subplot(5, 2, 10)
    plot(1:1:length(t1), h1, "-b", linewidth = 0.5)
    # yscale("log")
    xlabel("time (step)")
    ylabel("height (m)")

    suptitle(titlelabel)
    return nothing
end

function aginglaw(v, theta, dc)
    return 1 - theta * v / dc
end

function sliplaw(v, theta, dc)
    return -v * theta / dc * log(v * theta / dc)
end

# function derivs!(du, u, p, t)
#     dc, η, σn, a, b, μ, Vp, L, ρ, statelaw = p
#     θ = u[1]
#     v = u[2]
#     du[1] = statelaw(v, θ, dc)
#     du[2] = 1 / (η / σn + a / v) * (μ * (Vp - v) / (L * σn) - b * du[1] / θ)
#     return nothing
# end

function derivs!(dudt, u, p, t)
    dc, η, a, b, μ, Vp, L, ρ, statelaw, α, β, g = p
    θ = u[1]
    v = u[2]
    σ = u[3]
    τ = u[4]
    h = u[5] 
    dτdt = μ*(Vp-v)/L # + 0.1*dudt[3]
    # dτdt = μ*(Vp-0)/L # + 0.1*dudt[3]

    dudt[1] = statelaw(v, θ, dc)
    dudt[2] = 1/(η/σ + a/v) * (dτdt/σ - b*dudt[1]/θ + (η*v-τ)*dudt[3]/σ^2)
    dudt[3] = ρ*g*dudt[5]
    dudt[4] = dτdt
    dudt[5] = α*v - β*h^2
    return nothing
end


function rstopo()
    #! Model parameters
    siay = 365.25 * 24 * 60 * 60
    tspan = (0.0, siay * 5000.0)
    μ = 3e10
    ν = 0.25
    ρ = 2700.0
    η = μ / (2.0 * sqrt(μ/ρ))
    L = 60 * 1e3
    a = 0.015
    b = 0.02
    α = sind(45) # Should be given by trig function based on dip of fault
    # β = 1e-8 # Check Whipple and Meade (2006)...should be small  THIS WORKED!!!
    # β = 1e-17 # Check Whipple and Meade (2006)...should be small
    β = 1e-9 # Check Whipple and Meade (2006)...should be small

    vp = 1e-9
    σn = 30e6
    g = 9.81
    h = σn / (ρ*g)
    dc = 0.1
    abstol = 1e-4
    reltol = 1e-4
    
    #! Time integrate
    icsclassic = [1e8; vp / 1000; σn; 0; h]
    # pclassic = (dc, η, a, b, μ, vp, L, ρ, aginglaw, α, β, g)
    pclassic = (dc, η, a, b, μ, vp, L, ρ, aginglaw, α, β, g)

    probclassic = ODEProblem(derivs!, icsclassic, tspan, pclassic)
    # solclassic = solve(probclassic, RK4(), abstol=abstol, reltol=reltol)
    solclassic = solve(probclassic, DP8(), abstol=abstol, reltol=reltol)


    #! Plot results
    close("all")
    plottimeseries(solclassic, siay, "RSF (classic)")

    return nothing
end
rstopo()
