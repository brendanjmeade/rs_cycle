using DifferentialEquations
using PyCall
using PyPlot

function aginglaw(v, theta, dc)
    return 1 - theta * v / dc
end

function sliplaw(v, theta, dc)
    return -v * theta / dc * log(v * θ / dc)
end

function calcdvθclassic!(du, u, p, t)
    dc, η, σn, a, b, μ, Vp, L, ρ, statelaw = p
    θ = u[1]
    v = u[2]
    du[1] = statelaw(v, θ, dc)
    du[2] = 1 / (η / σn + a / v) * (μ * (Vp - v) / (L * σn) - b * du[1] / θ)
    return nothing
end

function plottimeseries(sol, siay, titlelabel)
    t1 = [x / siay for x in sol.t]
    θ1 = [x[1] for x in sol.u]
    v1 = [x[2] for x in sol.u]

    figure(figsize = (12, 6))
    subplot(2, 2, 1)
    plot(t1, θ1, "-b", linewidth = 0.5)
    yscale("log")
    ylabel(L"\theta")

    subplot(2, 2, 2)
    plot(1:1:length(t1), θ1, "-b", linewidth = 0.5)
    yscale("log")
    ylabel(L"\theta")

    subplot(2, 2, 3)
    plot(t1, v1, "-b", linewidth = 0.5)
    yscale("log")
    xlabel("time (years)")
    ylabel("v (m/s)")

    subplot(2, 2, 4)
    plot(1:1:length(t1), v1, "-b", linewidth = 0.5)
    yscale("log")
    xlabel("time (step)")
    ylabel("v (m/s)")
    suptitle(titlelabel)

    return nothing
end

function sliding()
    # Model parameters
    siay = 365.25 * 24 * 60 * 60
    tspan = (0.0, siay * 1000.0) # This will run
    # tspan = (0.0, siay * 302.0) # This will run but fail for the a, b case    
    μ = 3e10
    ν = 0.25
    ρ = 2700.0
    η = μ / (2.0 * sqrt(μ / ρ))
    L = 60 * 1e3
    a = 0.015
    b = 0.02
    vp = 1e-9
    σn = 30e6
    dc = 0.1
    abstol = 1e-4
    reltol = 1e-4

    # Parameters need for d(a, b σn)/dt
    fstar = 0.6 # ???
    vstar = vp # ???
    θstar = 1e9 # ???

    # Show Tushar scalings for a and b
    # tusharabmodels()
    
    # Time integrate - classic
    icsclassic = [1e8; vp / 1000]
    pclassic = (dc, η, σn, a, b, μ, vp, L, ρ, aginglaw)
    probclassic = ODEProblem(calcdvθclassic!, icsclassic, tspan, pclassic)
    solclassic = solve(probclassic, RK4(), abstol = abstol, reltol = reltol)
    plottimeseries(solclassic, siay, "RSF (classic)")

    # Time integrate - a, b can evolve
    # icsab = [1e8 ; vp / 1000 ; a ; b]
    # pab = (dc, η, σn, μ, vp, L, ρ, aginglaw, fstar, vstar, θstar)
    # probab = ODEProblem(calcdvθab!, icsab, tspan, pab)
    # solab = solve(probab, RK4(), abstol = abstol, reltol = reltol)
    # plottimeseries(solab, siay, "RSF (ab)")
    
    return nothing
end
sliding()
