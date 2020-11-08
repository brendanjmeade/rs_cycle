using DifferentialEquations
using PyCall
using PyPlot


function aginglaw(v, theta, dc)
    return 1 - theta * v / dc
end


function sliplaw(v, theta, dc)
    return -v * theta / dc * log(v * θ / dc)
end


function aofv(v, Va, a0, Sa)
    return a0 + Sa * log10((Va + v) / Va)
end


function dcofv(v, Vdc, dc0, Sdc)
    return dc0 + Sdc * log10((Vdc  +v) / Vdc)
end


function calcdvθclassic!(du, u, p, t)
    dc, η, σn, a, b, μ, Vp, L, ρ, statelaw = p
    θ = u[1]
    v = u[2]
    du[1] = statelaw(v, θ, dc)
    du[2] = 1 / (η / σn + a / v) * (μ * (Vp - v) / (L * σn) - b * du[1] / θ)
    return nothing
end


function calcdvθveldep!(du, u, p, t)
    η, σn, b, μ, vp, L, ρ, Va, a0, Sa, Vdc, dc0, Sdc, statelaw = p
    θ = u[1]
    v = u[2]
    a = aofv(v, Va, a0, Sa)
    dc = dcofv(v, Vdc, dc0, Sdc)
    du[1] = statelaw(v, θ, dc)
    du[2] = 1 / (η / σn + a / v) * (μ * (vp - v) / (L * σn) - b * du[1] / θ)
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
    close("all")
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
    
    # Time integrate - classic
    icsclassic = [1e8; vp / 1000]
    pclassic = (dc, η, σn, a, b, μ, vp, L, ρ, aginglaw)
    probclassic = ODEProblem(calcdvθclassic!, icsclassic, tspan, pclassic)
    solclassic = solve(probclassic, RK4(), abstol = abstol, reltol = reltol)
    # solclassic = solve(probclassic, RK4(), abstol = abstol, reltol = reltol)
    plottimeseries(solclassic, siay, "RSF (classic)")

    # Time integrate - a, dc evolved with v
    # Values and form from Im et al. 2020
    Va = 100e-6
    a0 = 0.005
    Sa = 0.0003 / (10 * siay)
    Vdc = 100e-6
    dc0 = 10e-6
    Sdc = 30e-6
    
    icsveldep = [1e8 ; vp / 1000]
    pveldep = (η, σn, b, μ, vp, L, ρ, Va, a0, Sa, Vdc, dc0, Sdc, aginglaw)
    probveldep = ODEProblem(calcdvθveldep!, icsveldep, tspan, pveldep)
    solveldep = solve(probveldep, RK4(), abstol = abstol, reltol = reltol)
    plottimeseries(solveldep, siay, "RSF (v)")

    return nothing
end
sliding()
