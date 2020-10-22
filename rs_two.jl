using DifferentialEquations
using PyCall
using PyPlot

function aginglaw(v, theta, dc)
    return 1 - theta * v / dc
end


function sliplaw(v, theta, dc)
    return -v * theta / dc * log(v * theta / dc)
end


function calcdvthetaclassic!(du, u, p, t)
    dc, eta, sigman, a, b, mu, Vp, L, rho, statelaw = p
    theta = u[1]
    v = u[2]
    dtau1_dt = mu * (Vp - v) / L
    du[1] = statelaw(v, theta, dc)
    du[2] = 1 / (eta / sigman + a / v) * (dtau1_dt / sigman - b * du[1] / theta)
    return nothing
end


function calcdvthetatwo!(du, u, p, t)
    dc, eta, sigman, a, b, mu, Vp, L, rho, statelaw = p

    theta1 = u[1]
    v1 = u[2]
    theta2 = u[3]
    v2 = u[4]
    theta3 = u[5]
    v3 = u[6]

    dtau1_dt = mu * (Vp - (v1 - 0.01*v2 - 0.04*v3)) / L
    dtau2_dt = mu * (Vp - (v2 - 0.05*v1 - 0.06*v3)) / L
    dtau3_dt = mu * (Vp - (v3 - 0.03*v1 - 0.01*v2)) / L
    
    du[1] = statelaw(v1, theta1, dc)
    du[2] = 1 / (eta / sigman + a / v1) * (dtau1_dt / sigman - b * du[1] / theta1)
    du[3] = statelaw(v2, theta2, dc)
    du[4] = 1 / (eta / sigman + a / v2) * (dtau2_dt / sigman - b * du[3] / theta2)
    du[5] = statelaw(v3, theta3, dc)
    du[6] = 1 / (eta / sigman + a / v3) * (dtau3_dt / sigman - b * du[5] / theta3)
    
    return nothing
end


function rs_two()
    close("all")
    
    # Model parameters
    siay = 365.25 * 24 * 60 * 60
    tspan = (0.0, siay * 20000.0) # This will run
    mu = 3e10
    nu = 0.25
    rho = 2700.0
    eta = mu / (2.0 * sqrt(mu / rho))
    L = 60 * 1e3
    a = 0.015
    b = 0.02
    vp = 1e-9
    sigman = 30e6
    dc = 0.1
    abstol = 1e-4
    reltol = 1e-4

    # Parameters need for d(a, b simgan)/dt
    fstar = 0.6 # ???
    vstar = vp # ???
    thetastar = 1e9 # ???

    # Time integrate - classic
    icsclassic = [1e8; vp / 1000]
    pclassic = (dc, eta, sigman, a, b, mu, vp, L, rho, aginglaw)
    probclassic = ODEProblem(calcdvthetaclassic!, icsclassic, tspan, pclassic)
    solclassic = solve(probclassic, RK4(), abstol = abstol, reltol = reltol)

    # Plot velocity time series
    t1 = [x / siay for x in solclassic.t]
    theta1 = [x[1] for x in solclassic.u]
    v1 = [x[2] for x in solclassic.u]

    figure()
    plot(t1, v1, "-b", linewidth = 0.5)
    yscale("log")
    xlabel("time (years)")
    ylabel("v (m/s)")

    # Time integrate - multi
    icstwo = [1e8; vp / 1000; 1e8; vp / 1000; 1e8 ; vp / 1000]

    ptwo = (dc, eta, sigman, a, b, mu, vp, L, rho, aginglaw)
    probtwo = ODEProblem(calcdvthetatwo!, icstwo, tspan, ptwo)
    soltwo = solve(probtwo, RK4(), abstol = abstol, reltol = reltol)

    # Plot velocity time series
    t1 = [x / siay for x in soltwo.t]
    theta1 = [x[1] for x in soltwo.u]
    v1 = [x[2] for x in soltwo.u]
    theta2 = [x[3] for x in soltwo.u]
    v2 = [x[4] for x in soltwo.u]
    v3 = [x[6] for x in soltwo.u]

    figure()
    plot(t1, v1, "-b", linewidth = 0.5)
    plot(t1, v2, "--r", linewidth = 0.5)
    plot(t1, v3, ":g", linewidth = 0.5)
    
    yscale("log")
    xlabel("time (years)")
    ylabel("v (m/s)")
    
    
    return nothing
end
rs_two()
