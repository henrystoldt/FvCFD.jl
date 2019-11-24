# julia shock tube solver
#Written with conservative Euler eqns, FD MacCormack method
#include PyPlot
using PyPlot
#using Plots.PlotMeasures
#pyplot()



function main()

    #General Constants
    #R = 8.314462618
    R = 287.05
    gamma = 1.4

    #Define the grid
    n = 1001
    n_half = 501
    x_pos = range(0,stop=1,length=n)

    #Time step size
    t_end = 0.14267
    #t_end = 0.5
    #dt = t_end / 10000
    #t_end = 0.14267
    dt = t_end/8000

    #Check time_step size is appropriate
    #For now, impose a soft cfl limit of 0.2
    a = speed_of_sound(0.15, R, gamma)
    dx = x_pos[2]-x_pos[1]
    CFL = (a+1) * dt/dx

    if CFL > sqrt(3)/2
        println("Max CFL on uniform grid is : ", CFL)
        println("You need a smaller time step bozo!")
        return
    else
        println("CFL is : ", CFL )
        println("CFL looks good!")
    end

    #Artificial Viscosity Constants
    C_x = 0.3
    

    

    #Arrays to store the conditions
    #First primitives
    rho = zeros(n)
    u_vel = zeros(n)
    ener = zeros(n)

    temp = zeros(n)
    pres = zeros(n)

    #Conserved quantities
    mass_flux = zeros(n)
    u_mom_flux = zeros(n)
    ener_flux = zeros(n)

    #Residuals

    #Time




    #Setup initial conditions
    for i in 1:n
        if i <= n_half
            #LHS of tube
            temp[i] = 0.00348432
            #temp[i] = 0.120279
            pres[i] = 1
        else
            #RHS of tube
            temp[i] = 0.00278746
            #temp[i] = 0.0962223
            pres[i] = 0.1
        end

        #Same for either side
        rho[i] = pres[i]/(R*temp[i])
        u_vel[i] = 0
        ener[i] = (1/(gamma-1)) * pres[i]/rho[i]

        #Not using these two
        #u_mom_flux[i] = 0
        #ener_flux[i] = rho[i] * ener[i]
    end

    #=

    #println("This is temp:", temp)
    speed_os = speed_of_sound(temp, R, gamma)

    min_sos = minimum(speed_os)

    println("This is sos ", min_sos)
    delta_x = x_pos[2]-x_pos[1]

    print("This is d_x ", delta_x, "\n")
    =#


    results = do_cfd(rho, u_vel, ener, pres, temp, R, gamma, t_end, n , dx, dt, C_x)


    println("Finshed solving CFD, about to plot!")

    rho_res = results[1]
    u_res = results[2]
    ener_res = results[3]
    pres_res = results[4]
    temp_res = results[5]

    #println("Rho results are: ", rho_res)

    plot(x_pos, rho_res, label="RHO")
    #title("Density Results at the final time")
    #show()
    
    #plot(x_pos, pres_res, marker="o", label="PRES")
    #title("Pressure Results at the final time")
    #show()

    
    #plot(x_pos, temp_res, marker="o", label="TEMP")
    #title("Temperature Results at the final time")
    #show()

    #plot(x_pos, u_res, marker="o", label="VELO")
    #title("Velocity Results at the final time")

    #show()
    


    #plot(x_pos)
    #gui()

end

function do_cfd(rho, u_vel, ener, pres, temp, R, gamma, t_end, n, dx, dt, C_x)
    #=Performs CFD
    Might be some machine learning involved too?
    =#

    #For the moment not storing a complete time history, only store the last time step

    #Loop through time
    t = 0
    while t<t_end
        #Predict step
        predicted_vals = predictor(rho, u_vel, ener, pres, R, gamma, dx, dt, C_x)
        #Unpack predicted vals
        drho_dt_pred = predicted_vals[1]
        dmom_dt_pred = predicted_vals[2]
        dener_dt_pred = predicted_vals[3]
        rho_pred = predicted_vals[4]
        u_pred = predicted_vals[5]
        e_pred = predicted_vals[6]
        pres_pred = predicted_vals[7]
        temp_pred = predicted_vals[8]

        mom_pred = rho_pred .* u_pred
        ener_flux_pred = rho_pred .* (e_pred + 0.5 .* u_pred .* u_pred)


        #Correct
        cor_slopes = corrector(rho_pred, u_pred, e_pred, pres_pred, dx)

        drho_dt_cor = cor_slopes[1]
        dmom_dt_cor = cor_slopes[2]
        dener_dt_cor = cor_slopes[3]


        #Average slopes
        drho_dt_avg = 0.5 .* (drho_dt_pred + drho_dt_cor)
        dmom_dt_avg = 0.5 .* (dmom_dt_pred + dmom_dt_cor)
        dener_dt_avg = 0.5 .* (dener_dt_pred + dener_dt_cor)

        #Calculate corrected artificial viscosity
        s_rho_cor = zeros(n)
        s_mom_cor = zeros(n)
        s_ener_cor = zeros(n)

        for i in 2:(n-1)
            s_rho_cor[i] = artificial_visc(pres_pred[i+1],pres_pred[i],pres_pred[i-1], rho_pred[i+1], rho_pred[i], rho_pred[i-1], C_x)
            s_mom_cor[i] = artificial_visc(pres_pred[i+1],pres_pred[i],pres_pred[i-1], mom_pred[i+1], mom_pred[i], mom_pred[i-1], C_x)
            s_ener_cor[i] = artificial_visc(pres_pred[i+1],pres_pred[i],pres_pred[i-1], ener_flux_pred[i+1], ener_flux_pred[i], ener_flux_pred[i-1], C_x)
        end

        #Update quantities
        #First conserved quantities
        new_rho = rho + drho_dt_avg .* dt + s_rho_cor
        new_u_mom = rho .* u_vel + dmom_dt_avg .* dt + s_mom_cor
        new_ener_flux = rho .* (ener + 0.5 .*u_vel.*u_vel) + dener_dt_avg .* dt + s_ener_cor

        #Now new primitives
        new_u_vel = new_u_mom ./ new_rho
        new_ener = (new_ener_flux ./ new_rho) - 0.5 .* new_u_vel .* new_u_vel
        new_pres = new_ener .* new_rho * (gamma-1)
        new_temp = new_ener .* (gamma-1) / R


        #Reset values for next time_step
        rho = new_rho
        u_vel = new_u_vel
        ener = new_ener
        pres = new_pres
        temp = new_temp

        t += dt
    end

    return [rho, u_vel, ener, pres, temp]


end


function predictor(rho, u_vel, ener, pres, R, gamma, dx, dt, C_x)
    #=Takes input of current variables for all cells
    Performs one explicit time step using forward difference methods

    Returns:
        -Predicted values of all primitives (rho, u, e) for all cells
        -Predicted values of pressure
        -Predicted slopes of fluxes(rho, mom, ener_flux)
    =#

    #Determine length of number of mesh cells
    n = length(rho)

    #Calculate the current flux values
    mom_flux = rho .* u_vel
    ener_flux = rho .* (ener + 0.5 .* u_vel .* u_vel)

    #Store predicted gradients
    drho_dt_pred = zeros(n)
    dmom_dt_pred = zeros(n)
    dener_dt_pred = zeros(n)

    #Store predicted updated conserved variables
    rho_pred = zeros(n)
    u_mom_pred = zeros(n)
    ener_flux_pred = zeros(n)

    #Store predicted updated primitives
    u_pred = zeros(n)
    e_pred = zeros(n)
    pres_pred = zeros(n)
    temp_pred = zeros(n)

    #For all cells, feed in input values to conservation Equations
    for i in 1:n #Loop over all cells
        if i==1
            #Use first BCs
            #Zero grad everything is hardcoded in for now
            drho_dt_pred[i] = mass_conserv(rho[i], rho[i], u_vel[i], u_vel[i], dx)
            dmom_dt_pred[i] = mom_conserv(rho[i], rho[i], u_vel[i], u_vel[i], pres[i], pres[i], dx)
            dener_dt_pred[i] = ener_conserv(rho[i], rho[i], u_vel[i], u_vel[i], pres[i], pres[i], ener[i], ener[i], dx)

            #No viscosity at boundary because dP/dx =0
            s_rho_pred = 0
            s_mom_pred = 0
            s_ener_pred = 0
        elseif i==n
            #Use last BCs
            #Zero grad everything is hardcoded in for now
            drho_dt_pred[i] = mass_conserv(rho[i], rho[i], u_vel[i], u_vel[i], dx)
            dmom_dt_pred[i] = mom_conserv(rho[i], rho[i], u_vel[i], u_vel[i], pres[i], pres[i], dx)
            dener_dt_pred[i] = ener_conserv(rho[i], rho[i], u_vel[i], u_vel[i], pres[i], pres[i], ener[i], ener[i], dx)

            s_rho_pred = 0
            s_mom_pred = 0
            s_ener_pred = 0
        else
            #Use actual conditions
            drho_dt_pred[i] = mass_conserv(rho[i], rho[i+1], u_vel[i], u_vel[i+1], dx)
            dmom_dt_pred[i] = mom_conserv(rho[i], rho[i+1], u_vel[i], u_vel[i+1], pres[i], pres[i+1], dx)
            dener_dt_pred[i] = ener_conserv(rho[i], rho[i+1], u_vel[i], u_vel[i+1], pres[i], pres[i+1], ener[i], ener[i+1],dx)

            s_rho_pred = artificial_visc(pres[i+1], pres[i], pres[i-1], rho[i+1], rho[i], rho[i-1], C_x)
            s_mom_pred = artificial_visc(pres[i+1], pres[i], pres[i-1], mom_flux[i+1], mom_flux[i], mom_flux[i-1], C_x)
            s_ener_pred = artificial_visc(pres[i+1], pres[i], pres[i-1], ener_flux[i+1], ener_flux[i], ener_flux[i-1], C_x)
        end

        

        #Now create predicted updated variables (conserved)
        rho_pred[i] = rho[i] + drho_dt_pred[i] * dt + s_rho_pred
        u_mom_pred[i] = rho[i]*u_vel[i] + dmom_dt_pred[i] * dt + s_mom_pred
        ener_flux_pred[i] = rho[i]*(ener[i] + 0.5*u_vel[i]*u_vel[i]) + dener_dt_pred[i] * dt + s_ener_pred

        #Convert conserved quantities into primitives
        u_pred[i] = u_mom_pred[i] / rho_pred[i]
        e_pred[i] = ener_flux_pred[i]/rho_pred[i] - 0.5 * (u_pred[i])^2
        pres_pred[i] = e_pred[i]*rho_pred[i]*(gamma-1)
        temp_pred[i] = e_pred[i]*(gamma-1)/R

    end

    return [drho_dt_pred, dmom_dt_pred, dener_dt_pred, rho_pred, u_pred, e_pred, pres_pred, temp_pred]

end

function corrector(rho_pred, u_pred, ener_pred, pres_pred, dx)
    #=Does the corrector step, using backwards differencing
    Takes input of predicted values

    Returns corrected flux slopes
    =#

    #Determine length of number of mesh cells
    n = length(rho_pred)

    #Store predicted gradients
    drho_dt_cor = zeros(n)
    dmom_dt_cor = zeros(n)
    dener_dt_cor = zeros(n)

    for i in 1:n
        if i==1
            #Use first BCs
            #Zero grad everything is hardcoded in for now
            drho_dt_cor[i] = mass_conserv(rho_pred[i], rho_pred[i], u_pred[i], u_pred[i], dx)
            dmom_dt_cor[i] = mom_conserv(rho_pred[i], rho_pred[i], u_pred[i], u_pred[i], pres_pred[i], pres_pred[i], dx)
            dener_dt_cor[i] = ener_conserv(rho_pred[i], rho_pred[i], u_pred[i], u_pred[i], pres_pred[i], pres_pred[i], ener_pred[i], ener_pred[i], dx)
        elseif i==n
            #Use last BCs
            #Zero grad everything is hardcoded in for now
            drho_dt_cor[i] = mass_conserv(rho_pred[i], rho_pred[i], u_pred[i], u_pred[i], dx)
            dmom_dt_cor[i] = mom_conserv(rho_pred[i], rho_pred[i], u_pred[i], u_pred[i], pres_pred[i], pres_pred[i], dx)
            dener_dt_cor[i] = ener_conserv(rho_pred[i], rho_pred[i], u_pred[i], u_pred[i], pres_pred[i], pres_pred[i], ener_pred[i], ener_pred[i], dx)
        else
            #Use actual conditions
            drho_dt_cor[i] = mass_conserv(rho_pred[i], rho_pred[i-1], u_pred[i], u_pred[i-1], dx, dif_direc=1)
            dmom_dt_cor[i] = mom_conserv(rho_pred[i], rho_pred[i-1], u_pred[i], u_pred[i-1], pres_pred[i], pres_pred[i-1], dx, dif_direc=1)
            dener_dt_cor[i] = ener_conserv(rho_pred[i], rho_pred[i-1], u_pred[i], u_pred[i-1], pres_pred[i], pres_pred[i-1], ener_pred[i], ener_pred[i-1], dx, dif_direc=1)
        end
    end

    return [drho_dt_cor, dmom_dt_cor, dener_dt_cor]

end

# -------------------------------------------------------------------
# Euler Equations
# 
# - Equations should only take scalar inputs, no vectors

function mass_conserv(rho_1, rho_2, u_1, u_2, dx; dif_direc=0)
    #=Performs the mass conservation equation 

    rho_1,u_1 should always be i
    rho_2,u_2 will be (i+1) if forward difference, and (i-1) if backward dif

    dif_direc is the flag indicating which direction the differencing should be
    if dif_direc==0
        forwards difference
    if dif_direc==1
        backwatds dif
    =#


    if dif_direc==0
        #Forward difference
        drho_dt = -1 * (rho_1 * (u_2-u_1) + u_1 * (rho_2-rho_1)) / dx
    elseif dif_direc==1
        #Backwards difference
        drho_dt = -1 * (rho_1 * (u_1-u_2) + u_1 * (rho_1-rho_2)) / dx
    else
        println("Can't take input other than 0 or 1 for dif direction!")
    end

    return drho_dt
end

function mom_conserv(rho_1, rho_2, u_1, u_2, p_1, p_2, dx; dif_direc=0)
    #=Solves conservative momentum conservation equation

    dif_direc is the flag indicating which direction the differencing should be
    if dif_direc==0
        forwards difference
    if dif_direc==1
        backwatds dif
    =#

    if dif_direc==0
        dmom_dt = -1 * ( (p_2-p_1) + 2*rho_1*u_1*(u_2-u_1) + u_1*u_1*(rho_2-rho_1) ) / dx
    elseif dif_direc==1
        dmom_dt = -1 * ( (p_1-p_2) + 2*rho_1*u_1*(u_1-u_2) + u_1*u_1*(rho_1-rho_2) ) / dx
    else
        println("Can't take input other than 0 or 1 for dif direction!")
    end

    return dmom_dt
end

function ener_conserv(rho_1, rho_2, u_1, u_2, p_1, p_2, e_1, e_2, dx; dif_direc=0)
    #= Solves energy conservation equation in conservative form
    dif_direc is the flag indicating which direction the differencing should be
    if dif_direc==0
        forwards difference
    if dif_direc==1
        backwatds dif
    =#

    if dif_direc==0
        dener_dt = -1 * ( p_1*(u_2-u_1) + u_1*(p_2-p_1) + rho_1*e_1*(u_2-u_1) + rho_1*u_1*(e_2-e_1) + u_1*e_1*(rho_2-rho_1) + 1.5*rho_1*u_1*u_1*(u_2-u_1) + 0.5*(u_1^3)*(rho_2-rho_1) ) / dx
    elseif dif_direc==1
        dener_dt = -1 * ( p_1*(u_1-u_2) + u_1*(p_1-p_2) + rho_1*e_1*(u_1-u_2) + rho_1*u_1*(e_1-e_2) + u_1*e_1*(rho_1-rho_2) + 1.5*rho_1*u_1*u_1*(u_1-u_2) + 0.5*(u_1^3)*(rho_1-rho_2) ) / dx
    else
        println("Can't take input other than 0 or 1 for dif direction!")
    end

    return dener_dt
end

function artificial_visc(p_iplus, p_i, p_imin, U_iplus, U_i, U_imin, C_x)
    #=This is the function that calculates the magnitude of the artificial viscosity
    =#

    visc = (C_x * abs(p_iplus - 2*p_i + p_imin)/(p_iplus + 2*p_i + p_imin)) * (U_iplus - 2*U_i + U_imin)
    return visc
end


# ----------------------------------------------------
# Helper functions

function speed_of_sound(temp, R, gamma)
    #=Takes input of temp(K), R(J/molK), gamma(unitless)
    Returns local speed of sound in m/s
    =#

    a = (gamma .* R .* temp ./ 0.02895) .^(0.5)
    return a
end

println("Got through code, about to call first time!")

main()