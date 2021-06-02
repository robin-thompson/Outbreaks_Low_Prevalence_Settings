%% Case Outbreak Risk ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dYdt = COR_ode(t,Y,beta,ti,mu)

dYdt = mu*(max(beta_int(t,beta,ti)/mu,1)*(Y-Y^2) + (Y-1));

function X = beta_int(t,beta,ti)

X = interp1(ti,beta,t);
    
end

end