%This is the function which sets up the ODE problem and also what we are minimizing

function j_func = group_functional(z,data,initialvalues,tspan,w,nu)     

%Run the ode45 solver

[t,y] = ode45(@(t,y) odesystem(t,y,z),tspan,initialvalues);
dt = t(2)-t(1);

I_tot = y(:,2) + y(:,5);

diff = I_tot - reshape(data,size(I_tot));
j_func = w*0.5*norm(diff,2)^2*dt + nu*0.5*norm(z,2)^2;

end