function dydt = odesystem(t,y,z)
%ODESYSTEM := this function defines the ode system at time t.
%   inputs:
%       - t <- time
%       - y <- states
%       - z <- parameters
%   outputs:
%       - dydt <- ode system at time t

% state order:
% y := [(1) S_1, (2) I_1, (3) R_1, (4) S_2, (5) I_2, (6) R_2, (7) H, (8) D];

% parameter order:
% z := [(1) beta, (2) alpha, (3) phi_12^max, (4) phi_21^max, (5) gamma_I, (6) eta, (7) delta_I, (8) p, (9) gamma_H, (10) delta_H, (11) a_1, (12) a_2]

    function phi_i = phi(z)
        phi_i = -1.0*ones(2,1);
        phi_i(1) = (z(3)*(y(4) + y(5) + y(6)))/(z(11) + (y(4) + y(5) + y(6))); %phi_12 using nonconstant Holling Type II
        phi_i(2) = (z(4)*y(8))/(z(12) + y(8)); %phi_21 using nonconstant Holling Type II
    end

I_M = (1-z(2))*y(2) + y(5);
phi_ij = phi(z);

dydt = -1.0*ones(8,1);

dydt(1) = -z(1)*(1-z(2))*I_M - phi_ij(1)*y(1) + phi_ij(2)*y(4); % dS_1
dydt(2) = z(1)*(1-z(2))*I_M - (z(5) + z(6) + z(7) + phi_ij(1))*y(2) + phi_ij(2)*y(5); % dI_1
dydt(3) = z(5)*y(2) + z(8)*z(9)*y(7) - phi_ij(1)*y(3) + phi_ij(2)*y(6); % dR_1
dydt(4) = -z(1)*y(4)*I_M + phi_ij(1)*y(1) - phi_ij(2)*y(4); % dS_2
dydt(5) = z(1)*y(4)*I_M - (z(5) + z(6) + z(7) + phi_ij(2))*y(5) + phi_ij(1)*y(2); % dI_2
dydt(6) = z(5)*y(5) + (1-z(8))*z(9)*y(7) + phi_ij(1)*y(3) - phi_ij(2)*y(6); % dR_2
dydt(7) = z(6)*(y(2) + y(5)) - z(9)*y(7) - z(10)*y(7); % dH
dydt(8) = z(7)*(y(2) + y(5)) + z(10)*y(7); % dD

end