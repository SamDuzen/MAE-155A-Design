
%Parameter input
delta = 10; %deg
Mo = 1.6;
gamma = 1.4;

%Evaluate Left Side
LHS = tand(delta);

%Initialize
error(1) = 20;
theta(1) = 0; %initial guess
convergence = 0.0001; %convergence threshold
i = 1; %initialize

while error(i) >= convergence
    %Compute Right Side
        RHS = (2*cotd(theta(i))*(Mo^2*sind(theta(i))^2 - 1))/(2 + Mo^2*(gamma + 1 - 2*sind(theta(i))^2));
    %Compute Error
        error(i+1) = abs(RHS-LHS);
    %Update Parameters
        if abs(RHS) < abs(LHS)
            theta(i+1) = theta(i) + convergence;
        end
        if abs(RHS) > abs(LHS)
            theta(i+1) = theta(i) - convergence;
        end
        i = i+1;

end

theta(end)
