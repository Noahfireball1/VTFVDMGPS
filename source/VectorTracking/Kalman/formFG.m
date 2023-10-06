function F = formFG(states,timestep)

F = [0 1 0 0 0 0 0 0;...
     0 0 0 1 0 0 0 0;
     0 0 0 0 0 1 0 0;...
     0 0 0 0 0 0 0 1];

F = expm(F*timestep);

end

