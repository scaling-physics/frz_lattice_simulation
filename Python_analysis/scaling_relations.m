
for Jb=1:6
    Ja = 4;
%     Jb = 6;
    N = 102;
    NB = [1:40];
    b = 2;
    
    E1 = -Ja*(3*N - b*sqrt(N)) - Jb.*NB.*heaviside(2*sqrt(pi*N)-NB) - Jb*(2*sqrt(pi*N)).*heaviside(NB-2*sqrt(pi*N));
    hold on
    % Jb*N.*heaviside(2*sqrt(pi*N)-NB)
    plot(NB,E1)
    xlabel('FrzB num')
    ylabel('Energy of cluster')
end

% E2 = -Ja*(3*N - b*sqrt(N)) - Jb.*NB.*heaviside(2*sqrt(pi*N)-NB);

