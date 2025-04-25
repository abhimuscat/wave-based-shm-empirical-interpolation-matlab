%% Double-Gaussian pulse function
function output = f(xx, yy, mumu)
    a_mu  = @(mu) 10 + 5*sin(mu) + 2*cos(2*mu) + 0.5*exp(0.1*mu) + 0.3*sin(xx) + 0.2*cos(yy.^2);
    b_mu  = @(mu) 10 + 5*cos(mu) + 2*sin(2*mu) + 0.5*exp(0.1*mu) + 0.2*cos(xx.^2) + 0.3*sin(yy);
    c_mu  = @(mu) 15 - 5*sin(mu) + 2*cos(2*mu) + 0.5*exp(0.1*mu) + 0.2*sin(xx) + 0.3*cos(yy.^2);
    d_mu  = @(mu) 15 - 5*cos(mu) + 2*sin(2*mu) + 0.5*exp(0.1*mu) + 0.3*cos(xx.^2) + 0.2*sin(yy);
    s1    = @(mu) 5 + 2*sin(mu) + cos(mu) + 0.5*exp(0.05*mu^2);
    s2    = @(mu) 5 + 2*cos(mu) + sin(mu) + 0.5*exp(0.05*mu^2);
    h1    = @(mu) abs(0.4 + 4.6*(sin(10*mu) + cos(8*mu)));
    h2    = @(mu) abs(0.4 + 4.6*(cos(12*mu) + sin(6*mu)));

    G1 = h1(mumu) .* exp(-((xx - a_mu(mumu)).^2 + (yy - b_mu(mumu)).^2) / (2 * s1(mumu)^2));
    G2 = h2(mumu) .* exp(-((xx - c_mu(mumu)).^2 + (yy - d_mu(mumu)).^2) / (2 * s2(mumu)^2));

    output = (G1 + G2) / max(G1 + G2, [], 'all');
end
