% SRG plotting functions

G = tf([1], [1 1 1])
SRG_LTI(G)

% Plots the SRG of an LTI transfer function by computing the convex hull in hyperbolic space, then transforming back.
% The SRG is the blue line, the Nyquist diagram is the red dots.
function p = SRG_LTI(G)
    [r, i] = nyquist(G);
    nyq = squeeze(r) + 1i * abs(squeeze(i));
    nyq_p = BeltramiKlein(nyq);
    srg_p = nyq_p(convhull(real(nyq_p), imag(nyq_p)));
    srg_p = [srg_p; srg_p(end)]; % Append the last point to close the loop
    srg = BeltramiKleinInverse(srg_p);
    plot(real(nyq), imag(nyq), 'r.');
    hold on;
    for i = 1:length(srg) - 1
        am = arc_min(srg(i), srg(i+1));
        plot(real(am), imag(am), 'b-');
    end
    hold off; axis equal;
end

% Plots the SRG of an LTI transfer function in the Poincare plane.  The
% unit circle is in black, the boundary of the SRG in blue and the Nyquist
% diagram in red points.
function p = SRG_LTI_Poincare(G)
    [r, i] = nyquist(G);
    nyq = squeeze(r) + 1i * abs(squeeze(i));
    nyq_p = BeltramiKlein(nyq);
    p = scatter(real(nyq_p), imag(nyq_p), 'r.');
    hold on;
    srg_p = nyq_p(convhull(real(nyq_p), imag(nyq_p)));
    srg_p = [srg_p; srg_p(end)]; % Append the last point to close the loop
    plot(real(srg_p), imag(srg_p), 'b-');
    plot_circle();
    hold off; axis equal;
end

% Utility functions

% Plot unit circle
function plot_circle()
    t = linspace(0, 2*pi, 100);
    x = cos(t);
    y = sin(t);
    plot(x, y, 'black-');
end

% Plot arc_min between two complex numbers
function am = arc_min(z1, z2)
    if z1 == z2
        am = z1; % arcmin is single point
    elseif real(z1) == real(z2)
        am = z1 + 1i * linspace(imag(z1), imag(z2), 5); % arcmin is line
    else
        xc = (imag(z1) - imag(z2)) * (imag(z1) + imag(z2)) / (2 * (real(z1) - real(z2))) + (real(z1) + real(z2)) / 2;
        r = sqrt((real(z1) - xc)^2 + imag(z1)^2);
        theta1 = angle(z1 - xc + 1i*0);
        theta2 = angle(z2 - xc + 1i*0);
        theta = linspace(theta1, theta2, 50);
        am = xc + r * exp(1i * theta);
    end
end
 

function f = BeltramiKlein(z)
    g = (z - 1i)./(z + 1i);
    f = 2*g./(1 + abs(g).^2);
end

function z = BeltramiKleinInverse(w)
    z = (imag(w) - 1i .* sqrt(1 - w .* conj(w))) ./ (real(w) - 1);
end