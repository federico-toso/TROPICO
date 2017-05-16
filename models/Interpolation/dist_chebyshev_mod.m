function [nodes] = dist_chebyshev_mod (ti, tf, nc)

% this chebyshev distribution shifts the start and end points instead of
% adding them after distributing the nodes included

if nc == 1
    nodes=0;
    return
end

k2 = nc:-1:1;
nodes = cos(((2*k2-1)./(2*nc))*pi);
nodes = ti+((nodes+1)/(nodes(end)-nodes(1))-(nodes(1)+1)/(nodes(end)-nodes(1)))*(tf-ti);