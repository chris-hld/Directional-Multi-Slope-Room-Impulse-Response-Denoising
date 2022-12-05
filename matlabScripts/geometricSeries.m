function gseries = geometricSeries(r, n)
%GEOMETRICSERIES The sum of the first n terms of a geometric series
gseries = (1 - r.^(n+1))./(1-r);
end