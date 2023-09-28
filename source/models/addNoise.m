function y = addNoise(u)


y = u + randn(size(u)).*[1000;1000;1000];
end

