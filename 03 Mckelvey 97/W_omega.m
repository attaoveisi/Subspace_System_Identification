function output = W_omega(omega,q)
output = zeros(1,q);
for i = 1:q
    output(1,i) = exp((i-1)*1i*omega);
end
output = output';

