function prime_number = nearestPrime(num)
% If the number is prime, just return
if(isprime(num))
    prime_number = num;
    return;
end

% Otherwise find the nearest prime number
i=1;
while (i>0)
    above_prime = isprime(num+i);
    below_prime = isprime(num-i);
    
    if(above_prime & below_prime)
        prime_number = num+i;
        return;
    elseif(above_prime)
        prime_number = num+i;
        return;
    elseif(below_prime)
        prime_number = num-i;
        return;
    end
    i=i+1;
end
