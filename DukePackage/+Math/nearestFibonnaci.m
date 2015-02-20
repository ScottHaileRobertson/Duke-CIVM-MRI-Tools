% Description: This function will return the closest fibonacci number to a
% given starting number. This is useful for picking a number of 3d golden
% means vectors as it will give almost uniform sampling
%
% Author: Scott Haile Robertson
function nearest_fibonnaci = nearestFibonnaci(startingNumber)
current_fibonacci = 1;
last_fibonacci = 0;

while current_fibonacci < startingNumber
    next_fibonacci=current_fibonacci+last_fibonacci;
    
    % Prepare for next itteration
    last_fibonacci = current_fibonacci;
    current_fibonacci = next_fibonacci;
    % 
end

% Figure out the closest fibonacci number
if(current_fibonacci == startingNumber || last_fibonacci == startingNumber)
    nearest_fibonnaci = startingNumber;
elseif(abs(current_fibonacci-startingNumber)<=abs(last_fibonacci-startingNumber))
    nearest_fibonnaci  = current_fibonacci;
else
    nearest_fibonnaci  = last_fibonacci;
end