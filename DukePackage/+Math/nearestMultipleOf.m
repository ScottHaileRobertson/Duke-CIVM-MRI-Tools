function multipleOf = nearestMultipleOf(x, multipleOf)
    multipleOf = round(x/multipleOf)*multipleOf; % Same as floor((x/multipleOf)+0.5)*multipleOf in EPIC
end