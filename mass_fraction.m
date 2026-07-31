function [combined] = mass_fraction(fuel_ratio, fuel_property, water_property)
    % uses mass fractions to approximate mixture coolant properties
    combined = fuel_ratio * fuel_property + (1 - fuel_ratio) * water_property;

end