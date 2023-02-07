% a method to measure the circularity of the complex data

function coefficient = circularityMeasure(z)
    % Covariance
    c = mean(abs(z.^2));

    % Pseudo-Covariance
    p = mean(z.^2);

    coefficient = abs(p) / c;
end