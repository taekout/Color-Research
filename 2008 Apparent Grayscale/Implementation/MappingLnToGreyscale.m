function Y = MappingLnToGreyscale(Lstar_N)
% input L*N in order to calculate the changed Y value...
% Y = Yn * L*n * (3 / 29)^3 (L*n <=8)
% Y = Yn * ((L8n + 16) / 116)^3 (L*n > 8)

% calculate the number of elements in Lstar_N.
nSize = numel(Lstar_N);
% for performance, allocate the variable Y
Y = ones(size(Lstar_N));

% Convert it to Y of XYZ format.
for i = 1 : nSize
    if(Lstar_N(i) <= 8)
        Y(i) = Lstar_N(i) * (3/29)^3;
    else
        Y(i) = ( (Lstar_N(i) + 16)/116 )^3;
    end
end

end
