function tangent = computeTangent(x)
nv = (numel(x) + 1)/4;
ne = nv - 1;
tangent = zeros(ne,3);
for c = 1:ne
    dx = x( 4*c+1:4*c+3 ) - x( 4*(c-1)+1:4*(c-1)+3 );
    tangent(c,:) = dx / norm(dx);
end
end
