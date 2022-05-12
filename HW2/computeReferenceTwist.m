function refTwist = computeReferenceTwist(u1, u2, t1, t2)
% refTwist is the reference twist to move u1 with tangent t1 to u2 with
% tangent t2.
ut = parallel_transport(u1, t1, t2);

w = cross(ut,u2);
angle = atan2( norm(w), dot(ut,u2) );
if (dot(t2,w) < 0) 
    angle = -angle;
end

refTwist = angle;
end
