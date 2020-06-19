function A_ij = area_covered_detector_circle_distance( c, t, xc, yc, Rc_A, R_d, theta_d )
% calculates A_ij for the current transducer --> defined as the intersection bet A & the circle centered at this
% transducer with radius ct_j 

% Rc_A --> % radius of circle of the ROI (circle A)

pos_det = [R_d*cos(theta_d),R_d*sin(theta_d)];            % position of the detector
distFromDet_v = c*t;                                        % c*t_j for all t_j values
d = sqrt( (pos_det(1)-xc)^2 + (pos_det(2)-yc)^2 );            % distance between the centres of the two circles (circle A & circle centred at the detector)

d1 = (d^2 - distFromDet_v.^2 + Rc_A^2) / (2*d);
d2 = d-d1;

theta1 = 2*acos(d1./Rc_A);
A1 = (Rc_A^2/2)*(theta1-sin(theta1));                       % Area of the first segment of circle

theta2 = 2*acos(d2./distFromDet_v);
A2 = (distFromDet_v.^2/2).*(theta2-sin(theta2));            % Area of the second segment of circle

A_ij = A1+A2;                                             % total Area

%% points for which A_ij equals 0 (i.e. no intersection bet A & the circle centered at this transducer with radius ct_j )
AreaMinLocs = (Rc_A+distFromDet_v)<=d;
A_ij(AreaMinLocs) = 0;

%% points for which A_ij equals the area of the whole ROI circle (A) (i.e. "the circle centered at this transducer with radius ct_j" is big enough to contain all of A)
AreaMaxLocs = abs(Rc_A-distFromDet_v)>=d;
A_ij(AreaMaxLocs) = pi*Rc_A^2;