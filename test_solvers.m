clear;
clc;
close all;

%% make random 3D points
X = rand(3,3)*2-1;

%% make random rotation and translation
R = rand(3)*2-1;
[U,S,V] = svd(R);
R = U*V.';
if det(R)<0,
    R = -R;
end
% outward-facing
t = R(:,3)-[0;0;1];

%% get points in second camera
PX = bsxfun(@plus,R*X,t);

%% get true E
t_x = [ 0 -t(3) t(2) ;
      t(3) 0 -t(1) ;
     -t(2) t(1) 0 ];
E = t_x*R;
E = E./norm(E(:));

%% get projections in cameras
u = bsxfun(@rdivide,X,X(3,:));
v = bsxfun(@rdivide,PX,PX(3,:));

%% get solutions
Esolns_AM = solve_spherical_action_matrix(u,v);
Esolns_poly = solve_spherical_polynomial(u,v);

%% get error in frobenius norm
error_AM = zeros(1,4);
error_poly = zeros(1,4);
for i=1:4,
    err1 = E-Esolns_AM(:,:,i);
    err2 = E+Esolns_AM(:,:,i);
    error_AM(i) = min( norm(err1(:)), norm(err2(:)) );
    
    err1 = E-Esolns_poly(:,:,i);
    err2 = E+Esolns_poly(:,:,i);
    error_poly(i) = min( norm(err1(:)), norm(err2(:)) );
end

%% display errors
fprintf('action matrix error: ');
fprintf('%0.2f ',error_AM);
fprintf('\n');
fprintf('polynomial error: ');
fprintf('%0.2f ',error_poly);
fprintf('\n');