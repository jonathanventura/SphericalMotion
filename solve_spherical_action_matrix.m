function E = solve_spherical_action_matrix(u,v)
% u is observations in first frame (3xN)
% v is observations in second frame (3XN)
% E is 3x3x4 matrix of essential matrix solutions

if size(u,1) ~= 3 || size(v,1) ~= 3 || size(u,2) ~= size(v,2),
    error('u and v must be size 3xN');
end

% build matrix of linear constraints
N = size(u,2);
coeffs = zeros(N,6);
for i=1:N,
    coeffs(i,:) = get_row(u(:,i),v(:,i));
end

% get null space
[~,~,V] = svd(coeffs);
B = V(:,4:6);

% build matrix A of non-linear constraints
t2 = B(1,1)^2;
t3 = 2*t2;
t4 = B(2,1)^2;
t5 = 2*t4;
t6 = B(4,1)^2;
t7 = 2*t6;
t8 = t3 + t5 + t7;
t9 = B(3,1)^2;
t10 = B(5,1)^2;
t11 = B(6,1)^2;
t12 = t3 + t5 + t6 + t9 + t10 + t11;
t13 = 4*B(1,1)*B(1,2);
t14 = 4*B(2,1)*B(2,2);
t15 = 2*B(1,1)*B(6,1);
t45 = 2*B(2,1)*B(5,1);
t16 = t15 - t45;
t17 = 2*B(3,1)*B(3,2);
t18 = 2*B(4,1)*B(4,2);
t19 = 2*B(5,1)*B(5,2);
t20 = 2*B(6,1)*B(6,2);
t21 = t13 + t14 + t17 + t18 + t19 + t20;
t22 = B(1,2)^2;
t23 = 2*t22;
t24 = B(2,2)^2;
t25 = 2*t24;
t26 = B(4,2)^2;
t27 = 2*B(1,1)*B(6,2);
t28 = 2*B(1,2)*B(6,1);
t51 = 2*B(2,1)*B(5,2);
t52 = 2*B(2,2)*B(5,1);
t29 = t27 + t28 - t51 - t52;
t30 = 4*B(4,1)*B(4,2);
t31 = t13 + t14 + t30;
t32 = B(1,1)*B(3,2);
t33 = B(1,2)*B(3,1);
t34 = t32 + t33;
t35 = 2*t26;
t36 = t23 + t25 + t35;
t37 = B(3,2)^2;
t38 = B(5,2)^2;
t39 = B(6,2)^2;
t40 = t23 + t25 + t26 + t37 + t38 + t39;
t41 = 2*B(1,2)*B(6,2);
t76 = 2*B(2,2)*B(5,2);
t42 = t41 - t76;
t43 = 4*B(1,1)*B(1,3);
t44 = 4*B(2,1)*B(2,3);
t46 = 2*B(3,1)*B(3,3);
t47 = 2*B(4,1)*B(4,3);
t48 = 2*B(5,1)*B(5,3);
t49 = 2*B(6,1)*B(6,3);
t50 = t43 + t44 + t46 + t47 + t48 + t49;
t53 = 2*B(1,1)*B(6,3);
t54 = 2*B(1,3)*B(6,1);
t82 = 2*B(2,1)*B(5,3);
t83 = 2*B(2,3)*B(5,1);
t55 = t53 + t54 - t82 - t83;
t56 = 4*B(4,1)*B(4,3);
t57 = t43 + t44 + t56;
t58 = 4*B(1,2)*B(1,3);
t59 = 4*B(2,2)*B(2,3);
t60 = B(1,1)*B(3,3);
t61 = B(1,3)*B(3,1);
t62 = t60 + t61;
t63 = 2*B(3,2)*B(3,3);
t64 = 2*B(4,2)*B(4,3);
t65 = 2*B(5,2)*B(5,3);
t66 = 2*B(6,2)*B(6,3);
t67 = t58 + t59 + t63 + t64 + t65 + t66;
t68 = 2*B(1,2)*B(6,3);
t69 = 2*B(1,3)*B(6,2);
t90 = 2*B(2,2)*B(5,3);
t91 = 2*B(2,3)*B(5,2);
t70 = t68 + t69 - t90 - t91;
t71 = 4*B(4,2)*B(4,3);
t72 = t58 + t59 + t71;
t73 = B(1,2)*B(3,3);
t74 = B(1,3)*B(3,2);
t75 = t73 + t74;
t77 = B(1,3)^2;
t78 = 2*t77;
t79 = B(2,3)^2;
t80 = 2*t79;
t81 = B(4,3)^2;
t84 = 2*t81;
t85 = t78 + t80 + t84;
t86 = B(3,3)^2;
t87 = B(5,3)^2;
t88 = B(6,3)^2;
t89 = t78 + t80 + t81 + t86 + t87 + t88;
t92 = 2*B(1,3)*B(6,3);
t94 = 2*B(2,3)*B(5,3);
t93 = t92 - t94;
t95 = 2*t10;
t96 = 2*t11;
t97 = t95 + t96;
t98 = 2*B(1,1)*B(5,1);
t99 = 2*B(2,1)*B(6,1);
t100 = t98 + t99;
t101 = 2*B(1,1)*B(5,2);
t102 = 2*B(1,2)*B(5,1);
t103 = 2*B(2,1)*B(6,2);
t104 = 2*B(2,2)*B(6,1);
t105 = t101 + t102 + t103 + t104;
t106 = 4*B(5,1)*B(5,2);
t107 = 4*B(6,1)*B(6,2);
t108 = t106 + t107;
t109 = 2*t38;
t110 = 2*t39;
t111 = t109 + t110;
t112 = 2*B(1,2)*B(5,2);
t113 = 2*B(2,2)*B(6,2);
t114 = t112 + t113;
t115 = 2*B(1,1)*B(5,3);
t116 = 2*B(1,3)*B(5,1);
t117 = 2*B(2,1)*B(6,3);
t118 = 2*B(2,3)*B(6,1);
t119 = t115 + t116 + t117 + t118;
t120 = 4*B(5,1)*B(5,3);
t121 = 4*B(6,1)*B(6,3);
t122 = t120 + t121;
t123 = 2*B(1,2)*B(5,3);
t124 = 2*B(1,3)*B(5,2);
t125 = 2*B(2,2)*B(6,3);
t126 = 2*B(2,3)*B(6,2);
t127 = t123 + t124 + t125 + t126;
t128 = 4*B(5,2)*B(5,3);
t129 = 4*B(6,2)*B(6,3);
t130 = t128 + t129;
t131 = 2*t87;
t132 = 2*t88;
t133 = t131 + t132;
t134 = 2*B(1,3)*B(5,3);
t135 = 2*B(2,3)*B(6,3);
t136 = t134 + t135;
t137 = B(2,1)*B(3,2);
t138 = B(2,2)*B(3,1);
t139 = t137 + t138;
t140 = B(2,1)*B(3,3);
t141 = B(2,3)*B(3,1);
t142 = t140 + t141;
t143 = B(2,2)*B(3,3);
t144 = B(2,3)*B(3,2);
t145 = t143 + t144;
A = [B(2,1)*t8 - B(2,1)*t12 - B(5,1)*t16 + 2*B(1,1)*B(3,1)*B(4,1), B(2,2)*t8 - B(2,2)*t12 - B(2,1)*t21 + B(2,1)*t31 - B(5,2)*t16 + 2*B(4,1)*t34 - B(5,1)*t29 + 2*B(1,1)*B(3,1)*B(4,2), B(2,2)*t31 - B(2,2)*t21 + B(2,1)*t36 - B(2,1)*t40 + 2*B(4,2)*t34 - B(5,2)*t29 - B(5,1)*t42 + 2*B(1,2)*B(3,2)*B(4,1), B(2,2)*t36 - B(2,2)*t40 - B(5,2)*t42 + 2*B(1,2)*B(3,2)*B(4,2), B(2,3)*t8 - B(2,3)*t12 - B(5,3)*t16 - B(2,1)*t50 + B(2,1)*t57 + 2*B(4,1)*t62 - B(5,1)*t55 + 2*B(1,1)*B(3,1)*B(4,3), B(2,3)*t31 - B(2,3)*t21 - B(2,2)*t50 + 2*B(4,3)*t34 + B(2,2)*t57 - B(5,3)*t29 - B(2,1)*t67 + B(2,1)*t72 + 2*B(4,2)*t62 - B(5,2)*t55 + 2*B(4,1)*t75 - B(5,1)*t70, B(2,3)*t36 - B(2,3)*t40 - B(2,2)*t67 + B(2,2)*t72 - B(5,3)*t42 + 2*B(4,2)*t75 - B(5,2)*t70 + 2*B(1,2)*B(3,2)*B(4,3), B(2,3)*t57 - B(2,3)*t50 + 2*B(4,3)*t62 + B(2,1)*t85 - B(5,3)*t55 - B(2,1)*t89 - B(5,1)*t93 + 2*B(1,3)*B(3,3)*B(4,1), B(2,3)*t72 - B(2,3)*t67 + B(2,2)*t85 - B(2,2)*t89 + 2*B(4,3)*t75 - B(5,3)*t70 - B(5,2)*t93 + 2*B(1,3)*B(3,3)*B(4,2), B(2,3)*t85 - B(2,3)*t89 - B(5,3)*t93 + 2*B(1,3)*B(3,3)*B(4,3);
B(1,1)*t100 - B(5,1)*t12 - B(2,1)*t16 + B(5,1)*t97, B(1,2)*t100 - B(2,1)*t29 - B(5,2)*t12 - B(5,1)*t21 - B(2,2)*t16 + B(1,1)*t105 + B(5,2)*t97 + B(5,1)*t108, B(1,2)*t105 - B(2,1)*t42 - B(5,2)*t21 - B(5,1)*t40 - B(2,2)*t29 + B(1,1)*t114 + B(5,2)*t108 + B(5,1)*t111, B(1,2)*t114 - B(5,2)*t40 - B(2,2)*t42 + B(5,2)*t111, B(1,3)*t100 - B(5,3)*t12 - B(2,1)*t55 - B(5,1)*t50 - B(2,3)*t16 + B(1,1)*t119 + B(5,3)*t97 + B(5,1)*t122, B(1,3)*t105 - B(5,3)*t21 - B(2,2)*t55 - B(2,1)*t70 - B(5,2)*t50 - B(2,3)*t29 - B(5,1)*t67 + B(1,2)*t119 + B(1,1)*t127 + B(5,3)*t108 + B(5,2)*t122 + B(5,1)*t130, B(1,3)*t114 - B(2,2)*t70 - B(5,3)*t40 - B(5,2)*t67 - B(2,3)*t42 + B(1,2)*t127 + B(5,3)*t111 + B(5,2)*t130, B(1,3)*t119 - B(5,3)*t50 - B(2,1)*t93 - B(2,3)*t55 - B(5,1)*t89 + B(1,1)*t136 + B(5,3)*t122 + B(5,1)*t133, B(1,3)*t127 - B(2,2)*t93 - B(5,3)*t67 - B(2,3)*t70 - B(5,2)*t89 + B(1,2)*t136 + B(5,3)*t130 + B(5,2)*t133, B(1,3)*t136 - B(5,3)*t89 - B(2,3)*t93 + B(5,3)*t133;
B(1,1)*t12 - B(1,1)*t8 - B(6,1)*t16 + 2*B(2,1)*B(3,1)*B(4,1), B(1,2)*t12 - B(1,2)*t8 + B(1,1)*t21 - B(1,1)*t31 - B(6,2)*t16 - B(6,1)*t29 + 2*B(4,1)*t139 + 2*B(2,1)*B(3,1)*B(4,2), B(1,2)*t21 - B(1,2)*t31 - B(1,1)*t36 + B(1,1)*t40 - B(6,2)*t29 - B(6,1)*t42 + 2*B(4,2)*t139 + 2*B(2,2)*B(3,2)*B(4,1), B(1,2)*t40 - B(1,2)*t36 - B(6,2)*t42 + 2*B(2,2)*B(3,2)*B(4,2), B(1,3)*t12 - B(1,3)*t8 + B(1,1)*t50 - B(1,1)*t57 - B(6,3)*t16 - B(6,1)*t55 + 2*B(4,1)*t142 + 2*B(2,1)*B(3,1)*B(4,3), B(1,3)*t21 - B(1,3)*t31 + B(1,2)*t50 - B(1,2)*t57 + B(1,1)*t67 - B(1,1)*t72 - B(6,3)*t29 - B(6,2)*t55 - B(6,1)*t70 + 2*B(4,3)*t139 + 2*B(4,2)*t142 + 2*B(4,1)*t145, B(1,3)*t40 - B(1,3)*t36 + B(1,2)*t67 - B(1,2)*t72 - B(6,3)*t42 - B(6,2)*t70 + 2*B(4,2)*t145 + 2*B(2,2)*B(3,2)*B(4,3), B(1,3)*t50 - B(1,3)*t57 - B(1,1)*t85 + B(1,1)*t89 - B(6,3)*t55 - B(6,1)*t93 + 2*B(4,3)*t142 + 2*B(2,3)*B(3,3)*B(4,1), B(1,3)*t67 - B(1,3)*t72 - B(1,2)*t85 + B(1,2)*t89 - B(6,3)*t70 - B(6,2)*t93 + 2*B(4,3)*t145 + 2*B(2,3)*B(3,3)*B(4,2), B(1,3)*t89 - B(1,3)*t85 - B(6,3)*t93 + 2*B(2,3)*B(3,3)*B(4,3);
B(1,1)*t16 - B(6,1)*t12 + B(2,1)*t100 + B(6,1)*t97, B(1,2)*t16 + B(1,1)*t29 - B(6,2)*t12 - B(6,1)*t21 + B(2,2)*t100 + B(2,1)*t105 + B(6,2)*t97 + B(6,1)*t108, B(1,2)*t29 + B(1,1)*t42 - B(6,2)*t21 - B(6,1)*t40 + B(2,2)*t105 + B(2,1)*t114 + B(6,2)*t108 + B(6,1)*t111, B(1,2)*t42 - B(6,2)*t40 + B(2,2)*t114 + B(6,2)*t111, B(1,3)*t16 + B(1,1)*t55 - B(6,3)*t12 - B(6,1)*t50 + B(2,3)*t100 + B(2,1)*t119 + B(6,3)*t97 + B(6,1)*t122, B(1,3)*t29 + B(1,2)*t55 + B(1,1)*t70 - B(6,3)*t21 - B(6,2)*t50 + B(2,3)*t105 - B(6,1)*t67 + B(2,2)*t119 + B(2,1)*t127 + B(6,3)*t108 + B(6,2)*t122 + B(6,1)*t130, B(1,3)*t42 + B(1,2)*t70 - B(6,3)*t40 - B(6,2)*t67 + B(2,3)*t114 + B(2,2)*t127 + B(6,3)*t111 + B(6,2)*t130, B(1,3)*t55 + B(1,1)*t93 - B(6,3)*t50 + B(2,3)*t119 - B(6,1)*t89 + B(2,1)*t136 + B(6,3)*t122 + B(6,1)*t133, B(1,3)*t70 + B(1,2)*t93 - B(6,3)*t67 + B(2,3)*t127 - B(6,2)*t89 + B(2,2)*t136 + B(6,3)*t130 + B(6,2)*t133, B(1,3)*t93 - B(6,3)*t89 + B(2,3)*t136 + B(6,3)*t133;
B(4,1)*t8 + 2*B(4,1)*t9 - B(4,1)*t12, B(4,2)*t8 + 2*B(4,2)*t9 - B(4,2)*t12 - B(4,1)*t21 + B(4,1)*t31 + 4*B(3,1)*B(3,2)*B(4,1), B(4,2)*t31 - B(4,2)*t21 + B(4,1)*t36 + 2*B(4,1)*t37 - B(4,1)*t40 + 4*B(3,1)*B(3,2)*B(4,2), B(4,2)*t36 + 2*B(4,2)*t37 - B(4,2)*t40, B(4,3)*t8 + 2*B(4,3)*t9 - B(4,3)*t12 - B(4,1)*t50 + B(4,1)*t57 + 4*B(3,1)*B(3,3)*B(4,1), B(4,3)*t31 - B(4,3)*t21 - B(4,2)*t50 + B(4,2)*t57 - B(4,1)*t67 + B(4,1)*t72 + 4*B(3,1)*B(3,2)*B(4,3) + 4*B(3,1)*B(3,3)*B(4,2) + 4*B(3,2)*B(3,3)*B(4,1), B(4,3)*t36 + 2*B(4,3)*t37 - B(4,3)*t40 - B(4,2)*t67 + B(4,2)*t72 + 4*B(3,2)*B(3,3)*B(4,2), B(4,3)*t57 - B(4,3)*t50 + B(4,1)*t85 + 2*B(4,1)*t86 - B(4,1)*t89 + 4*B(3,1)*B(3,3)*B(4,3), B(4,3)*t72 - B(4,3)*t67 + B(4,2)*t85 + 2*B(4,2)*t86 - B(4,2)*t89 + 4*B(3,2)*B(3,3)*B(4,3), B(4,3)*t85 + 2*B(4,3)*t86 - B(4,3)*t89;
B(3,1)*t100 - B(4,1)*t16, B(3,2)*t100 - B(4,1)*t29 - B(4,2)*t16 + B(3,1)*t105, B(3,2)*t105 - B(4,1)*t42 - B(4,2)*t29 + B(3,1)*t114, B(3,2)*t114 - B(4,2)*t42, B(3,3)*t100 - B(4,1)*t55 - B(4,3)*t16 + B(3,1)*t119, B(3,3)*t105 - B(4,2)*t55 - B(4,1)*t70 - B(4,3)*t29 + B(3,2)*t119 + B(3,1)*t127, B(3,3)*t114 - B(4,2)*t70 - B(4,3)*t42 + B(3,2)*t127, B(3,3)*t119 - B(4,1)*t93 - B(4,3)*t55 + B(3,1)*t136, B(3,3)*t127 - B(4,2)*t93 - B(4,3)*t70 + B(3,2)*t136, B(3,3)*t136 - B(4,3)*t93];

% gaussian elimination
G = rref(A);

% build action matrix on x
Ax = [ -G([3 5 6],7:10); 0 1 0 0 ];

% eigenvalue decomposition of action matrix
[V,~] = eig(Ax);
xsolns = V(2,:)./V(4,:);
ysolns = V(3,:)./V(4,:);

% extract essential matrix solutions
E = zeros(3,3,4);
for i=1:4,
    xysoln = [xsolns(i);ysolns(i);1];
    e = B*xysoln;
    myE = [e(1) e(2) e(3) ; e(2) -e(1) e(4) ; e(5) e(6) 0 ];
    E(:,:,i) = myE./norm(myE(:));
end

function row = get_row(u,v)
row = [ u(1)*v(1) - u(2)*v(2), u(1)*v(2) + u(2)*v(1), u(3)*v(1), u(3)*v(2), u(1)*v(3), u(2)*v(3) ];