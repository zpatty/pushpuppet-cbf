function puppet_cl_ref(u,p,t)
    q = u[SA[1,2,3]]
    dq = u[SA[4,5,6]]
    qr = u[SA[7,8,9]];
    dqr = u[SA[10,11,12]];

    dx = q[1];
    dy = q[2];
    dL = q[3];

    I = [1/2*m*d^2; 1/2*m*d^2; 1/12*m*(3*d^2 + hp^2)];
    M = szeros(3,3);
    C = szeros(3,3);
    G = szeros(3);
    dM = szeros(3,3);

    # i = findfirst(x -> x >= t, td)
    # qd = qd_v[i]
    # dqd = dqd_v[i]
    # ddqd = ddqd_v[i]

    for i = 1:3
        Mi, Ci, Gi, dMi = mass(SA[dx/3*i;dy/3*i;dL/3*i],dq,m,I,L0/3*i,d);
        M = M + Mi;
        C = C + Ci;
        G = G + Gi*m*9.81;
        dM = dM + dMi;
    end

    K = diagm(SA[kb;kb;ks]);
    D = diagm(SA[bb;bb;bs]);

    del = sqrt(dx^2+dy^2);
    theta = del/d;

    if theta < 1e-6
        c = (L0+dL)/3;
        dtheta_dx = 0;
        dtheta_dy = 0;
        dc_dx = dtheta_dx;
        dc_dy = dtheta_dy;
        dc_dL = 1/3;
        dJ = SA[0, 0, 0];
    else
        c = 2*((L0 + dL)/(theta) - d)*sin(theta/6);
        dc_dtheta = - (cos(theta/6)*(2*d - (2*(L0 + dL))/theta))/6 - (2*sin(theta/6)*(L0 + dL))/theta^2;
        dtheta_dx = dx/del/d;
        dtheta_dy = dy/del/d;
        dc_dx = dc_dtheta*dtheta_dx;
        dc_dy = dc_dtheta*dtheta_dy;
        dc_dL = 2/theta*sin(theta/6);
        dJ1_dq = SA[(dx^2*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - ((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2))/(d*(dx^2 + dy^2)^(1/2)) - (dx*((d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) (dx*dy*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - (dx*((d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) -(dx*((2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2) - (d*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)^(1/2))))/(d*(dx^2 + dy^2)^(1/2))];
        dJ2_dq = SA[(dx*dy*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - (dy*((d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) (dy^2*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - ((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2))/(d*(dx^2 + dy^2)^(1/2)) - (dy*((d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) -(dy*((2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2) - (d*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)^(1/2))))/(d*(dx^2 + dy^2)^(1/2))];
        dJ3_dq = SA[(dx*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)) - (2*d*dx*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2)^(3/2) (dy*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)) - (2*d*dy*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2)^(3/2) 0];
        # dqi = SA[ddx; ddy; ddL];
        dJ = SA[(dJ1_dq*dq)[1], (dJ2_dq*dq)[1], (dJ3_dq*dq)[1]];
    end

    # if dx < 10e-6 && dy < 10e-6
    #     dJ = SA[0, 0, 0];
    # else
    #     dqi = dq;
    #     qi = q;
    #     dJ = SA[2*dqi[2]*((6*d*qi[1]*qi[2]*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(5/2) - (2*qi[1]*qi[2]*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(3*(qi[1]^2 + qi[2]^2)^2) + (qi[1]*qi[2]*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(6*d*(qi[1]^2 + qi[2]^2)^(3/2)) + (qi[1]*qi[2]*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(36*d^2*(qi[1]^2 + qi[2]^2))) - 2*dqi[1]*((2*d*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(3/2) + (2*qi[1]^2*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(3*(qi[1]^2 + qi[2]^2)^2) + (cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(6*d*(qi[1]^2 + qi[2]^2)^(1/2)) - (qi[1]^2*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(6*d*(qi[1]^2 + qi[2]^2)^(3/2)) - (qi[1]^2*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(36*d^2*(qi[1]^2 + qi[2]^2)) - (6*d*qi[1]^2*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(5/2)), 2*dqi[1]*((6*d*qi[1]*qi[2]*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(5/2) - (2*qi[1]*qi[2]*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(3*(qi[1]^2 + qi[2]^2)^2) + (qi[1]*qi[2]*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(6*d*(qi[1]^2 + qi[2]^2)^(3/2)) + (qi[1]*qi[2]*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(36*d^2*(qi[1]^2 + qi[2]^2))) - 2*dqi[2]*((2*d*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(3/2) + (2*qi[2]^2*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(3*(qi[1]^2 + qi[2]^2)^2) + (cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(6*d*(qi[1]^2 + qi[2]^2)^(1/2)) - (qi[2]^2*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(6*d*(qi[1]^2 + qi[2]^2)^(3/2)) - (qi[2]^2*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(36*d^2*(qi[1]^2 + qi[2]^2)) - (6*d*qi[2]^2*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(5/2)), 0];
    
    # end

    dc_dq = SA[dc_dx, dc_dy, dc_dL];

    ur_nom = (-kpu*(qr - qd) - kdu*dqr);

    # tau = (G + C*dq + M*ddqd + K*qd + D*dqd + Kp*(qd-q) + KD*(dqd - dq));
    tau = (G + C*dqr + M*ur_nom + K*qr + D*dqr + Kp*(qr-q) + KD*(dqr - dq));


    qdd = M\(-C*dq - K*q - D*dq - G + tau);

    J = dc_dq;
    h = c - e;
    h2 = -dL + 0;
    J2 = SA[0,0,-1];
    dJ2 = SA[0,0,0];
    N = length(q);
    w = 1;
    # H = diag([zeros(N,1); ones(3,1)*w]);
    # f = [zeros(N,1); -2*tau];
    
    # Aeq = [M -I3 szeros(3,3)];
    # beq = -C*dq - K*q - D*dq - G;

    Ain = [(dq-dqr)' -(dq-dqr)'; 
            -J'/M szeros(1,N); 
            -J2'/M szeros(1,N)]

    bin = [(dq-dqr)'*(KD + D)*(dq-dqr) + (dq-dqr)'*(tau - M*ur_nom); 
        dJ' * dq + p1^2*h + (2*p1)*J'*dq - J'*(M\(C*dq + K*q + D*dq + G)); 
        dJ2' * dq + p1^2*h2 + (2*p1)*J2'*dq - J2'*(M\(C*dq + K*q + D*dq + G))]

    x = Convex.Variable(2*N)
    problem = minimize(sumsquares(tau - M*ur_nom + M*x[4:6] - x[1:3]) + 10*sumsquares(x[4:6] - ur_nom)) # objective
    problem.constraints += [bin - (Ain*x) >= 0.0]; # constraint
    # problem.constraints += Aeq*x == beq; 
    
    Convex.solve!(problem, solver; silent_solver = true)

    tau_c = x.value[1:3]
    ur = x.value[4:6]
    push!(tau_hist, tau_c)
    push!(Vdot_hist, -(bin - Ain*x.value)[1])
    push!(t_hist, t)
    # if -(dq-dqr)'*(KD + D)*(dq-dqr) +(dq-dqr)'*(tau_c - tau + M*ur_nom - M*ur) > 0.001
    #     @show (bin - Ain*x.value)[1]
    # end

    qdd = M\(-C*dq - K*q - D*dq - G + tau_c)
    return [dq; qdd; dqr; ur]
end


function puppet_cl(u,p,t)
    q = u[SA[1,2,3]]
    dq = u[SA[4,5,6]]

    dx = q[1];
    dy = q[2];
    dL = q[3];

    I = [1/2*m*d^2; 1/2*m*d^2; 1/12*m*(3*d^2 + hp^2)];
    M = szeros(3,3);
    C = szeros(3,3);
    G = szeros(3);
    dM = szeros(3,3);

    i = findfirst(x -> x >= t, td)
    qd = qd_v[i]
    dqd = dqd_v[i]
    ddqd = ddqd_v[i]

    for i = 1:3
        Mi, Ci, Gi, dMi = mass(SA[dx/3*i;dy/3*i;dL/3*i],dq,m,I,L0/3*i,d);
        M = M + Mi;
        C = C + Ci;
        G = G + Gi*m*9.81;
        dM = dM + dMi;
    end

    K = diagm(SA[kb;kb;ks]);
    D = diagm(SA[bb;bb;bs]);

    del = sqrt(dx^2+dy^2);
    theta = del/d;

    if theta < 1e-6
        c = (L0+dL)/3;
        dtheta_dx = 0;
        dtheta_dy = 0;
        dc_dx = dtheta_dx;
        dc_dy = dtheta_dy;
        dc_dL = 1/3;
        dJ = SA[0, 0, 0];
    else
        c = 2*((L0 + dL)/(theta) - d)*sin(theta/6);
        dc_dtheta = - (cos(theta/6)*(2*d - (2*(L0 + dL))/theta))/6 - (2*sin(theta/6)*(L0 + dL))/theta^2;
        dtheta_dx = dx/del/d;
        dtheta_dy = dy/del/d;
        dc_dx = dc_dtheta*dtheta_dx;
        dc_dy = dc_dtheta*dtheta_dy;
        dc_dL = 2/theta*sin(theta/6);
        dJ1_dq = SA[(dx^2*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - ((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2))/(d*(dx^2 + dy^2)^(1/2)) - (dx*((d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) (dx*dy*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - (dx*((d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) -(dx*((2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2) - (d*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)^(1/2))))/(d*(dx^2 + dy^2)^(1/2))];
        dJ2_dq = SA[(dx*dy*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - (dy*((d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) (dy^2*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - ((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2))/(d*(dx^2 + dy^2)^(1/2)) - (dy*((d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) -(dy*((2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2) - (d*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)^(1/2))))/(d*(dx^2 + dy^2)^(1/2))];
        dJ3_dq = SA[(dx*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)) - (2*d*dx*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2)^(3/2) (dy*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)) - (2*d*dy*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2)^(3/2) 0];
        # dqi = SA[ddx; ddy; ddL];
        dJ = SA[(dJ1_dq*dq)[1], (dJ2_dq*dq)[1], (dJ3_dq*dq)[1]];
    end

    # if dx < 10e-6 && dy < 10e-6
    #     dJ = SA[0, 0, 0];
    # else
    #     dqi = dq;
    #     qi = q;
    #     dJ = SA[2*dqi[2]*((6*d*qi[1]*qi[2]*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(5/2) - (2*qi[1]*qi[2]*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(3*(qi[1]^2 + qi[2]^2)^2) + (qi[1]*qi[2]*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(6*d*(qi[1]^2 + qi[2]^2)^(3/2)) + (qi[1]*qi[2]*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(36*d^2*(qi[1]^2 + qi[2]^2))) - 2*dqi[1]*((2*d*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(3/2) + (2*qi[1]^2*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(3*(qi[1]^2 + qi[2]^2)^2) + (cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(6*d*(qi[1]^2 + qi[2]^2)^(1/2)) - (qi[1]^2*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(6*d*(qi[1]^2 + qi[2]^2)^(3/2)) - (qi[1]^2*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(36*d^2*(qi[1]^2 + qi[2]^2)) - (6*d*qi[1]^2*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(5/2)), 2*dqi[1]*((6*d*qi[1]*qi[2]*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(5/2) - (2*qi[1]*qi[2]*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(3*(qi[1]^2 + qi[2]^2)^2) + (qi[1]*qi[2]*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(6*d*(qi[1]^2 + qi[2]^2)^(3/2)) + (qi[1]*qi[2]*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(36*d^2*(qi[1]^2 + qi[2]^2))) - 2*dqi[2]*((2*d*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(3/2) + (2*qi[2]^2*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(3*(qi[1]^2 + qi[2]^2)^2) + (cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(6*d*(qi[1]^2 + qi[2]^2)^(1/2)) - (qi[2]^2*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(6*d*(qi[1]^2 + qi[2]^2)^(3/2)) - (qi[2]^2*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(36*d^2*(qi[1]^2 + qi[2]^2)) - (6*d*qi[2]^2*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(5/2)), 0];
    
    # end

    dc_dq = SA[dc_dx, dc_dy, dc_dL];


    tau = (G + C*dq + M*ddqd + K*qd + D*dqd + Kp*(qd-q) + KD*(dqd - dq));
    # tau = (G + C*dqr + M*ur_nom + K*qr + D*dqr + Kp*(qr-q) + KD*(dqr - dq));


    qdd = M\(-C*dq - K*q - D*dq - G + tau);

    J = dc_dq;
    h = c - e;
    h2 = -dL + 0;
    J2 = SA[0,0,-1];
    dJ2 = SA[0,0,0];
    N = length(q);
    w = 1;
    # H = diag([zeros(N,1); ones(3,1)*w]);
    # f = [zeros(N,1); -2*tau];
    
    # Aeq = [M -I3 szeros(3,3)];
    # beq = -C*dq - K*q - D*dq - G;

    Ain = [-J'/M; 
            -J2'/M]

    bin = [dJ' * dq + p1^2*h + (2*p1)*J'*dq - J'*(M\(C*dq + K*q + D*dq + G)); 
        dJ2' * dq + p1^2*h2 + (2*p1)*J2'*dq - J2'*(M\(C*dq + K*q + D*dq + G))]

    x = Convex.Variable(N)
    problem = minimize(sumsquares(tau - x)) # objective
    problem.constraints += [bin - (Ain*x) >= 0.0]; # constraint
    # problem.constraints += Aeq*x == beq; 
    
    Convex.solve!(problem, solver; silent_solver = true)

    tau_c = vec(x.value)
    push!(tau_hist, tau_c)
    push!(Vdot_hist, (dq-dqd)'*tau_c - (dq-dqd)'*(KD + D)*(dq-dqd) - (dq-dqd)'*(tau))
    push!(t_hist, t)
    # if -(dq-dqr)'*(KD + D)*(dq-dqr) +(dq-dqr)'*(tau_c - tau + M*ur_nom - M*ur) > 0.001
    #     @show (bin - Ain*x.value)[1]
    # end

    qdd = M\(-C*dq - K*q - D*dq - G + tau_c)
    return [dq; qdd]
end



function puppet_calc(u,p,t)
    q = u[SA[1,2,3]]
    dq = u[SA[4,5,6]]
    qr = u[SA[7,8,9]];
    dqr = u[SA[10,11,12]];

    dx = q[1];
    dy = q[2];
    dL = q[3];

    I = [1/2*m*d^2; 1/2*m*d^2; 1/12*m*(3*d^2 + hp^2)];
    M = szeros(3,3);
    C = szeros(3,3);
    G = szeros(3);
    dM = szeros(3,3);

    # i = findfirst(x -> x => t, td)
    # qd = qd_array(:,i)

    for i = 1:3
        Mi, Ci, Gi, dMi = mass(SA[dx/3*i;dy/3*i;dL/3*i],dq,m,I,L0/3*i,d);
        M = M + Mi;
        C = C + Ci;
        G = G + Gi*m*9.81;
        dM = dM + dMi;
    end

    K = diagm(SA[kb;kb;ks]);
    D = diagm(SA[bb;bb;bs]);

    del = sqrt(dx^2+dy^2);
    theta = del/d;

    if theta < 1e-6
        c = (L0+dL)/3;
        dtheta_dx = 0;
        dtheta_dy = 0;
        dc_dx = dtheta_dx;
        dc_dy = dtheta_dy;
        dc_dL = 1/3;
        dJ = SA[0, 0, 0];
    else
        c = 2*((L0 + dL)/(theta) - d)*sin(theta/6);
        dc_dtheta = - (cos(theta/6)*(2*d - (2*(L0 + dL))/theta))/6 - (2*sin(theta/6)*(L0 + dL))/theta^2;
        dtheta_dx = dx/del/d;
        dtheta_dy = dy/del/d;
        dc_dx = dc_dtheta*dtheta_dx;
        dc_dy = dc_dtheta*dtheta_dy;
        dc_dL = 2/theta*sin(theta/6);
        dJ1_dq = SA[(dx^2*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - ((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2))/(d*(dx^2 + dy^2)^(1/2)) - (dx*((d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) (dx*dy*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - (dx*((d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) -(dx*((2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2) - (d*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)^(1/2))))/(d*(dx^2 + dy^2)^(1/2))];
        dJ2_dq = SA[(dx*dy*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - (dy*((d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) (dy^2*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - ((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2))/(d*(dx^2 + dy^2)^(1/2)) - (dy*((d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) -(dy*((2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2) - (d*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)^(1/2))))/(d*(dx^2 + dy^2)^(1/2))];
        dJ3_dq = SA[(dx*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)) - (2*d*dx*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2)^(3/2) (dy*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)) - (2*d*dy*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2)^(3/2) 0];
        # dqi = SA[ddx; ddy; ddL];
        dJ = SA[(dJ1_dq*dq)[1], (dJ2_dq*dq)[1], (dJ3_dq*dq)[1]];
    end

    # if dx < 10e-6 && dy < 10e-6
    #     dJ = SA[0, 0, 0];
    # else
    #     dqi = dq;
    #     qi = q;
    #     dJ = SA[2*dqi[2]*((6*d*qi[1]*qi[2]*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(5/2) - (2*qi[1]*qi[2]*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(3*(qi[1]^2 + qi[2]^2)^2) + (qi[1]*qi[2]*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(6*d*(qi[1]^2 + qi[2]^2)^(3/2)) + (qi[1]*qi[2]*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(36*d^2*(qi[1]^2 + qi[2]^2))) - 2*dqi[1]*((2*d*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(3/2) + (2*qi[1]^2*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(3*(qi[1]^2 + qi[2]^2)^2) + (cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(6*d*(qi[1]^2 + qi[2]^2)^(1/2)) - (qi[1]^2*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(6*d*(qi[1]^2 + qi[2]^2)^(3/2)) - (qi[1]^2*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(36*d^2*(qi[1]^2 + qi[2]^2)) - (6*d*qi[1]^2*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(5/2)), 2*dqi[1]*((6*d*qi[1]*qi[2]*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(5/2) - (2*qi[1]*qi[2]*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(3*(qi[1]^2 + qi[2]^2)^2) + (qi[1]*qi[2]*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(6*d*(qi[1]^2 + qi[2]^2)^(3/2)) + (qi[1]*qi[2]*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(36*d^2*(qi[1]^2 + qi[2]^2))) - 2*dqi[2]*((2*d*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(3/2) + (2*qi[2]^2*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(3*(qi[1]^2 + qi[2]^2)^2) + (cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(6*d*(qi[1]^2 + qi[2]^2)^(1/2)) - (qi[2]^2*cos((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(6*d*(qi[1]^2 + qi[2]^2)^(3/2)) - (qi[2]^2*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(2*d - (2*d*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(1/2)))/(36*d^2*(qi[1]^2 + qi[2]^2)) - (6*d*qi[2]^2*sin((qi[1]^2 + qi[2]^2)^(1/2)/(6*d))*(L0 + dL))/(qi[1]^2 + qi[2]^2)^(5/2)), 0];
    
    # end

    dc_dq = SA[dc_dx, dc_dy, dc_dL];

    ur_nom = (-10*(qr - qd) - 5*dqr);

    # tau = (G + C*dq + M*ddqd + K*qd + D*dqd + Kp*(qd-q) + KD*(dqd - dq));
    tau = (G + C*dqr + M*ur_nom + K*qr + D*dqr + Kp*(qr-q) + KD*(dqr - dq));


    qdd = M\(-C*dq - K*q - D*dq - G + tau);

    J = dc_dq;
    h = c - e;
    h2 = -dL + 0;
    J2 = SA[0,0,-1];
    dJ2 = SA[0,0,0];
    N = length(q);
    w = 1;
    # H = diag([zeros(N,1); ones(3,1)*w]);
    # f = [zeros(N,1); -2*tau];
    
    # Aeq = [M -I3 szeros(3,3)];
    # beq = -C*dq - K*q - D*dq - G;

    Ain = [(dq-dqr)' -(dq-dqr)'; 
            -J'/M szeros(1,N); 
            -J2'/M szeros(1,N)]

    bin = [(dq-dqr)'*(KD + D)*(dq-dqr) + (dq-dqr)'*(tau - M*ur_nom); 
        dJ' * dq + p1^2*h + (2*p1)*J'*dq - J'*(M\(C*dq + K*q + D*dq + G)); 
        dJ2' * dq + p1^2*h2 + (2*p1)*J2'*dq - J2'*(M\(C*dq + K*q + D*dq + G))]

    x = Convex.Variable(2*N)
    problem = minimize(sumsquares(tau - M*ur_nom + M*x[4:6] - x[1:3]) + sumsquares(x[4:6] - ur_nom)) # objective
    problem.constraints += [bin - (Ain*x) >= 0.0]; # constraint
    # problem.constraints += Aeq*x == beq; 
    
    Convex.solve!(problem, solver; silent_solver = true)

    tau_c = x.value[1:3]
    ur = x.value[4:6]

    # push!(t_hist, t)

    qdd = M\(-C*dq - K*q - D*dq - G + tau_c)

    push!(tau_hist, tau_c)
    push!(Vdot_hist, -(dq-dqr)'*(KD + D)*(dq-dqr) +(dq-dqr)'*(tau_c - tau + M*ur_nom - M*ur))


    return [dq; qdd; dqr; ur]
end