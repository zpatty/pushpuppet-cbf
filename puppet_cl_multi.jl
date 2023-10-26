function puppet_cl_ref(u,p,t)
    q = u[SVector{3*N}(collect(1:3*N))]
    dq = u[SVector{3*N}(collect(3*N+1:2*3*N))]
    qr = u[SVector{3*N}(collect(2*3*N+1:3*3*N))]
    dqr = u[SVector{3*N}(collect(3*3*N+1:4*3*N))]
    inds = @SVector [SVector{3}(collect((i-1)*3 + 1:i*3)) for i=1:N];

    M,C = MC(q,dq,m,r,L0,d,N)
    
    K = diagm(SVector{N*3}(repeat([kb;kb;ks],N)));
    D = diagm(SVector{N*3}(repeat([bb;bb;bs],N)));
    c = zeros(N)
    J = zeros(N,N*3)
    dJ = zeros(N,N*3)
    for i  = 1:N
        dx = q[inds[i][1]];
        dy = q[inds[i][2]];
        dL = q[inds[i][3]];
        del = sqrt(dx^2+dy^2);
        theta = del/d;

        if theta < 1e-6
            c[i] = (L0+dL)/3;
            dtheta_dx = 0;
            dtheta_dy = 0;
            dc_dx = dtheta_dx;
            dc_dy = dtheta_dy;
            dc_dL = 1/3;
            J[i,inds[i]] = SA[dc_dx, dc_dy, dc_dL];
            dJ[i,inds[i]] = SA[0, 0, 0];
        else
            c[i] = 2*((L0 + dL)/(theta) - d)*sin(theta/6);
            dc_dtheta = - (cos(theta/6)*(2*d - (2*(L0 + dL))/theta))/6 - (2*sin(theta/6)*(L0 + dL))/theta^2;
            dtheta_dx = dx/del/d;
            dtheta_dy = dy/del/d;
            dc_dx = dc_dtheta*dtheta_dx;
            dc_dy = dc_dtheta*dtheta_dy;
            dc_dL = 2/theta*sin(theta/6);
            J[i,inds[i]] = SA[dc_dx, dc_dy, dc_dL];
            dJ1_dq = SA[(dx^2*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - ((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2))/(d*(dx^2 + dy^2)^(1/2)) - (dx*((d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) (dx*dy*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - (dx*((d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) -(dx*((2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2) - (d*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)^(1/2))))/(d*(dx^2 + dy^2)^(1/2))];
            dJ2_dq = SA[(dx*dy*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - (dy*((d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) (dy^2*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - ((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2))/(d*(dx^2 + dy^2)^(1/2)) - (dy*((d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) -(dy*((2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2) - (d*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)^(1/2))))/(d*(dx^2 + dy^2)^(1/2))];
            dJ3_dq = SA[(dx*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)) - (2*d*dx*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2)^(3/2) (dy*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)) - (2*d*dy*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2)^(3/2) 0];
            dJ[i,inds[i]] = SA[(dJ1_dq*dq[inds[i]])[1], (dJ2_dq*dq[inds[i]])[1], (dJ3_dq*dq[inds[i]])[1]];
        end
    end


    ur_nom = (-kpu*(qr - qd) - kdu*dqr);

    # tau = (G + C*dq + M*ddqd + K*qd + D*dqd + Kp*(qd-q) + KD*(dqd - dq));
    tau = (C + M*ur_nom + K*qr + D*dqr + Kp*(qr-q) + KD*(dqr - dq));



    h = c .- e;
    h2 = -[q[3]; q[6]];
    J2 = SA[0 0 -1 0 0 0;
            0 0 0 0 0 -1];
    dJ2 = szeros(2,6);
    n = length(q);
    w = 1;
    # H = diag([zeros(N,1); ones(3,1)*w]);
    # f = [zeros(N,1); -2*tau];
    
    # Aeq = [M -I3 szeros(3,3)];
    # beq = -C*dq - K*q - D*dq - G;
    Ain = [(dq-dqr)' -(dq-dqr)'; 
            -J/M szeros(2,n); 
            -J2/M szeros(2,n)]

    bin = [(dq-dqr)'*(KD + D)*(dq-dqr) + (dq-dqr)'*(tau - M*ur_nom); 
        dJ * dq + p1^2*h + (2*p1)*J*dq - J*(M\(C + K*q + D*dq)); 
        dJ2 * dq + p1^2*h2 + (2*p1)*J2*dq - J2*(M\(C + K*q + D*dq))]

    x = Convex.Variable(2*n)
    problem = minimize(sumsquares(tau - M*ur_nom + M*x[SVector{n}(collect(n+1:2*n))] - x[SVector{n}(collect(1:n))]) + sumsquares(x[SVector{n}(collect(n+1:2*n))] - ur_nom)) # objective
    problem.constraints += [bin - (Ain*x) >= 0.0]; # constraint
    
    Convex.solve!(problem, solver; silent_solver = true)

    tau_c = x.value[SVector{n}(collect(1:n))]
    ur = x.value[SVector{n}(collect(n+1:2*n))]
    push!(tau_hist, tau_c)
    push!(Vdot_hist, -(bin - Ain*x.value)[1])
    push!(t_hist, t)
    # if -(dq-dqr)'*(KD + D)*(dq-dqr) +(dq-dqr)'*(tau_c - tau + M*ur_nom - M*ur) > 0.001
    #     @show (bin - Ain*x.value)[1]
    # end

    qdd = M\(-C - K*q - D*dq + tau_c)
    return [dq; qdd; dqr; ur]
end


function puppet_cl(u,p,t)
    q = u[SVector{3*N}(collect(1:3*N))]
    dq = u[SVector{3*N}(collect(3*N+1:2*3*N))]
    inds = @SVector [SVector{3}(collect((i-1)*3 + 1:i*3)) for i=1:N];

    M,C = MC(q,dq,m,r,L0,d,N)
    
    K = diagm(SVector{N*3}(repeat([kb;kb;ks],N)));
    D = diagm(SVector{N*3}(repeat([bb;bb;bs],N)));
    c = zeros(N)
    J = zeros(N,N*3)
    dJ = zeros(N,N*3)
    for i  = 1:N
        dx = q[inds[i][1]];
        dy = q[inds[i][2]];
        dL = q[inds[i][3]];
        del = sqrt(dx^2+dy^2);
        theta = del/d;

        if theta < 1e-6
            c[i] = (L0+dL)/3;
            dtheta_dx = 0;
            dtheta_dy = 0;
            dc_dx = dtheta_dx;
            dc_dy = dtheta_dy;
            dc_dL = 1/3;
            J[i,inds[i]] = SA[dc_dx, dc_dy, dc_dL];
            dJ[i,inds[i]] = SA[0, 0, 0];
        else
            c[i] = 2*((L0 + dL)/(theta) - d)*sin(theta/6);
            dc_dtheta = - (cos(theta/6)*(2*d - (2*(L0 + dL))/theta))/6 - (2*sin(theta/6)*(L0 + dL))/theta^2;
            dtheta_dx = dx/del/d;
            dtheta_dy = dy/del/d;
            dc_dx = dc_dtheta*dtheta_dx;
            dc_dy = dc_dtheta*dtheta_dy;
            dc_dL = 2/theta*sin(theta/6);
            J[i,inds[i]] = SA[dc_dx, dc_dy, dc_dL];
            dJ1_dq = SA[(dx^2*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - ((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2))/(d*(dx^2 + dy^2)^(1/2)) - (dx*((d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) (dx*dy*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - (dx*((d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) -(dx*((2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2) - (d*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)^(1/2))))/(d*(dx^2 + dy^2)^(1/2))];
            dJ2_dq = SA[(dx*dy*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - (dy*((d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dx*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dx*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) (dy^2*((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)))/(d*(dx^2 + dy^2)^(3/2)) - ((cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/6 + (2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2))/(d*(dx^2 + dy^2)^(1/2)) - (dy*((d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(2*L0 + 2*dL))/(6*(dx^2 + dy^2)^(3/2)) + (d*dy*cos((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(3*(dx^2 + dy^2)^(3/2)) - (dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(2*d - (d*(2*L0 + 2*dL))/(dx^2 + dy^2)^(1/2)))/(36*d*(dx^2 + dy^2)^(1/2)) - (4*d^2*dy*sin((dx^2 + dy^2)^(1/2)/(6*d))*(L0 + dL))/(dx^2 + dy^2)^2))/(d*(dx^2 + dy^2)^(1/2)) -(dy*((2*d^2*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2) - (d*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)^(1/2))))/(d*(dx^2 + dy^2)^(1/2))];
            dJ3_dq = SA[(dx*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)) - (2*d*dx*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2)^(3/2) (dy*cos((dx^2 + dy^2)^(1/2)/(6*d)))/(3*(dx^2 + dy^2)) - (2*d*dy*sin((dx^2 + dy^2)^(1/2)/(6*d)))/(dx^2 + dy^2)^(3/2) 0];
            dJ[i,inds[i]] = SA[(dJ1_dq*dq[inds[i]])[1], (dJ2_dq*dq[inds[i]])[1], (dJ3_dq*dq[inds[i]])[1]];
        end
    end



    tau = (C + M*ddqd + K*qd + D*dqd + Kp*(qd-q) + KD*(dqd - dq));

    h = c .- e;
    h2 = -[q[3]; q[6]];
    J2 = SA[0 0 -1 0 0 0;
            0 0 0 0 0 -1];
    dJ2 = szeros(2,6);
    n = length(q);
    w = 1;
    # H = diag([zeros(N,1); ones(3,1)*w]);
    # f = [zeros(N,1); -2*tau];
    
    # Aeq = [M -I3 szeros(3,3)];
    # beq = -C*dq - K*q - D*dq - G;
    Ain = [-J/M; 
            -J2/M]

    bin = [dJ * dq + p1^2*h + (2*p1)*J*dq - J*(M\(C + K*q + D*dq)); 
        dJ2 * dq + p1^2*h2 + (2*p1)*J2*dq - J2*(M\(C + K*q + D*dq))]

    x = Convex.Variable(n)
    problem = minimize(sumsquares(tau - x[SVector{n}(collect(1:n))])) # objective
    problem.constraints += [bin - (Ain*x) >= 0.0]; # constraint
    
    Convex.solve!(problem, solver; silent_solver = true)

    tau_c = x.value[SVector{n}(collect(1:n))]
    push!(tau_hist, tau_c)
    # push!(Vdot_hist, -(bin - Ain*x.value)[1])
    push!(t_hist, t)
    # if -(dq-dqr)'*(KD + D)*(dq-dqr) +(dq-dqr)'*(tau_c - tau + M*ur_nom - M*ur) > 0.001
    #     @show (bin - Ain*x.value)[1]
    # end

    qdd = M\(-C - K*q - D*dq + tau_c)
    return [dq; qdd]
end