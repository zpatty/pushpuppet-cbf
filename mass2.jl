########################################################
# This function calculates the Mass and Coriolis + Gravity Matrices
# given the current state and geometric parameters of the pushpuppet robot
########################################################

function MC(q,dq,m,r,L0,d,N)
    # q = state
    # dq = velocity
    # ddq = acceleration
    # m = mass of a module;
    # r = radius of the hex plate;
    # L0 = Length of a module;
    # d = distance to cable;
    # N = number of links (number of modules plus number of motors)

    
    a_grav = SA[0;0;-9.81;0;0;0];

    I_i = zeros(6,6,N);
    # inertia = SA[1/4*m*r^2; 1/4*m*r^2; 1/2*m*r^2];
    I_i[:,:,1] = diagm(SA[m;m;m;1/4*m*r^2; 1/4*m*r^2; 1/2*m*r^2]);
    I_i[:,:,2] = I_i[:,:,1];


    Smod = zeros(6,3,N);
    Xup = zeros(6,6,N);
    v = zeros(6,N);
    a = zeros(6,N);
    f = zeros(6,N);
    C = zeros(3*N);

    Xtree = SMatrix{6,6}(1I);
    i = 1;
    qi = SA[1,2,3];
    XJ, Smodi,dummy,dJ = PCC_jacobian(q[qi],d,L0,dq[qi]);
    Smod[:,:,i] = Smodi;
    Xup[:,:,1] = XJ*Xtree;
    Xtree = SMatrix{6,6}(1I);
    vJ = Smod[:,:,i]*dq[qi];
    v[:,i] = vJ;
    a[:,i] = Xup[:,:,1]*(-a_grav) + spatial_cross(v[:,i])*vJ + dJ * dq[qi];
    f[:,i] = I_i[:,:,i]*a[:,i] + -spatial_cross(v[:,i])'*I_i[:,:,i]*v[:,i];

    # Recursive Newton Euler to Calculate C+G
    for i=2:N
        qi = (i-1)*3 .+ SA[1,2,3]
        XJ, Smod[:,:,i],dummy,dJ = PCC_jacobian(q[qi],d,L0,dq[qi]);
        Xup[:,:,i] = XJ*Xtree;
        Xtree = SMatrix{6,6}(1I);
        vJ = Smod[:,:,i]*dq[qi];
        v[:,i] = Xup[:,:,i]*v[:,i-1] + vJ;
        a[:,i] = Xup[:,:,i]*a[:,i-1] + dJ * dq[qi] + spatial_cross(v[:,i])*vJ;
        f[:,i] = I_i[:,:,i]*a[:,i] + -spatial_cross(v[:,i])'*I_i[:,:,i]*v[:,i];
    end
    for i = N:-1:1
        qi = (i-1)*3 .+ SA[1,2,3]#(i-1)*3+SA[1,2,3]*i;
        C[qi,1] = 2 * Smod[:,:,i]' * f[:,i];
        if i != 1
            f[:,i-1] = f[:,i-1] + Xup[:,:,i]'*f[:,i];
        end
    end

    # Composite Rigid Body Algorithm to calculate M
    I_iC = I_i;				# composite inertia calculation

    for i = N:-1:1
        if i != 1
            I_iC[:,:,i-1] = I_iC[:,:,i-1] + Xup[:,:,i]'*I_iC[:,:,i]*Xup[:,:,i];
        end
    end

    H = zeros(N*3,N*3);

    for i = 1:N
        qi = (i-1)*3 .+ SA[1,2,3];
        fh3 = I_iC[:,:,i] * Smod[:,:,i];
        H[qi,qi] = Smod[:,:,i]' * fh3;
        j = i;
        while j < N+1 && j > 1
            fh3 = Xup[:,:,j]' * fh3;
            j = j - 1;
            qj = (j-1)*3 .+ SA[1,2,3];
            H[qi,qj] = fh3' * Smod[:,:,j]; #(S{j}' * fh)';
            H[qj,qi] = Smod[:,:,j]' * fh3;
        end
    end

    M = H;
    return (M, C)
end

function PCC_jacobian(q,d,L0,dq) 

    dx = q[1];
    dy = q[2];
    dL = q[3];
    dx_dot = dq[1];
    dy_dot = dq[2];
    dL_dot = dq[3];
    del = sqrt(dx^2 + dy^2);
    I3 = SMatrix{3,3}(1I)
    if del < 10^-6

        T_q = SA[1 0 0 0;
                0 1 0 0;
                0 0 1 L0+dL; 
                0 0 0 1];
        J = SA[(L0 + dL)/(2*d)               0 0
            0 (L0 + dL)/(2*d) 0;
            0               0 1;
            0            -1/d 0;
            1/d               0 0;
            0               0 0];
        dJ = SA[              dL_dot/(2*d)                               0 -dx_dot/(2*d)
            0                    dL_dot/(2*d) -dy_dot/(2*d)
            (dx_dot*(L0 + dL))/(6*d^2) (L0*dy_dot + dL*dy_dot)/(6*d^2)             0
            0                               0             0
            0                               0             0
            dy_dot/(2*d^2)                 -dx_dot/(2*d^2)             0];
    else
    
        Rq = SA[1 + dx^2/del^2*(cos(del/d) - 1) dx*dy/del^2*(cos(del/d) - 1) -dx/del*sin(del/d);
            dx*dy/del^2*(cos(del/d) - 1) 1+dy^2/del^2*(cos(del/d) - 1) -dy/del*sin(del/d);
            dx/del*sin(del/d) dy/del*sin(del/d) cos(del/d)]';
    
        tq = d*(L0 + dL)/del^2*[dx*(1 - cos(del/d)); dy*(1 - cos(del/d)); del*sin(del/d)];
    
        T_q = [Rq tq; 0 0 0 1];
        J = [((dy^2 + dx^2*cos((dx^2 + dy^2)^(1/2)/d))*((dx^2*sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL))/(dx^2 + dy^2)^(3/2) - (d*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2) + (2*d*dx^2*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2)^2))/(dx^2 + dy^2) - (dx^2*sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - d*sin((dx^2 + dy^2)^(1/2)/d)))/(dx^2 + dy^2)^2 + (dx^2*dy^2*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*(sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - 2*d + 2*d*cos((dx^2 + dy^2)^(1/2)/d)))/(dx^2 + dy^2)^3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0 (d*dx*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2);
            0 ((dx^2 + dy^2*cos((dx^2 + dy^2)^(1/2)/d))*((dy^2*sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL))/(dx^2 + dy^2)^(3/2) - (d*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2) + (2*d*dy^2*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2)^2))/(dx^2 + dy^2) - (dy^2*sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - d*sin((dx^2 + dy^2)^(1/2)/d)))/(dx^2 + dy^2)^2 + (dx^2*dy^2*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*(sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - 2*d + 2*d*cos((dx^2 + dy^2)^(1/2)/d)))/(dx^2 + dy^2)^3 (d*dy*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2);
            (dx*(L0 + dL)*((dx^2 + dy^2)^(3/2) - d*sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)))/(dx^2 + dy^2)^(5/2)                                                                                                                                                                                                                                                                                                                                                                                                                                                                            (dy*(L0 + dL)*((dx^2 + dy^2)^(3/2) - d*sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)))/(dx^2 + dy^2)^(5/2)    (d*sin((dx^2 + dy^2)^(1/2)/d))/(dx^2 + dy^2)^(1/2);
            -(dx*dy*((dx^2 + dy^2)^(1/2) - d*sin((dx^2 + dy^2)^(1/2)/d)))/(d*(dx^2 + dy^2)^(3/2))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          -(dy^2*(dx^2 + dy^2)^(1/2) + d*dx^2*sin((dx^2 + dy^2)^(1/2)/d))/(d*(dx^2 + dy^2)^(3/2))                                                     0;
            (dx^2*(dx^2 + dy^2)^(1/2) + d*dy^2*sin((dx^2 + dy^2)^(1/2)/d))/(d*(dx^2 + dy^2)^(3/2))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             (dx*dy*((dx^2 + dy^2)^(1/2) - d*sin((dx^2 + dy^2)^(1/2)/d)))/(d*(dx^2 + dy^2)^(3/2))                                                     0;
            -(dy*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              (dx*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2)                                                     0];
    
        dJ = [((dy^2 + dx^2*cos((dx^2 + dy^2)^(1/2)/d))*((sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*(dx^2 + dy^2)^(3/2)) - (d*dL_dot*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2) + (dL_dot*dx^2*sin((dx^2 + dy^2)^(1/2)/d))/(dx^2 + dy^2)^(3/2) - (5*dx^2*sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*(dx^2 + dy^2)^(5/2)) + (2*dx*dx_dot*sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL))/(dx^2 + dy^2)^(3/2) + (2*d*dL_dot*dx^2*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2)^2 + (d*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*(2*dx*dx_dot + 2*dy*dy_dot))/(dx^2 + dy^2)^2 + (4*d*dx*dx_dot*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2)^2 - (4*d*dx^2*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*(2*dx*dx_dot + 2*dy*dy_dot))/(dx^2 + dy^2)^3 + (dx^2*cos((dx^2 + dy^2)^(1/2)/d)*(L0 + dL)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*d*(dx^2 + dy^2)^2)))/(dx^2 + dy^2) + (((dx^2*sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL))/(dx^2 + dy^2)^(3/2) - (d*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2) + (2*d*dx^2*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2)^2)*(2*dy*dy_dot + 2*dx*dx_dot*cos((dx^2 + dy^2)^(1/2)/d) - (dx^2*sin((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*d*(dx^2 + dy^2)^(1/2))))/(dx^2 + dy^2) - ((dy^2 + dx^2*cos((dx^2 + dy^2)^(1/2)/d))*(2*dx*dx_dot + 2*dy*dy_dot)*((dx^2*sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL))/(dx^2 + dy^2)^(3/2) - (d*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2) + (2*d*dx^2*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2)^2))/(dx^2 + dy^2)^2 - (dL_dot*dx^2*sin((dx^2 + dy^2)^(1/2)/d)*(cos((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - d*sin((dx^2 + dy^2)^(1/2)/d)))/(dx^2 + dy^2)^2 + (dx^2*sin((dx^2 + dy^2)^(1/2)/d)^2*(L0 + dL)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*d*(dx^2 + dy^2)^2) + (2*dx^2*sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - d*sin((dx^2 + dy^2)^(1/2)/d))*(2*dx*dx_dot + 2*dy*dy_dot))/(dx^2 + dy^2)^3 - (dx^2*dy^2*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*((sin((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*(dx^2 + dy^2)^(1/2)) - (cos((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*d)))/(dx^2 + dy^2)^3 - (2*dx*dx_dot*sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - d*sin((dx^2 + dy^2)^(1/2)/d)))/(dx^2 + dy^2)^2 + (dL_dot*dx^2*dy^2*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*(sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - 2*d + 2*d*cos((dx^2 + dy^2)^(1/2)/d)))/(dx^2 + dy^2)^3 - (dx^2*cos((dx^2 + dy^2)^(1/2)/d)*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - d*sin((dx^2 + dy^2)^(1/2)/d))*(2*dx*dx_dot + 2*dy*dy_dot))/(2*d*(dx^2 + dy^2)^(5/2)) - (3*dx^2*dy^2*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*(2*dx*dx_dot + 2*dy*dy_dot)*(sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - 2*d + 2*d*cos((dx^2 + dy^2)^(1/2)/d)))/(dx^2 + dy^2)^4 + (2*dx*dx_dot*dy^2*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*(sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - 2*d + 2*d*cos((dx^2 + dy^2)^(1/2)/d)))/(dx^2 + dy^2)^3 + (2*dx^2*dy*dy_dot*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*(sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - 2*d + 2*d*cos((dx^2 + dy^2)^(1/2)/d)))/(dx^2 + dy^2)^3 - (dx^2*dy^2*sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL)*(2*dx*dx_dot + 2*dy*dy_dot)*(sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - 2*d + 2*d*cos((dx^2 + dy^2)^(1/2)/d)))/(2*d*(dx^2 + dy^2)^(7/2)) 0 (d*dx_dot*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2) - (dx*sin((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*(dx^2 + dy^2)^(3/2)) - (d*dx*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*(2*dx*dx_dot + 2*dy*dy_dot))/(dx^2 + dy^2)^2; 0 ((dx^2 + dy^2*cos((dx^2 + dy^2)^(1/2)/d))*((sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*(dx^2 + dy^2)^(3/2)) - (d*dL_dot*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2) + (dL_dot*dy^2*sin((dx^2 + dy^2)^(1/2)/d))/(dx^2 + dy^2)^(3/2) - (5*dy^2*sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*(dx^2 + dy^2)^(5/2)) + (2*dy*dy_dot*sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL))/(dx^2 + dy^2)^(3/2) + (2*d*dL_dot*dy^2*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2)^2 + (d*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*(2*dx*dx_dot + 2*dy*dy_dot))/(dx^2 + dy^2)^2 + (4*d*dy*dy_dot*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2)^2 - (4*d*dy^2*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*(2*dx*dx_dot + 2*dy*dy_dot))/(dx^2 + dy^2)^3 + (dy^2*cos((dx^2 + dy^2)^(1/2)/d)*(L0 + dL)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*d*(dx^2 + dy^2)^2)))/(dx^2 + dy^2) + (((dy^2*sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL))/(dx^2 + dy^2)^(3/2) - (d*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2) + (2*d*dy^2*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2)^2)*(2*dx*dx_dot + 2*dy*dy_dot*cos((dx^2 + dy^2)^(1/2)/d) - (dy^2*sin((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*d*(dx^2 + dy^2)^(1/2))))/(dx^2 + dy^2) - ((dx^2 + dy^2*cos((dx^2 + dy^2)^(1/2)/d))*(2*dx*dx_dot + 2*dy*dy_dot)*((dy^2*sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL))/(dx^2 + dy^2)^(3/2) - (d*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2) + (2*d*dy^2*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2)^2))/(dx^2 + dy^2)^2 - (dL_dot*dy^2*sin((dx^2 + dy^2)^(1/2)/d)*(cos((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - d*sin((dx^2 + dy^2)^(1/2)/d)))/(dx^2 + dy^2)^2 + (dy^2*sin((dx^2 + dy^2)^(1/2)/d)^2*(L0 + dL)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*d*(dx^2 + dy^2)^2) + (2*dy^2*sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - d*sin((dx^2 + dy^2)^(1/2)/d))*(2*dx*dx_dot + 2*dy*dy_dot))/(dx^2 + dy^2)^3 - (dx^2*dy^2*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*((sin((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*(dx^2 + dy^2)^(1/2)) - (cos((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*d)))/(dx^2 + dy^2)^3 - (2*dy*dy_dot*sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - d*sin((dx^2 + dy^2)^(1/2)/d)))/(dx^2 + dy^2)^2 + (dL_dot*dx^2*dy^2*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*(sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - 2*d + 2*d*cos((dx^2 + dy^2)^(1/2)/d)))/(dx^2 + dy^2)^3 - (dy^2*cos((dx^2 + dy^2)^(1/2)/d)*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - d*sin((dx^2 + dy^2)^(1/2)/d))*(2*dx*dx_dot + 2*dy*dy_dot))/(2*d*(dx^2 + dy^2)^(5/2)) - (3*dx^2*dy^2*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*(2*dx*dx_dot + 2*dy*dy_dot)*(sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - 2*d + 2*d*cos((dx^2 + dy^2)^(1/2)/d)))/(dx^2 + dy^2)^4 + (2*dx*dx_dot*dy^2*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*(sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - 2*d + 2*d*cos((dx^2 + dy^2)^(1/2)/d)))/(dx^2 + dy^2)^3 + (2*dx^2*dy*dy_dot*(L0 + dL)*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*(sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - 2*d + 2*d*cos((dx^2 + dy^2)^(1/2)/d)))/(dx^2 + dy^2)^3 - (dx^2*dy^2*sin((dx^2 + dy^2)^(1/2)/d)*(L0 + dL)*(2*dx*dx_dot + 2*dy*dy_dot)*(sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)^(1/2) - 2*d + 2*d*cos((dx^2 + dy^2)^(1/2)/d)))/(2*d*(dx^2 + dy^2)^(7/2)) (d*dy_dot*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2) - (dy*sin((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*(dx^2 + dy^2)^(3/2)) - (d*dy*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*(2*dx*dx_dot + 2*dy*dy_dot))/(dx^2 + dy^2)^2; (dx_dot*(L0 + dL)*((dx^2 + dy^2)^(3/2) - d*sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)))/(dx^2 + dy^2)^(5/2) - (dx*(L0 + dL)*((cos((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot)*(dx^2 + dy^2)^(1/2))/2 - (3*(2*dx*dx_dot + 2*dy*dy_dot)*(dx^2 + dy^2)^(1/2))/2 + d*sin((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot)))/(dx^2 + dy^2)^(5/2) + (dL_dot*dx*((dx^2 + dy^2)^(3/2) - d*sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)))/(dx^2 + dy^2)^(5/2) - (5*dx*(L0 + dL)*((dx^2 + dy^2)^(3/2) - d*sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2))*(2*dx*dx_dot + 2*dy*dy_dot))/(2*(dx^2 + dy^2)^(7/2)) (dy_dot*(L0 + dL)*((dx^2 + dy^2)^(3/2) - d*sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)))/(dx^2 + dy^2)^(5/2) - (dy*(L0 + dL)*((cos((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot)*(dx^2 + dy^2)^(1/2))/2 - (3*(2*dx*dx_dot + 2*dy*dy_dot)*(dx^2 + dy^2)^(1/2))/2 + d*sin((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot)))/(dx^2 + dy^2)^(5/2) + (dL_dot*dy*((dx^2 + dy^2)^(3/2) - d*sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2)))/(dx^2 + dy^2)^(5/2) - (5*dy*(L0 + dL)*((dx^2 + dy^2)^(3/2) - d*sin((dx^2 + dy^2)^(1/2)/d)*(dx^2 + dy^2))*(2*dx*dx_dot + 2*dy*dy_dot))/(2*(dx^2 + dy^2)^(7/2)) (cos((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*(dx^2 + dy^2)) - (d*sin((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*(dx^2 + dy^2)^(3/2)); (3*dx*dy*((dx^2 + dy^2)^(1/2) - d*sin((dx^2 + dy^2)^(1/2)/d))*(2*dx*dx_dot + 2*dy*dy_dot))/(2*d*(dx^2 + dy^2)^(5/2)) - (dx*dy_dot*((dx^2 + dy^2)^(1/2) - d*sin((dx^2 + dy^2)^(1/2)/d)))/(d*(dx^2 + dy^2)^(3/2)) - (dx_dot*dy*((dx^2 + dy^2)^(1/2) - d*sin((dx^2 + dy^2)^(1/2)/d)))/(d*(dx^2 + dy^2)^(3/2)) - (dx*dy*((2*dx*dx_dot + 2*dy*dy_dot)/(2*(dx^2 + dy^2)^(1/2)) - (cos((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*(dx^2 + dy^2)^(1/2))))/(d*(dx^2 + dy^2)^(3/2)) (3*(dy^2*(dx^2 + dy^2)^(1/2) + d*dx^2*sin((dx^2 + dy^2)^(1/2)/d))*(2*dx*dx_dot + 2*dy*dy_dot))/(2*d*(dx^2 + dy^2)^(5/2)) - ((dy^2*(2*dx*dx_dot + 2*dy*dy_dot))/(2*(dx^2 + dy^2)^(1/2)) + 2*dy*dy_dot*(dx^2 + dy^2)^(1/2) + 2*d*dx*dx_dot*sin((dx^2 + dy^2)^(1/2)/d) + (dx^2*cos((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*(dx^2 + dy^2)^(1/2)))/(d*(dx^2 + dy^2)^(3/2)) 0; ((dx^2*(2*dx*dx_dot + 2*dy*dy_dot))/(2*(dx^2 + dy^2)^(1/2)) + 2*dx*dx_dot*(dx^2 + dy^2)^(1/2) + 2*d*dy*dy_dot*sin((dx^2 + dy^2)^(1/2)/d) + (dy^2*cos((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*(dx^2 + dy^2)^(1/2)))/(d*(dx^2 + dy^2)^(3/2)) - (3*(dx^2*(dx^2 + dy^2)^(1/2) + d*dy^2*sin((dx^2 + dy^2)^(1/2)/d))*(2*dx*dx_dot + 2*dy*dy_dot))/(2*d*(dx^2 + dy^2)^(5/2)) (dx*dy*((2*dx*dx_dot + 2*dy*dy_dot)/(2*(dx^2 + dy^2)^(1/2)) - (cos((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*(dx^2 + dy^2)^(1/2))))/(d*(dx^2 + dy^2)^(3/2)) + (dx*dy_dot*((dx^2 + dy^2)^(1/2) - d*sin((dx^2 + dy^2)^(1/2)/d)))/(d*(dx^2 + dy^2)^(3/2)) + (dx_dot*dy*((dx^2 + dy^2)^(1/2) - d*sin((dx^2 + dy^2)^(1/2)/d)))/(d*(dx^2 + dy^2)^(3/2)) - (3*dx*dy*((dx^2 + dy^2)^(1/2) - d*sin((dx^2 + dy^2)^(1/2)/d))*(2*dx*dx_dot + 2*dy*dy_dot))/(2*d*(dx^2 + dy^2)^(5/2)) 0; (dy*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*(2*dx*dx_dot + 2*dy*dy_dot))/(dx^2 + dy^2)^2 - (dy_dot*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2) + (dy*sin((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*d*(dx^2 + dy^2)^(3/2)) (dx_dot*(cos((dx^2 + dy^2)^(1/2)/d) - 1))/(dx^2 + dy^2) - (dx*(cos((dx^2 + dy^2)^(1/2)/d) - 1)*(2*dx*dx_dot + 2*dy*dy_dot))/(dx^2 + dy^2)^2 - (dx*sin((dx^2 + dy^2)^(1/2)/d)*(2*dx*dx_dot + 2*dy*dy_dot))/(2*d*(dx^2 + dy^2)^(3/2)) 0];
    
    end
    X = adj_calc(T_q,1);
    return (X,J,T_q,dJ)
end
    
function hat(ω)
    return SA[0 -ω[3] ω[2];
            ω[3] 0 -ω[1];
            -ω[2] ω[1] 0]
end

function adj_calc(g,inv)
    # convert transform to rotation and translation
    R = g[SA[1,2,3],SA[1,2,3]];
    p = g[SA[1,2,3],4];
    # get skew symmetric matrix of translation
    p_hat = hat(p);
    if inv == 0
        # package into adjoint
        adj = [R p_hat*R; szeros(3,3) R];
    else
        adj = [R' -R'*p_hat; szeros(3,3) R'];
    end
    return adj
end

function vee(g)
    R = g[SA[1,2,3],SA[1,2,3]];
    v = g[SA[1,2,3],4];
    return SA[v; R[3,2]; R[1,3]; R[2,1]]
end

function spatial_cross(V)

    v = V[SA[1,2,3]];
    w = V[SA[4,5,6]];
    Vx = zeros(6,6);
    skw = hat(w);
    Vx[SA[1,2,3],SA[1,2,3]] = skw;
    Vx[SA[4,5,6],SA[4,5,6]] = skw;
    Vx[SA[1,2,3],SA[4,5,6]] = hat(v);
    return Vx
end