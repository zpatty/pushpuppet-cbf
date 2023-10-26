
# qd,dqd,ddqd,d,hp,m,kb,ks,bb,bs,L0,e,Kp,KD,kc,ka,offset,contact)
using Mosek, Convex, StaticArrays, SCS, Plots, DifferentialEquations, MosekTools, Interpolations, MAT, LinearAlgebra

szeros(::Type{T}, N) where T = @SVector zeros(T, N)
szeros(N)= @SVector zeros(N)
szeros(::Type{T}, N1, N2) where T = @SMatrix zeros(T, N1, N2)
szeros(N1, N2)= @SMatrix zeros(N1, N2)

N = 2
include("mass2.jl")
include("puppet_cl_multi.jl")

d = 0.04;

td = collect(0:0.01:1)

# qd_v = [SA[0.08, 0.0, -0.09] for t in td]
# dqd_v = [szeros(3) for t in td];
# ddqd_v = [szeros(3) for t in td];

qd = SA[0.08, 0.0, -0.05, 0.0, -0.06, -0.1]


# qd_v = [SA[0.1*sin(t*2*pi*0.5); 0.0; -0.1*sin(t*2*pi*0.5)] for t in td]
# dqd_v = [SA[0.1*2*pi*0.5*cos(t*2*pi*0.5), 0.0, -0.1*2*pi*0.5*cos(t*2*pi*0.5)] for t in td];
# ddqd_v = [SA[-0.1*(2*pi*0.5)^2*sin(t*2*pi*0.5), 0.0, 0.1*(2*pi*0.5)^2*sin(t*2*pi*0.5)] for t in td];

L0 = 0.1;


q0 = szeros(3*N);
dq0 = szeros(3*N);
dqd = szeros(3*N)
ddqd = szeros(3*N)

kb = 1;
ks = 10;
bb = 5;
bs = 5;
e = 0.002;
m = 0.015;
r = 0.04;
Kp = 5;
KD = diagm(SA[0.1; 0.1; 0.1; 0.1; 0.1; 0.1]*10);
kpu = 10;
kdu = 5;
hp = 0.005;
t = 0;
qhist = [q0];
tauhist = [];
tvec = [t];
chist = [];
tspan = (0.0, 2)
p = 0.0
p1 = 50;
tau_hist = Vector{Vector{Float64}}()
Vdot_hist = Vector{Float64}()
t_hist = Vector{Float64}()
solver = SCS.Optimizer
# u0 = [q0; dq0; qd; dqd];
# prob = ODEProblem(puppet_cl_ref, u0, tspan)
u0 = [q0; dq0];
prob = ODEProblem(puppet_cl, u0, tspan)
sol = DifferentialEquations.solve(prob);

qs = reduce(hcat,[q[1:6] for q in sol.u]);
vs = reduce(hcat,[q[7:12] for q in sol.u]);
# ur = reduce(hcat,[q[7:9] for q in sol.u]);

plot(sol.t,qs')
plot!(td,reduce(hcat,qd_v)')



inds = [length(t_hist) - findfirst(x -> x == t,reverse(t_hist)) for t in sol.t]
tvec = t_hist[inds]
Vdot_p = Vdot_hist[inds]
tau_p = tau_hist[inds]
plot(tvec,Vdot_p)



chist = zeros(length(sol.t),2)
for i = 1:length(sol.t)
    del = sqrt(qs[1,i]^2+qs[2,i]^2);
    theta = del/d;
    if theta < 1e-6
        c = (L0+qs[3,i])/3;
    else
        c = 2*((L0 + qs[3,i])/(theta) - d)*sin(theta/6);
    end
    del2 = sqrt(qs[4,i]^2+qs[5,i]^2);
    theta2 = del2/d;
    if theta2 < 1e-6
        c2 = (L0+qs[6,i])/3;
    else
        c2 = 2*((L0 + qs[6,i])/(theta2) - d)*sin(theta2/6);
    end
    chist[i,:] = [c- e, c2- e]
end
plot(sol.t,chist)

path = "/home/zach/Documents/git-repos/Dojo_sim/puppet-cbf/"
file = matopen(path*"2_cbf_not_passive.mat", "w")
write(file, "state", qs)
write(file, "velocity", vs)
write(file, "Vdot", Vdot_p)
write(file, "tau", reduce(hcat,tau_p))
write(file, "h", [chist .- e -qs[3,:]])
write(file, "t", sol.t)
write(file, "td", td)
write(file, "qd", reduce(hcat,qd_v))
write(file, "dqd", reduce(hcat,dqd_v))
write(file, "ddqd", reduce(hcat,ddqd_v))
close(file)


# tau_hist = Vector{Vector{Float64}}()
# Vdot_hist = Vector{Float64}()
# t_hist = Vector{Float64}()
# for i=1:length(sol.t)
#     puppet_calc(sol.u[i],p,sol.t[i]);
# end
# plot(sol.t,Vdot_hist)
# plot(sol.t,reduce(hcat,[tau for tau in tau_hist])')




