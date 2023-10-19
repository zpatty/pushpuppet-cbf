
# qd,dqd,ddqd,d,hp,m,kb,ks,bb,bs,L0,e,Kp,KD,kc,ka,offset,contact)
using Mosek, Convex, StaticArrays, SCS, Plots, DifferentialEquations, MosekTools, Interpolations, MAT

szeros(::Type{T}, N) where T = @SVector zeros(T, N)
szeros(N)= @SVector zeros(N)
szeros(::Type{T}, N1, N2) where T = @SMatrix zeros(T, N1, N2)
szeros(N1, N2)= @SMatrix zeros(N1, N2)

I3 = SMatrix{3,3}(1,0,0,0,1,0,0,0,1)

include("mass.jl")
include("puppet_cl.jl")

d = 0.04;
theta = [0 0.1 0.3 0.5];
phi = [0 0 0 0];
dL = [-0.0 -0.01 -0.07 -0.09];

theta = [0, 1, 1, 1]*0.5;
phi = [0, 0, 0, 0];
dL = [-0.0, -1, -1, -1]*0.088;


dx = theta.*d.*cos.(phi);
dy = theta.*d.*sin.(phi);


td = collect(0:0.01:1)

qd_v = [SA[0.08, 0.0, -0.09] for t in td]
dqd_v = [szeros(3) for t in td];
ddqd_v = [szeros(3) for t in td];

qd = SA[0.08, 0.0, -0.09]
dqd = szeros(3)
ddqd = szeros(3)

# qd_v = [SA[0.1*sin(t*2*pi*0.5); 0.0; -0.1*sin(t*2*pi*0.5)] for t in td]
# dqd_v = [SA[0.1*2*pi*0.5*cos(t*2*pi*0.5), 0.0, -0.1*2*pi*0.5*cos(t*2*pi*0.5)] for t in td];
# ddqd_v = [SA[-0.1*(2*pi*0.5)^2*sin(t*2*pi*0.5), 0.0, 0.1*(2*pi*0.5)^2*sin(t*2*pi*0.5)] for t in td];

L0 = 0.1;

q0 = SA[0; 0; 0];
dq0 = SA[0;0;0];

kb = 1;
ks = 1;
bb = 1;
bs = 1;
e = 0.002;
m = 0.015;
Kp = 5;
KD = diagm(SA[0.1; 0.1; 0.1]);
kpu = 10;
kdu = 5;
hp = 0.005;
t = 0;
qhist = [q0];
tauhist = [];
tvec = [t];
chist = [];
tspan = (0.0, td[end])
p = 0.0
p1 = 100;
tau_hist = Vector{Vector{Float64}}()
Vdot_hist = Vector{Float64}()
t_hist = Vector{Float64}()
solver = SCS.Optimizer
# u0 = [q0; dq0; qd_v[1]; dqd_v[1]];
# prob = ODEProblem(puppet_cl_ref, u0, tspan)
u0 = [q0; dq0];
prob = ODEProblem(puppet_cl, u0, tspan)
sol = DifferentialEquations.solve(prob);

qs = reduce(hcat,[q[1:3] for q in sol.u]);
vs = reduce(hcat,[q[4:6] for q in sol.u]);
# ur = reduce(hcat,[q[7:9] for q in sol.u]);

plot(sol.t,qs')
plot!(td,reduce(hcat,qd_v)')



inds = [length(t_hist) - findfirst(x -> x == t,reverse(t_hist)) for t in sol.t]
tvec = t_hist[inds]
Vdot_p = Vdot_hist[inds]
tau_p = tau_hist[inds]
plot(tvec,Vdot_p)



chist = Vector{Float64}()
for i = 1:length(sol.t)
    del = sqrt(qs[1,i]^2+qs[2,i]^2);
    theta = del/d;
    if theta < 1e-6
        c = (L0+qs[3,i])/3;
    else
        c = 2*((L0 + qs[3,i])/(theta) - d)*sin(theta/6);
    end
    append!(chist, c)
end
plot(sol.t,[chist .- e -qs[3,:]])

path = "/home/zach/Documents/git-repos/Dojo_sim/puppet-cbf/"
file = matopen(path*"cbf_passive.mat", "w")
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




