#--------------------------------------------------
#--- Variable Selection for High-dimensional Nodal Attributes 
#--- in Social Networks with Degree Heterogeneity
#--------------------------------------------------
#--------- Simulation Study
#--------- BSM-net: written in Julia 1.3.1
#--------------------------------------------------

iter  = ARGS[1]; # get the number of replication
iter = parse(Int, iter);
# Pkg.status(); # you may need to download packages first by using "Pkg.add"
using Distributions, Roots, StatsBase, CPUTime, Clustering, Random, LinearAlgebra, DelimitedFiles # call the packages
Random.seed!(2021*iter); # set seed

#--------------------------------------------------
#----------------- set up the simulation setting
#--------------------------------------------------
Xtype = 1; # type of X {1,2}
corr = 0.75; # correlation {0.25, 0.75}
sp = 1; # value of s for sparsity {0, 1}
K = 10; # preliminary guess of the model size
code = "jia_example";
#--------------------------------------------------
#----------------- set-up 
#--------------------------------------------------
d = 1000; # dimension
n = 50; # network size
n_a = 30; # number of active nodes
M = 5000; # number of iterations in mcmc
option = [1, 2 , 5, 8]; # the true model
signal = [2.1, -2.4, 2.7, -3.0]; # the true signals

#--------------------------------------------------
#----------------- functions
#--------------------------------------------------
# a function to get column means for a matrix
function colmean(x)
    y = ones(size(x)[2]);
    for i in 1: size(x)[2]
        y[i] = mean(x[:,i]);
    end
    return y
end

# a function to get the row with highest freq in a matrix
function MAP(mat) 
    row_list = [mat[1,:]];
    freq = [1];
    for i in 2:size(mat)[1]
        if mat[i,:] in row_list
            ind = findall(x -> x==mat[i,:], row_list)[1];
            freq[ind] = freq[ind]+1;
        else
            row_list = append!(row_list, [mat[i,:]]);
            freq = append!(freq, [1]);
        end
    end
    freq_max = maximum(freq);
    max_indx = findall(x->x==freq_max, freq);
    return findall(x->x>0, row_list[max_indx][1])
end

function g(x,y)
    return x.*y
end
#--------------------------------------------------
#----------------- read X, Y, delta
#--------------------------------------------------
X = readdlm(join(["/simulation_data/$code","_X_K=$K", "_n=$n", "_d=$d","_Xtype=$Xtype","_corr=$corr","_sp=$sp","_loop=$iter",".csv"]));
Y = readdlm(join(["/simulation_data/$code","_Y_K=$K", "_n=$n", "_d=$d","_Xtype=$Xtype","_corr=$corr","_sp=$sp","_loop=$iter",".csv"]));
delta_true = readdlm(join(["/simulation_data/$code","_delta_K=$K", "_n=$n", "_d=$d","_Xtype=$Xtype","_corr=$corr","_sp=$sp","_loop=$iter",".csv"]));
delta_true = convert(Array{Float64,1},  delta_true[1:(size(delta_true)[2]-1)]);

## get the design matrix T from X
T = fill(0.0, n+d);
for i = 2:n
    for j = 1:(i-1)
        b = fill(0.0, n);
        b[i] = 1;
        b[j] = 1;
        a = vcat(g(X[i,:], X[j,:]), b);
        global T = hcat(T,a);
    end
end
T = T[:,2:size(T)[2]]';

#--------------------------------------------------
#------------------ set up the parameters for the MCMC 
#--------------------------------------------------
f(x) = K/d-x-1.29*sqrt(x*(1-x)/d);
p = fzero(f,[0,1]); # q_1: preliminary inclusion probability
sigma_1 = 0.5; # sd for spike prior (beta)
sigma_2 = sqrt(exp(2.1*n*log(n)^(0.01)/100)/n/(n-1)*2)*200; # sd for spike prior (delta)
sigma_d = log(n)*sqrt(n)/20; # sd for slab prior (beta)
v = 7.3;
w = pi*sqrt((v-2)/3/v);
# initialization
spot = 1;
m = 0;
psi_path = [zeros(Float64, M+1, d+n), ones(Float64,Int.(n*(n-1)/2)), zeros(Int64, M+1, d), ones(Float64,Int.(n*(n-1)/2))];

#--------------------------------------------------
#------------------ start the iteration to update :  
#------------------      beta, delta, z (g in the manuscript), psi, gamma
#--------------------------------------------------

@CPUtime for t in 1:M # record the CPU time

    #--------------------------------------------------
    #------------------ update beta and delta
    #--------------------------------------------------
    Active = findall(x -> x == 1, psi_path[3][t,:]);
    Inactive = findall(x -> x == 0, psi_path[3][t,:]);
    D = Diagonal(1 ./vcat(fill(sigma_2^2, length(Active)), fill(sigma_d^2, n)));
    V = D+T[:,vcat(Active,(d+1):(d+n))]'*
    Diagonal(1 ./psi_path[4])*T[:,vcat(Active,(d+1):(d+n))];
    Sigma = inv(V);
    m = Sigma*T[:,vcat(Active,(d+1):(d+n))]'*
    Diagonal(1 ./psi_path[4])*psi_path[2];
    Sigma = (Sigma+Sigma')/2;
    psi_path[1][t+1,vcat(Active,(d+1):(d+n))] = (Sigma^(0.5)*rand(Normal(0,1), length(vcat(Active,(d+1):(d+n))))+m)';

    if length(Inactive)>0
        psi_path[1][t+1,Inactive] = rand(Normal(0,sigma_1), length(Inactive))';
    end

    #--------------------------------------------------
    #------------------ update z and psi
    #--------------------------------------------------
    m = T[:,vcat(Active,(d+1):(d+n))]*psi_path[1][t+1,vcat(Active,(d+1):(d+n))];
    global spot = 1;
    for i = 2:n
        for j = 1:(i-1)
            if Y[i,j]==1
                psi_path[2][spot] = rand(truncated(Normal(m[spot],sqrt(psi_path[4][spot])) ,0, Inf), 1)[1];
            else
                psi_path[2][spot] = rand(truncated(Normal(m[spot],sqrt(psi_path[4][spot])),-Inf,0), 1)[1];
            end
            psi_path[4][spot] = rand(InverseGamma(v/2+1/2, w^2*v/2+(psi_path[2][spot]-m[spot])^2/2) ,1)[1];
            global spot = spot +1;
        end
    end

    #--------------------------------------------------
    #------------------ update gamma
    #--------------------------------------------------
    for k in 1:d
        p1 = log(p/sigma_2)-psi_path[1][t+1,k]^2/(2*sigma_2^2);
        p1 = p1 + psi_path[1][t+1,k]*T[:,k]'*Diagonal(1 ./psi_path[4])*
        (psi_path[2]-T[:,(d+1):(d+n)]*psi_path[1][t+1,(d+1):(d+n)]-T[:,k]*psi_path[1][t+1,k]/2);

        C = intersect(Active, setdiff(1:d, k));
        if length(C)>0
            p1 = p1 - psi_path[1][t+1,k]*T[:,k]'*Diagonal(1 ./psi_path[4])*T[:,C]*psi_path[1][t+1, C];
        end
        p2 = log((1-p)/sigma_1)-psi_path[1][t+1,k]^2/(2*sigma_1^2);
        dif = p1-p2;
        if dif>10
            prob = 1;
        elseif dif< -10
            prob = 0;
        else
            prob = exp(dif)/(1+exp(dif));
        end
        psi_path[3][t+1,k] = rand(Bernoulli(prob),1)[1];
    end
end

#--------------------------------------------------
#------------------ get values for criteria:
#------------------   1: the minimal of marginal posterior probabilities on the true active covariates
#------------------   2: the maximal of marginal posterior probabilities on the true inactive covariates
#--------------------------------------------------                                                        
#------------------ model selected by MPM: 
#------------------   3: whether the true models is selected
#------------------   4: whether the true models is contained in the selected set
#------------------   5: true positive rate (TPR)
#------------------   6: true positive rate (TPR) given model size <=4
#------------------   7: false discovery rate (FDR)
#------------------   8: model error (beta)
#------------------   9: true positive rate (TPR) for identifying active nodes
#------------------   10: false discovery rate (FDR) for identifying active nodes
#--------------------------------------------------
#------------------ model selected by MAP:
#------------------   11: whether the true models is selected
#------------------   12: whether the true models is contained in the selected set
#------------------   13: true positive rate (TPR)
#------------------   14: true positive rate (TPR) given model size <=4
#------------------   15: false discovery rate (FDR)
#------------------   16: model error (beta)
#------------------   17: true positive rate (TPR) for identifying active nodes
#------------------   18: false discovery rate (FDR) for identifying active nodes
#--------------------------------------------------
p_min = minimum(colmean(psi_path[3][Int(M*0.6):M, option]));
p_max = maximum(colmean(psi_path[3][Int(M*0.6):M, setdiff(1:d, option)]));
                                                        
#--------------------------------------------------
#------------------   model selected by MPM
#--------------------------------------------------
rank_beta = findall(x -> x > 0.5, colmean(psi_path[3][Int(M*0.6):M,1:d]));
rank_top_4 = sortperm(colmean(psi_path[3][Int(M*0.6):M,1:d]), rev = true)[1:length(option)];
if (length(intersect(option, rank_beta ))==length(option))&(length(setdiff(rank_beta , option))==0)
    exact_active_mpm = 1;
else 
    exact_active_mpm = 0;
end
if length(intersect(option, rank_beta ))==length(option)
    include_active_mpm = 1;
else 
    include_active_mpm = 0;
end
tp_beta_mpm = length(intersect(rank_beta , option));
tp_top_4_beta_mpm = length(intersect(rank_top_4, option));
fp_beta_mpm = length(setdiff(rank_beta , option));
beta_true = zeros(Float64, d);
beta_true[option] = signal;
beta_hat = zeros(Float64, d);
beta_hat[rank_beta] = colmean(psi_path[1][Int(M*0.6):M,rank_beta ]);
me_mpm = (beta_true-beta_hat)'*(beta_true-beta_hat);

#--------------------------------------------------
#------------------   model selected by MAP
#--------------------------------------------------
rank_beta = MAP(psi_path[3][Int(M*0.6):M,1:d]);
index_model_size_4 = findall(x->x==length(option) , vec(sum(psi_path[3][Int(M*0.6):M,1:d], dims=2)));
if length(index_model_size_4) > 0
    rank_top_4 = MAP(psi_path[3][Int(M*0.6):M,1:d][index_model_size_4,:]);
    tp_top_4_beta_map = length(intersect(rank_top_4, option));
else
    tp_top_4_beta_map = 0;
end
if (length(intersect(option, rank_beta ))==length(option))&(length(setdiff(rank_beta , option))==0)
    exact_active_map = 1;
else 
    exact_active_map = 0;
end
if length(intersect(option, rank_beta ))==length(option)
    include_active_map = 1;
else 
    include_active_map = 0;
end
tp_beta_map = length(intersect(rank_beta , option));
fp_beta_map = length(setdiff(rank_beta , option));
beta_true = zeros(Float64, d);
beta_true[option] = signal;
beta_hat = zeros(Float64, d);
beta_hat[rank_beta] = colmean(psi_path[1][Int(M*0.6):M,rank_beta ]);
me_map = (beta_true-beta_hat)'*(beta_true-beta_hat);

count = [p_max, p_min, 
        # MPM
        exact_active_mpm, include_active_mpm, tp_beta_mpm / length(option), tp_top_4_beta_mpm / length(option), 
        fp_beta_mpm / (tp_beta_mpm + fp_beta_mpm), me_mpm,
        # MAP
        exact_active_map, include_active_map, tp_beta_map / length(option), tp_top_4_beta_map / length(option), 
        fp_beta_map / (tp_beta_map + fp_beta_map), me_map];
#--------------------------------------------------
#----------------- write results
#--------------------------------------------------
open(join(["/simulation_result/$code","_nonsp_K=$K", "_n=$n", "_d=$d","_Xtype=$Xtype","_corr=$corr","_sp=$sp","_loop=$iter",".csv"])
      ,"w") do f
    for i in 1:length(count)
        a=count[i];
        write(f, "$a ");
    end
end
