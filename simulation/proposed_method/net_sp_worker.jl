#--------------------------------------------------
#--- Variable Selection for High-dimensional Nodal Attributes 
#--- in Social Networks with Degree Heterogeneity
#--------------------------------------------------
#--------- Simulation Study
#--------- BSM-net.sp: written in Julia 1.3.1, last edited in Aug 27, 2022
#--------------------------------------------------

iter = ARGS[1]; # get the number of replication
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

# a function to get the row with the highest frequence in a matrix
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

# a function to reshape the design matrix to a vector based on the current node set
function reshape_T_vec(T, Model)
    l = length(Model);
    temp = T_full[Model, Model, vcat(1:d, Model .+d)];
    temp = reshape(temp, (l*l, l+d));
    temp = temp[.!vec(any(isnan.(temp), dims=2)), : ];
    return temp
end

# define covariance matrix for X
if Xtype==2 # compound symmetric matrix
    sigma = reshape(fill(corr, d*d), d,d)+Diagonal(fill(1-corr, d));
else  # autocorrelation
    sigma = zeros(d,d);
    for i in 1:d
        for j in 1:d
            sigma[i,j]=corr^(abs(i-j));
        end
    end
end

function g(x,y)
    return x.*y
end
#--------------------------------------------------
#----------------- simulate X and Y
#--------------------------------------------------
X = rand(MultivariateNormal(zeros(d), sigma), n)';
influence = sortperm(X[:,1]+ rand(Normal(0,4),n), rev = true)[1:n_a];
delta_true = rand(Uniform(-1, 1),n);
delta_true[setdiff(1:n, influence)] =  delta_true[setdiff(1:n, influence)] + sp * 2 * rand(Uniform(-4,-2), n-length(influence));

Y = reshape(fill(0,n*n), n, n);
for i = 2:n
    for j = 1:(i-1)
        pp = exp(delta_true[i]+delta_true[j] + dot(g(X[i,option], X[j,option]),signal));
        Y[i,j] = rand(Bernoulli(pp/(1+pp)),1)[1];
        Y[j,i] = Y[i,j];
    end
end

#--------------------------------------------------
#----------------- write simulated X, Y, delta
#--------------------------------------------------
data = X;
file_name = join(["/storage/home/j/jzw88/work/Research_NET_revision/simulation_data/$code","_X_K=$K", "_n=$n", "_d=$d","_Xtype=$Xtype","_corr=$corr","_sp=$sp","_loop=$iter",".csv"]);
open(file_name,"w") do f
    for i in 1:size(data)[1]
        for j in 1:size(data)[2]
            a = data[i,j];
            write(f, "$a ");
        end
        write(f, "\n");
    end
end
data = Y;
file_name = join(["/storage/home/j/jzw88/work/Research_NET_revision/simulation_data/$code","_Y_K=$K", "_n=$n", "_d=$d","_Xtype=$Xtype","_corr=$corr","_sp=$sp","_loop=$iter",".csv"]);
open(file_name,"w") do f
    for i in 1:size(data)[1]
        for j in 1:size(data)[2]
            a = data[i,j];
            write(f, "$a ");
        end
        write(f, "\n");
    end
end
data = delta_true;
file_name = join(["/storage/home/j/jzw88/work/Research_NET_revision/simulation_data/$code","_delta_K=$K", "_n=$n", "_d=$d","_Xtype=$Xtype","_corr=$corr","_sp=$sp","_loop=$iter",".csv"]);
open(file_name,"w") do f
    for i in 1:size(data)[1]
        a = data[i];
        write(f, "$a ");
    end
end
#--------------------------------------------------
#------------------ set up the parameters for the MCMC 
#--------------------------------------------------
function f_p(x)
    K/d-x-1.29*sqrt(x*(1-x)/d)
end
p_1 = find_zero(f_p, (0,1)); # q_1n: preliminary inclusion probability
n_a_hat = 30;
function f_p(x)
    n_a_hat/n-x-1.29*sqrt(x*(1-x)/n)
end
p_2 = find_zero(f_p, (0,1)); # q_2n: preliminary probability of being an active node
# get a rough restimation for tau
p_vec = sort(colmean(Y));
p_hat = percentile(p_vec, 80);
tau = log((1-p_hat)/p_hat)/log(n); 
N = n*(n-1)/2;
sigma_1 = 0.5; # sd for spike prior (beta)
s = log(n_a_hat)/log(n); 
sigma_2 = max( exp((2.01+2*s)*n^s*log(n)/600)/n, log(n) ) ; # sd for slab prior (beta)
sigma_d = log(n)*sqrt(n)/20; # sd for slab prior (delta)
sigma_1d = 1.96*exp(-tau*log(n)/2); # sd for slab prior (delta)
v = 7.3;
w = pi*sqrt((v-2)/3/v);

# initialization
T_full = reshape(fill(NaN, n*n*(n+d)), (n,n,d+n));
for i = 2:n
    for j = 1:(i-1)
        b = fill(0.0, n);
        b[i] = 1;
        b[j] = 1;
        a = vcat(g(X[i,:],X[j,:]), b);
        T_full[j,i,:] = a ;
    end
end

z_init = reshape(fill(NaN, n*n), (n,n));
psi_init = reshape(fill(NaN, n*n), (n,n));
z_default =  log(percentile(p_vec, 90)/(1-percentile(p_vec, 90)));
for i in 2:n
    for j in 1:(i-1)
        z_init[j,i] = z_default;
        psi_init[j,i] = 1.0;
    end 
end 
psi_path=[zeros(Float64, M+1, d+n), z_init, zeros(Int32, M+1, d+n), psi_init];
psi_path[3][1,(d+1):(d+n)] = ones(Float64, n);

#--------------------------------------------------
#------------------ start the iteration to update :  
#------------------      beta, delta, z (g in the paper), psi, gamma, eta
#--------------------------------------------------

@CPUtime for t in 1:M # record the CPU time

    Active = findall(x -> x == 1, psi_path[3][t, 1:d]);
    Model = findall(x -> x == 1, psi_path[3][t, (d+1):(d+n)]); 
    l = length(Model);
    T = reshape_T_vec(T_full, Model);
    psi_vec_model = vec(reshape(psi_path[4][Model, Model], l*l, 1));
    psi_vec_model = psi_vec_model[.!vec(isnan.(psi_vec_model))];
    z_vec_model = vec(reshape(psi_path[2][Model, Model], l*l, 1));
    z_vec_model = z_vec_model[.!vec(isnan.(z_vec_model))];

    #--------------------------------------------------
    #------------------ update beta and delta
    #--------------------------------------------------
    D = Diagonal(1 ./vcat(fill(sigma_2^2, length(Active)), fill(sigma_d^2, length(Model))));
    V = D+T[:,vcat(Active, (d+1):(d+l))]'*Diagonal(1 ./psi_vec_model)*T[:,vcat(Active, (d+1):(d+l))];
    Sigma = inv(V);
    m = Sigma*T[:,vcat(Active, (d+1):(d+l))]'*Diagonal(1 ./psi_vec_model)*z_vec_model;
    Sigma = (Sigma+Sigma')/2;
    psi_path[1][t+1,vcat(Active,d .+ Model)] = (Sigma^(0.5)*rand(Normal(0,1), length(Active)+ length(Model))+m)';

    # update inactive beta
    Inactive = findall(x -> x == 0, psi_path[3][t,1:d]);
    if length(Inactive)>0
        psi_path[1][t+1,Inactive] = rand(Normal(0,sigma_1), length(Inactive))';
    end

    # update inactive delta
    Inactive = findall(x -> x == 0, psi_path[3][t,(d+1):(d+n)]);
    if length(Inactive)>0
        psi_path[1][t+1,d .+ Inactive] = log.(rand(truncated(Normal(0,sigma_1d), 0, Inf), length(Inactive))');
    end

    #--------------------------------------------------
    #------------------ update z and psi
    #--------------------------------------------------
    m_mat = reshape(fill(-tau*log(n), n*n) , n, n );
    if length(Active)>0
        for i in 2:l
            for j in 1:(i-1)
                m_mat[Model[i],Model[j]] = 1* (psi_path[1][t+1,(d+Model[j])]+psi_path[1][t+1,(d+Model[i])])+ dot(g(X[Model[i],Active],X[Model[j],Active]), psi_path[1][t+1, Active]);
                m_mat[Model[j],Model[i]] = m_mat[Model[i],Model[j]];
            end
        end
        else
        for i in 2:l
            for j in 1:(i-1)
                m_mat[Model[i],Model[j]] = 1* (psi_path[1][t+1,(d+Model[j])]+psi_path[1][t+1,(d+Model[i])]);
                m_mat[Model[j],Model[i]] = m_mat[Model[i],Model[j]] ;
            end
        end
    end

    for i = 2:n
        for j = 1:(i-1)
            if Y[i,j]==1
                psi_path[2][j,i] = rand(truncated(Normal(m_mat[j,i],sqrt(psi_path[4][j,i])),0, Inf), 1)[1];
            else
                psi_path[2][j,i] = rand(truncated(Normal(m_mat[j,i],sqrt(psi_path[4][j,i])),-Inf,0), 1)[1];
            end
            psi_path[4][j,i] = rand(InverseGamma(v/2+1/2, w^2*v/2+(psi_path[2][j,i]-m_mat[j,i])^2/2) ,1)[1];
        end
    end
    z_mat = UpperTriangular(psi_path[2])+UpperTriangular(psi_path[2])';
    psi_mat = UpperTriangular(psi_path[4])+UpperTriangular(psi_path[4])';

    #--------------------------------------------------
    #------------------ update gamma
    #--------------------------------------------------
    for k in 1: d
        p1 = log(p_1/sigma_2)-psi_path[1][t+1,k]^2/(2*sigma_2^2);
        p1 = p1 + psi_path[1][t+1,k]*T[:, k]'*Diagonal(1 ./psi_vec_model)*
        (z_vec_model- T[:,(d+1):(d+l)] * psi_path[1][t+1, d .+ Model] -T[:, k]*psi_path[1][t+1,k]/2);

        C = intersect(Active, setdiff(1:d, k));
        if length(C)>0
            p1 = p1 - psi_path[1][t+1,k]*T[:, k]'*Diagonal(1 ./psi_vec_model)*T[:, C ]*psi_path[1][t+1, C];
        end
        p2 = log((1-p_1)/sigma_1)-psi_path[1][t+1,k]^2/(2*sigma_1^2);
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

    #--------------------------------------------------
    #------------------ update eta
    #--------------------------------------------------
    for k in 1: n
        C = intersect(Model, setdiff(1:n, k));
        a = 0;
        if length(Active)>0
            for i in C
                a = a +  ((z_mat[i,k]+tau*log(n))^2 - (z_mat[i,k]- psi_path[1][t+1,(d+k)] - psi_path[1][t+1,(d+i)] - dot( g(X[i,Active], X[k,Active]), psi_path[1][t+1, Active]))^2 )/(2*psi_mat[i,k]);
            end
            else
            for i in C
                a = a +  ((z_mat[i,k]+tau*log(n))^2 - (z_mat[i,k]- psi_path[1][t+1,(d+k)] - psi_path[1][t+1,(d+i)])^2 )/(2*psi_mat[i,k]);
            end
        end 
        p1 = log(p_2/sigma_d)-psi_path[1][t+1,k+d]^2/(2*sigma_d^2);
        p1 = p1 + a ;
        p2 = log((1-p_2)*2/sigma_1d)- exp(2*psi_path[1][t+1,k+d])/(2*sigma_1d^2) + psi_path[1][t+1,k+d];
        dif = p1-p2;
        if dif>10
            prob = 1;
        elseif dif< -10
            prob = 0;
        else
            prob = exp(dif)/(1+exp(dif));
        end
        psi_path[3][t+1,k+d] = rand(Bernoulli(prob),1)[1];
    end
end

#--------------------------------------------------
#------------------ get values for criteria:
#------------------   1: the minimal of marginal posterior probabilities on the true active covariates
#------------------   2: the maximal of marginal posterior probabilities on the true inactive covariates
#------------------ model selected by MPM: 
#------------------   3: whether the true models is selected
#------------------   4: whether the true models is contained in the selected set
#------------------   5: true positive rate (TPR)
#------------------   6: true positive rate (TPR) given model size <=4
#------------------   7: false discovery rate (FDR)
#------------------   8: model error (beta)
#------------------   9: true positive rate (TPR) for identifying active nodes
#------------------   10: alse discovery rate (FDR) for identifying active nodes
#------------------ model selected by MAP:
#------------------   11: whether the true models is selected
#------------------   12: whether the true models is contained in the selected set
#------------------   13: true positive rate (TPR)
#------------------   14: true positive rate (TPR) given model size <=4
#------------------   15: false discovery rate (FDR)
#------------------   16: model error (beta)
#------------------   17: true positive rate (TPR) for identifying active nodes
#------------------   18: alse discovery rate (FDR) for identifying active nodes
#--------------------------------------------------
p_min = minimum(colmean(psi_path[3][Int(M*0.6):M, option]));
p_max = maximum(colmean(psi_path[3][Int(M*0.6):M, setdiff(1:d, option)]));

#--------------------------------------------------
#------------------   model selected by MPM
#--------------------------------------------------
rank_beta = findall(x -> x > 0.5, colmean(psi_path[3][Int(M*0.6):M,1:d]));
rank_top_4 = sortperm(colmean(psi_path[3][Int(M*0.6):M,1:d]), rev = true)[1:length(option)];
rank_delta = findall(x -> x > 0.5, colmean(psi_path[3][Int(M*0.6):M,(d+1):(d+n)]));
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
tp_delta_mpm = length(intersect(rank_delta, influence));
fp_delta_mpm = length(setdiff(rank_delta, influence));

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
rank_delta = MAP(psi_path[3][Int(M*0.6):M,(d+1):(d+n)]);
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
tp_beta_map = length(intersect(rank_beta, option));
fp_beta_map = length(setdiff(rank_beta, option));
beta_true = zeros(Float64, d);
beta_true[option] = signal;
beta_hat = zeros(Float64, d);
beta_hat[rank_beta] = colmean(psi_path[1][Int(M*0.6):M,rank_beta]);
me_map = (beta_true-beta_hat)'*(beta_true-beta_hat);
tp_delta_map = length(intersect(rank_delta, influence));
fp_delta_map = length(setdiff(rank_delta, influence));

count = [p_min, p_max, 
        # MPM
        exact_active_mpm, include_active_mpm, tp_beta_mpm / length(option), tp_top_4_beta_mpm / length(option), 
        fp_beta_mpm / (tp_beta_mpm + fp_beta_mpm), me_mpm, tp_delta_mpm/length(influence), fp_delta_mpm / (tp_delta_mpm + fp_delta_mpm),
        # MAP
        exact_active_map, include_active_map, tp_beta_map / length(option), tp_top_4_beta_map / length(option), 
        fp_beta_map / (tp_beta_map + fp_beta_map), me_map, tp_delta_map/length(influence), fp_delta_map / (tp_delta_map + fp_delta_map)];
#--------------------------------------------------
#----------------- write results
#--------------------------------------------------
# write result
open(join(["/storage/home/j/jzw88/work/Research_NET_revision/simulation_result/$code","_sp_K=$K", "_n=$n", "_d=$d", "_Xtype=$Xtype","_corr=$corr","_sp=$sp","_loop=$iter",".csv"])
      ,"w") do f
    for i in 1:length(count)
        a=count[i];
        write(f, "$a ");
    end
end
