#--------------------------------------------------
#--- Variable Selection for High-dimensional Nodal Attributes 
#--- in Social Networks with Degree Heterogeneity
#--------------------------------------------------
#--------- Real data analysis (paper citation network)
#--------- BSM-net.sp: written in Julia 1.3.1
#--------------------------------------------------

#--------------------------------------------------
#----------------- set up functions and parameters
#--------------------------------------------------
using Distributions, Roots, StatsBase, CPUTime, Clustering, Random, LinearAlgebra, DelimitedFiles
M = 1000;

function g(x,y) 
    return x.*y
end

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

#--------------------------------------------------
#----------------- read X and Y
#--------------------------------------------------
code = 1207;
file_name = join(["/citation_X_$code","_n=476.csv"]);
X = readdlm(file_name, ',');
X = X[2:size(X)[1], 2:size(X)[2]];
file_name = join(["/citation_Y_$code","_n=476.csv"]);
Y = readdlm(file_name, ',');
Y = Y[2:size(Y)[1], 2:size(Y)[1]];
n,d = size(X);

# standardized X
for i in 1:size(X)[2]
    m = mean(X[:,i]);
    sd = sqrt(sum((X[:,i] .- m) .^2) / (length(X[:,i]) - 1));
    X[:,i]= (X[:,i] .- m)/sd;
end

#--------------------------------------------------
#------------------ set up the parameters for the MCMC 
#--------------------------------------------------
function f_p(x)
    10/d-x-1.29*sqrt(x*(1-x)/d)
end
p_1 = find_zero(f_p, (0,1)); # q_1n: preliminary inclusion probability
function f_p(x)
    50/n-x-1.29*sqrt(x*(1-x)/n)
end
p_2 = find_zero(f_p, (0,1)); # q_2n: preliminary probability of being an active node
# get a rough restimation for tau
p_vec = sort(colmean(Y));
p_hat = percentile(p_vec, 80);
tau = log((1-p_hat)/p_hat)/log(n); 
N = n*(n-1)/2;
sigma_1 = 0.5; # sd for spike prior (beta)
s = log(n/2)/log(n); 
sigma_2 = max( exp((2.01+2*s)*n^s*log(n)/600)/n, log(n) ) ; # sd for spike prior (delta)
sigma_d = log(n)*sqrt(n)/20; # sd for slab prior (beta)
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
z_default =  log(maximum(p_vec)/(1-maximum(p_vec))); 
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
#----------------- write results
#--------------------------------------------------
gamma = psi_path[3];
file_name = join(["/gamma_sp_product_standardized_$code","_n=$n","_M=$M.csv"]);
open(file_name,"w") do f
    for i in 1:size(gamma)[1]
        for j in 1:size(gamma)[2]
            a = gamma[i,j];
            write(f, "$a ");
        end
        write(f, "\n");
    end
end

beta_delta = psi_path[1];
file_name = join(["/beta_delta_sp_product_standardized_$code","_n=$n","_M=$M.csv"]);
open(file_name,"w") do f
    for i in 1:size(beta_delta)[1]
        for j in 1:size(beta_delta)[2]
            a = beta_delta[i,j];
            write(f, "$a ");
        end
        write(f, "\n");
    end
end
