function [gIdx,c, closest]=rvq(X,k)
%RVQ, also uniform density k-means (uniform in each cluster).
% note that iterations are hardcoded to 100.
% and that the delta checks for less than 5% change in cluster membership
% inputs:
% X - data set of points to be reduced
% k - the number of points to reduce to
% outputs:
% gIdx - the cluster membership of each input in the original data set
% c - the centroid locations of the reduced points
% counts - the number of points in each center
% closest - the actual clostest center to each input in the original data


itermax = 30;

[n,m]=size(X);

c=X(ceil(rand(k,1)*n),:);
% allocating variables
g0=ones(n,1);
gIdx=zeros(n,1);
D=zeros(n,k);
counts = zeros(k,1);
for t=1:k
    counts(t) = size(X(gIdx==t,:),1);
end
err = [];

iter = 0;
% Main loop converge if previous partition is the same as current
while sum(sum(g0~=gIdx))/n > .05 && iter < itermax
    %     disp(sum(g0~=gIdx))
    g0=gIdx;
    iter = iter + 1;
    % Loop for each centroid
    for t=1:k
        d=zeros(n,1);
        % Loop for each dimension
        for s=1:m
            d=d+(X(:,s)-c(t,s)).^2;
        end
        D(:,t)=d;
    end
    
    %new lines
    [~ ,indexSorted] = sort(D);
    %end new lines
    ixLoc = ones(k,1);
    isused = gIdx *0;
    
    % pick closest point not already assigned for each cluster
    while min(isused) < 1
        for i=1:k
            loc = indexSorted(ixLoc(i), i);
            while isused(loc) == 1
                ixLoc(i) = ixLoc(i) + 1;
                loc = indexSorted(ixLoc(i), i);
            end
            gIdx(loc) = i;
            isused(loc) = 1;
        end
    end
    
    %     [~,gIdx]=min(D,[],2);
    % Update centroids using means of partitions
    %     err = [err; sum(sum(gIdx ~= g0))];
    for t=1:k
        c(t,:)=mean(X(gIdx==t,:));
    end
    
    
%     for t=1:k
%         counts(t) = size(X(gIdx==t,:),1);
%     end
    
    %     if mod(iter,10) == 0
    %        figure(iter);scatter(c(:,1),c(:,2),counts*5);
    %     end
    
end
% closest = 0;
% debug info
for t=1:k
    d=zeros(n,1);
    % Loop for each dimension
    for s=1:m
        d=d+(X(:,s)-c(t,s)).^2;
    end
    D(:,t)=d;
end
[~,closest]=min(D,[],2);
% disp('iters')
% disp(iter);