% single-trial mediation
clear
close all
clc

% define subject names
allsubj = {'blah1','blah2','blah3'};

% these are dummy variables for mediation, single-trial values 
subjects = zeros(1,1000); % this will be a trial-by-trial index of the subject
sides = zeros(1,1000); % this will be a trial-by-trial index of the target side
N2cs = zeros(1,1000); % this will be trial-by-trial N2c amplitude
CPPslopes = zeros(1,1000); % this will be trial-by-trial CPP slope
CPPamps = zeros(1,1000); % this will be trial-by-trial CPP amplitude
RTs = zeros(1,1000); % this will be trial-by-trial RT

% zscore these values inside subject and condition
for s = 1:length(allsubj)
    for side = 1:2
        indx = find(subjects==s & sides==side);
        N2cs(indx) = zscore(N2cs(indx));
        CPPslopes(indx) = zscore(CPPslopes(indx));
        CPPamps(indx) = zscore(CPPamps(indx));
        RTs(indx) = zscore(RTs(indx));
    end
end
        
% perform single-trial mediation, N2c amp -> CPP slope -> RT, covarying for
% CPP amplitude. Make sure the vectors are in the correct dimensions! Might
% take a little messing around. It will be some combination of either
% "variable" and "variable'". See for example below, I've inverted N2cs to
% N2cs'
[paths, toplevelstats, firstlevelstats] = mediation(N2cs',RTs',CPPslopes,'stats','covs',CPPamps');

% x = N2c
% y = RT
% m = CPP slope
% cov = CPP amp
% 1 a   X -> M relationship
% 2 b   M -> Y relationship
% 3 cp  unmediated X -> Y relationship (residual)
% 4 c   X -> Y relationship
% 5 ab  mediated X -> Y by M (a * b)
% [paths, toplevelstats, 1stlevelstats] = mediation(X, Y, M, [stats], [plots], [other optional args])
