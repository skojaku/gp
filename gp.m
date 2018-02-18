function [cids,Qs] = gp(varargin)
	Defaults = {[],'dcsbm',2,1,'kl'};
	Defaults(1:nargin) = varargin;

	A = Defaults{1};
	qfunc = Defaults{2};
	K = Defaults{3};
	numRun = Defaults{4};
	algorithm = Defaults{5};
		
	[r,c,v] = find(triu(A,1));
	[cids, Qs] = gp_mex([r,c,v], size(A,1), length(r), K, numRun, qfunc, algorithm);
end
