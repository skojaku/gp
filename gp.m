classdef gp 
 
	methods (Access = public)
	
		function [cids, Qs, h, pvals, param] = detect(self, A, param)
			
			qfunc = param.qfunc;
			K = param.K;
			num_of_runs = param.num_of_runs;
			algorithm = param.algorithm;
		        sfunc = param.sfunc;
		        num_of_rand_nets = param.num_of_rand_nets;
			
			% community detection	
			[r,c,v] = find(triu(A,1));
			edges = [r,c,v];
			N = size(A, 1);
			
			[cids, Qs] = gp_mex(edges, N, length(r), K, num_of_runs, qfunc, algorithm);
			
			% significance test
			if param.significance_level < 1
		        	sfunc = param.sfunc;
		        	num_of_rand_nets = param.num_of_rand_nets;
				alpha = param.significance_level;
			
		        	[pvals, s, q, shat, qhat] = qstest_mex( edges, N, size(edges,1), cids, qfunc, sfunc, algorithm, num_of_runs, num_of_rand_nets );
				
				corrected_alpha = 1.0 - (1.0 - alpha)^(1.0 / max(cids) ); 
				h = full(sparse(1:N, 1, pvals(cids) <=corrected_alpha ));
			
				param.corrected_alpha = corrected_alpha;
				param.shat = shat;
				param.qhat = qhat;
				param.s = s;
				param.q = q;
			end	
			
		end
		
		function param = init(self)
			param = struct;
			param.qfunc = 'dcsbm';	
			param.algorithm = 'kl';	
			param.K = 2;	
			param.num_of_runs = 1;	
			param.significance_level = 0.05;
			param.sfunc = 'edges';	
			param.num_of_rand_nets = 500;
		end
	end
end
