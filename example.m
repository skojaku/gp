function example()

	T = readtable('links_karate.dat', 'Delimiter', '\t', 'HeaderLines',0);
	T = table2array(T);
	N = max([max(T(:,2)),max(T(:,1))]);
	A = sparse(T(:,1), T(:,2), T(:,3), N,N);
	[r,c,v] = find(triu(A,1));

	cids = gp(A, 'mod',2, 1, 'louvain');
	
        [r,c,v] = find(triu(A,1));
	edges = [r,c,v];

        qfunc = 'mod';
        sfunc = 'nodes';
        cmd_alg = 'louvain';
        num_of_runs = 1;
        num_of_rand_nets = 100;

        [pvals, s, q, shat, qhat] = qstest_mex( edges, size(A,1), length(r), cids, qfunc, sfunc, cmd_alg, num_of_runs, num_of_rand_nets );
	pvals
		
	
end
