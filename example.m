function test()

	T = readtable('links_karate.dat', 'Delimiter', '\t');
	N = max([max(T.source),max(T.target)]);
	A = sparse(T.source, T.target, T.value, N,N);
	 
	[r,c,v] = find(triu(A,1));

	
	numRun = 30;
	K = 2;
	qfunc = 'dcsbm';
	algorithm = 'kl';
 	
	[cids, Qs] = graph_partitioning_mex([r,c,v], size(A,1), length(r), K, numRun, qfunc, algorithm);
	sum(Qs)
	Qs
	C = sparse(1:length(cids),cids,1);
	Q = sum(Qs);
	score = ones(N,1) * Q/N;
	
	full(C)
end
