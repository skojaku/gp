function example()

	T = readtable('links_karate.dat', 'Delimiter', '\t', 'HeaderLines',0);
	T = table2array(T);
	N = max([max(T(:,2)),max(T(:,1))]);
	A = sparse(T(:,1), T(:,2), T(:,3), N,N);
	[r,c,v] = find(triu(A,1));

	cids = gp(A, 'dcsbm',2)
	
end
