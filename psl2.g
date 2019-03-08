#This method returns a StraightLineProgram with relations for a presentation of the projective special linear group on p^e elements, more exactly the relations (t*u)^3, (u*d)^2, u^2, (t1*u*d)^3, d^((p^e-1)/2), t^p, [t, t1], [t1, t^d], t^(d^m) or t^(d^(m+1)) (depends on the boolean input isSquare), t*t1, t1^(d^m), t1*t^d for inputs t, d, u, t1, d^m. The element t^(f), where f is the minimal polynomial of a primitive element over GF(q), is not computed in this SLP.
GetPSL2OddPresAsSLP := function(p, e, isSquare)
	local lines, returnValues;
	
	#create a list of the elements of the SLP; the element t will be at the first position d at the second, u at third, t1 at fourth and d^m at fifth
	lines := [[[1, 1, 3, 1], 6], #save t*u as the sixth element of the SLP
	[[3, 1, 2, 1], 7], #save u*d as the seventh element of the SLP
	[[4, 1, 3, 1, 2, 1], 8], #t1*u*d
	[[2, -1], 9], #d^-1
	[[9, 1, 1, 1, 2, 1], 10], #t^d
	[[4, -1], 11], #t1^-1
	[[5, -1], 12]]; #d^-m
	
	if isSquare then
		Add(lines, [[5, 1, 2, 1], 13]); #d^(m+1)
	fi;
	
	#(t*u)^3, (u*d)^2, u^2, (t1*u*d)^3, d^((p^e-1)/2), t^p, [t, t1], [t1, t^d]
	returnValues := [[6, 3], [7, 2], [3, 2], [8, 3], [2, (p^e-1)/2], [1, p], [1, -1, 11, 1, 1, 1, 4, 1], [11, 1, 10, -1, 4, 1, 10, 1]];
	
	if isSquare then
		Add(returnValues, [12, 1, 1, 1, 5, 1]); #t^(d^m)
	else
		Add(returnValues, [13, -1, 1, 1, 13, 1]); #t^(d^(m+1))
	fi;
	
	#t*t1, t1^(d^m), t1*t^d
	Append(returnValues, [[1, 1, 4, 1], [12, 1, 4, 1, 5, 1], [4, 1, 10, 1]]);
	
	Add(lines, returnValues);
	
	return StraightLineProgram(lines);
end;

BruteforceCalculationOfM := function(w)
	local oddIsRight, evenIsRight, m, save;
	
	m := 0;
	
	while true do
		if w^(2*m) = w+1 then
			return [m, true];
		fi;
		if w^(2*m-1) = w+1 then
			return [m, false];
		fi;
		m := m+1;
	od;
end;

##########Not executable yet

#This method returns true if the inputs t,d,u fulfill the presentation used in GetPSL2OddPresAsSLP for the projective special linear group on p^e elements for p odd and e > 1, or false if not.
#The presentation can be found in LEEDHAM-GREEN, C. R.; O'BRIEN, E. A. Presentations on standard generators for classical groups. arXiv preprint arXiv:1811.00685, 2018, page 9
#in this method t equals tau, d equals delta and u equals U
CheckGeneratorsForPSL2 := function(p, e, t, d, u)
	local q, w, m, ind, minpo, coeffsMinpo, coeffsRepresW, b, wbas, mat, dstore, conj, result, gens, elm, isSquare, SLP, t1, dm, tf, i, list, wstore;
	
	if e <= 1 then
		Error("e shall be larger than 1");
	fi;
	
	if not (IsPrime(p) and p > 2) then
		Error("p shall be an odd prime");
	fi;
	
	q := p^e;
	w := d[2][2];
	ind := Indeterminate(GF(p));
	minpo := MinimalPolynomial(GF(p), w);
	coeffsMinpo := PolynomialCoefficientsOfPolynomial(minpo, ind);
	
	#do a basis transformation to get the coefficients of the linear combination of w in the elements {w^0, w^2, w^4, ..., w^(2*(e-1))}
	b := Basis(GF(q));
	wbas := [];
	for i in [0..(e-1)] do
		Add(wbas, w^(2*i));
	od;
	mat := List(wbas, y -> Coefficients(b, y));
	coeffsRepresW := Coefficients(b, w) * mat^-1;
	
	#get m and the boolean if 1+w is a square in GF(q) or not
	list := BruteforceCalculationOfM(w);
	
	m := list[1];
	isSquare := list[2];
	
	#compute the element t1 := prod_{i=0}^{e-1} t^{a_i d^i}
	#dstore stores the element d^i
	t1 := [[1, 0], [0, 1]];
	dstore := [[1, 0], [0, 1]];
	
	#for performance test if d^m is already computed during comuting t1
	if m <= e-1 then
		for i in [0..e-1] do
			#right now dstore equals d^i
			if m = i then
				dm := dstore;
			fi;
			if not IsZero(coeffsRepresW[i+1]) then
				conj := coeffsRepresW[i+1]*dstore;
				t1 := t1 * (conj)^-1 * t * conj;
			fi;
			dstore := dstore * d;
		od;
	else
		for i in [0..e-1] do
			if not IsZero(coeffsRepresW[i+1]) then
				conj := coeffsRepresW[i+1]*dstore;
				t1 := t1 * (conj)^-1 * t * conj;
			fi;
			dstore := dstore * d; #dstore is set as d^(i+1)
		od;
		
		dm := dstore * d^(m-e);
	fi;
	
	SLP := GetPSL2OddPresAsSLP(p, e, isSquare);
	
	gens := [t, d, u, t1, dm];
	
	result := ResultOfStraightLineProgram(SLP, gens);
	
	#the first 8 elements have to equal 1
	for elm in result{[1..8]} do
		if not IsOne(elm) then
			return false;
		fi;
	od;
	
	if isSquare then
		if not (result[9] = result[10] or result[11] = result[12]) then
			return false;
		fi;
	else
		if not (result[11] = result[10] or result[9] = result[12]) then
			return false;
		fi;
	fi;
	
	#compute t^(minpo) and check if it is != 1
	tf := t^coeffsMinpo[1];
	wstore := 1;
	
	for i in [1..e-1] do
		wstore := wstore * w;
		tf := tf * ([[1, wstore], [0, 1]])^(coeffsMinpo[i+1]);
	od;
	
	if not IsOne(tf) then
		return false;
	fi;
	
	return true;
end;














