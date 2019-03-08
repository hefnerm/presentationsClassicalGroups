#This method returns a StraightLineProgram with relations for a presentation of the alternating group on n elements if n is odd, more exactly the relations s^(n-2), t^3, (st)^n, (ts^(-k)ts^k)^2, 1 <= k <= (n-3)/2, for two inputs s,t.
GetAnPresAsSLPodd := function(n)
	local lines, returnValues, k;
	
	#create a list of the first elements of the SLP; the element s will be at the first position and t at the second
	lines := [[[1, 1, 2, 1], 3], #save st as the third element of the SLP
	[[2, 1, 1, -1], 4], #save ts^(-1) as the fourth element of the SLP
	[[4, 1, 2, -1], 5], #save ts^(-1)t^(-1) as the fifth element of the SLP
	[[4, 1, 2, 1, 1, 1], 6]]; #ts^(-1)ts
	
	#create a list, where the relations to return shall be stored
	#s^(n-2), t^3, (st)^n, (ts^(-1)ts^1)^2
	returnValues := [[1, n-2], [2, 3], [3, n], [6, 2]];
	
	#append the lines for creating the relations (ts^(-k)ts^k)^2, 2 <= k <= (n-3)/2
	for k in [2..(n-3)/2] do
		Add(lines, [[5, 1, 4+k, 1, 1, 1], 5+k]); #ts^(-k)ts^k
		
		#store the relation (ts^(-k)ts^k)^2
		Add(returnValues, [5+k, 2]);
	od;
	
	Add(lines, returnValues);
	
	return StraightLineProgram(lines);
end;

#This method returns a StraightLineProgram with relations for a presentation of the alternating group on n elements if n is even, more exactly the relations s^(n-2), t^3, (st)^(n-1), (t^((-1)^k)s^(-k)ts^k)^2, 1 <= k <= (n-3)/2, for two inputs s,t.
GetAnPresAsSLPeven := function(n)
	local lines, returnValues, k;
	
	#create a list of the first elements of the SLP; the element s will be at the first position and t at the second
	lines := [[[1, 1, 2, 1], 3], #save st as the third element of the SLP
	[[2, 1, 1, -1, 2, 1], 4], #save ts^(-1)t as the fourth element of the SLP
	[[3, -1, 2, -1], 5], #save t^(-1)s^(-1)t^(-1) as the fifth element of the SLP
	[[3, -1, 2, 1, 1, 1], 6]]; #t^(-1)s^(-1)ts
	
	#create a list, where the relations to return shall be stored
	#s^(n-2), t^3, (st)^(n-1), (t^(-1)s^(-1)ts^1)^2
	returnValues := [[1, n-2], [2, 3], [3, n-1], [6, 2]];
	
	#append the lines for creating the relations (t^((-1)^k)s^(-k)ts^k)^2, 2 <= k <= (n-3)/2
	for k in [2..(n-2)/2] do
		#distinguish two cases, either the exponent of t in the relation shall be -1 or 1 and the previously stored element is ts^(-(k-1))ts^(k-1) or t^(-1)s^(-(k-1))ts^(k-1)
		if IsEvenInt(k) then
			Add(lines, [[4, 1, 4+k, 1, 1, 1], 5+k]); #ts^(-k)ts^k
		else
			Add(lines, [[5, 1, 4+k, 1, 1, 1], 5+k]); #t^(-1)s^(-k)ts^k
		fi;
		
		#store the relation (t^((-1)^k)s^(-k)ts^k)^2
		Add(returnValues, [5+k, 2]);
	od;
	
	Add(lines, returnValues);
	
	return StraightLineProgram(lines);
end;

#This method returns a StraightLineProgram with relations for a presentation of the alternating group on d elements.
GetAnPresAsSLP := function(d)
	if d > 3 then
		if IsOddInt(d) then
			return GetAnPresAsSLPodd(d);
		else
			return GetAnPresAsSLPeven(d);
		fi;
	else
		Error("argument must be larger than 3!");
	fi;
end;

#This method returns true if two inputs s,t fulfill the presentation used in GetAnPresAsSLP for the alternating group on d elements, or false if not.
#The presentation can be found in COXETER, Harold SM; MOSER, William OJ. Generators and relations for discrete groups. Springer Science & Business Media, 2013, page 67
CheckGeneratorsForAn := function(s, t, d)
	local SLP, ResultSLP, elm;
	SLP := GetAnPresAsSLP(d);
	
	ResultSLP := ResultOfStraightLineProgram(SLP, [s, t]);
	for elm in ResultSLP do
		
		#the relations returned by the StraightLineProgram have to be fulfilled
		if not IsOne(elm) then
			return false;
		fi;
	od;
	
	return true;
end;

#This method does tests if the CheckGeneratorsForAn and GetAnPresAsSLP methods work correctly for s = (1,2)*(3,..,d) (d even) or s = (3,..,d) (d odd)  and t = (1,2,3) for d in [4..500]. 
#This method returns nothing if the Methods work correctly or throws an Error otherwise
CheckPresAn := function()
	local s, t, d, sList;
	
	t := (1,2,3);
	for d in [4..500] do
		sList := [1,2];
		Append(sList, [4..d]);
		Add(sList, 3);
		if IsOddInt(d) then
			s := PermList(sList);
		else
			s := (1,2)*PermList(sList);
		fi;
		if not CheckGeneratorsForAn(s,t,d) then
			Error("FALSE");
		fi;
	od;
end;

#This method compares the runtime of the GetAnPresAsSLP method for computing with cycles or computing with matrices of a subgroup of some GL(d,q) by using a random transformation matrix
TestRunTimeAn := function(d, q)
	local rand, s, t, permMats, permMatt, transformMats, transformMatt, sList, start_time;
	rand := PseudoRandom(GL(d, q));
	t := (1,2,3);
	sList := [1,2];
	Append(sList, [4..d]);
	Add(sList, 3);
	if IsOddInt(d) then
		s := PermList(sList);
	else
		s := (1,2)*PermList(sList);
	fi;
	permMats := PermutationMat(s, d, GF(q));
	permMatt := PermutationMat(t, d, GF(q));
	transformMats := rand^-1*permMats*rand;
	transformMatt := rand^-1*permMatt*rand;
	start_time := Runtime();
	
	CheckGeneratorsForAn(transformMats, transformMatt, d);
	Print("time elapsed while computing as a subgroup of GL(d,q): ", Runtime() - start_time, "ms\n");
	
	start_time := Runtime();
	CheckGeneratorsForAn(s, t, d);
	Print("time elapsed while computing with cycles: ", Runtime() - start_time, "ms\n");
end;

