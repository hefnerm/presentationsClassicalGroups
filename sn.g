#This method returns a StraightLineProgram with relations for a presentation of the symmetric group on d elements, more exactly the relations u^2, v^d, (uv)^(d-1), (uu^v)^3, (uu^(v^j))^2, 2 <= j <= d/2, for two inputs u,v.
GetSnPresAsSLP := function(d)
	local lines, j, listJ, returnValues;
	
	if d > 2 then
		
		#create a list of the first elements of the SLP; the element u will be at the first position and v at the second
		lines := [[[2, -1], 3], #save v^(-1) as the third element of the SLP
 [[1, 1, 3, 1], 4], #save uv^(-1) as the fourth element of the SLP
 [[1, 1, 2, 1], 5], #uv
 [[4, 1, 5, 1], 6]]; #uu^v
		
		#create a list, where the relations to return shall be stored
		#u^2, v^d, (uv)^(d-1), (uu^v)^3
		returnValues := [[1, 2], [2, d], [5, d-1], [6, 3]];
		
		#append the lines for creating the relations (uu^(v^j))^2, 2 <= j <= d/2
		if d > 3 then
			
			#first create the lines needed for the relation (uu^(v^2))^2
			listJ := [[[4, 1, 3, 1], 7], #uv^(-2)
			[[5, 1, 2, 1], 8], #uv^2
			[[7, 1, 8, 1], 9]]; #uu^(v^2)
			
			Append(lines, listJ);
			
			#store the relation (uu^(v^2))^2
			Add(returnValues, [9, 2]);
			
			#append the lines for creating the relations (uu^(v^j))^2, 3 <= j <= d/2
			for j in [3..Int(d/2)] do
				
				#first create the lines needed for the relation (uu^(v^j))^2
				listJ := [[[7 + (j-3)*3, 1, 3, 1], 10 + (j-3)*3], #uv^(-j)
				[[8 + (j-3)*3, 1, 2, 1], 11 + (j-3)*3], #uv^j 
				[[10 + (j-3)*3, 1, 11 + (j-3)*3, 1], 12 + (j-3)*3]]; #uu^(v^j)
				
				Append(lines, listJ);
				
				#store the relation (uu^(v^j))^2
				Add(returnValues, [12 + (j-3)*3, 2]);
			od;
		fi;
		
		Add(lines, returnValues);
	else
		Error("argument must be larger than 2!");
	fi;
	
	return StraightLineProgram(lines);
end;

#This method returns true if two inputs u,v fulfill the presentation used in GetSnPresAsSLP for the symmetric group on d elements, or false if not.
#The presentation can be found in LEEDHAM-GREEN, C. R.; O'BRIEN, E. A. Presentations on standard generators for classical groups. arXiv preprint arXiv:1811.00685, 2018, page 7
CheckGeneratorsForSn := function(u, v, d)
	local SLP, elm;
	SLP := GetSnPresAsSLP(d);
	
	for elm in ResultOfStraightLineProgram(SLP, [u, v]) do
		
		#the relations returned by the StraightLineProgram have to be fulfilled
		if not IsOne(elm) then
			return false;
		fi;
	od;
	
	return true;
end;

#This method does tests if the CheckGeneratorsForSn and GetSnPresAsSLP methods work correctly for u = (1,2) and v = (1,2,3,...,d) for d in [3..500]. 
#This method returns nothing if the Methods work correctly or throws an Error otherwise
CheckPresSn := function()
	local u, v, d, vList;
	
	vList := [];
	
	u := (1,2);
	
	for d in [3..500] do
		vList := [2..d];
		Add(vList, 1);
		v := PermList(vList);
		
		if not CheckGeneratorsForSn(u,v,d) then
			Error("FALSE");
		fi;
	od;
end;


#This method compares the runtime of the GetSnPresAsSLP method for computing with cycles or computing with matrices of a subgroup of some GL(d,q) by using a random transformation matrix
TestRunTimeSn := function(d, q)
	local rand, u, v, permMatu, permMatv, transformMatu, transformMatv, vList, start_time;
	rand := PseudoRandom(GL(d, q));
	u := (1,2);
	vList := [2..d];
	Add(vList, 1);
	v := PermList(vList);
	permMatu := PermutationMat(u, d, GF(q));
	permMatv := PermutationMat(v, d, GF(q));
	transformMatu := rand^-1*permMatu*rand;
	transformMatv := rand^-1*permMatv*rand;
	
	start_time := Runtime();
	CheckGeneratorsForSn(transformMatu, transformMatv, d);
	Print("time elapsed while computing as a subgroup of GL(d,q): ", Runtime() - start_time, "ms\n");
	
	start_time := Runtime();
	CheckGeneratorsForSn(u, v, d);
	Print("time elapsed while computing with cycles: ", Runtime() - start_time, "ms\n");
end;

