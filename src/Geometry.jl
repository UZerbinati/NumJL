function oddDomain(x,y)
#	
#               (0,2) -- (2,2)
#                 |        |
#    (-2,0) --- (0,0)    (2,0)
#       |                  |
#    (-2,-2) ---------- (2,-2)

	if (-2 <= x <= 0)
		if (-2 <= y <= 0)
			return 1;
		end
	end 
	if (0 <= x <= 2)
		if (-2 <= y <= 2)
			return 1;
		end
	end
	return 0;
end
function UnitCircle(x,y)
	if sqrt(x^2+y^2) <= 1
		return 1;
	end
	return 0;
end
function UnitSquare(x,y)
	if -1 <= x <= 1
		if -1 <= y <= 1
			return 1;
		end
	end
	return 0;
end
