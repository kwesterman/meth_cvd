with open("../int/unprunedGenos.csv", "r") as inMat, open("../int/sexGenos.csv", "w") as outMat:
	print(next(inMat).strip().split(',')[1:5])
