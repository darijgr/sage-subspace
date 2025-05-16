attach("combinatorial_vector_tuple.py")

@cached_function
def QS(field, n):
    return SymmetricGroupAlgebra(field, n)

def I22(field):
    SGA4 = QS(field, 4)
    a = SGA4.sum(SGA4(w) for w in Tableau([[1,2],[3,4]]).row_stabilizer())
    return VectorTuple([SGA4(u)*a*SGA4(v) for u in Permutations(4) for v in Permutations(4)]).span_dimension()

print(I22(QQ))
print(I22(GF(2)))