r"""
Various bases of the symmetric group algebra.
The field can be provided, but needs to have
characteristic 0 for many of the bases.

"The notes" will always refer to my notes
arXiv:2507.20706v1.
"""

attach("/mnt/c/worksnap/TEMPrepth/sage-subspace/combinatorial_vector_tuple.py")
# You'll need to change this to the correct path to
# combinatorial_vector_tuple.py.
# Or import the file properly (I'm too lazy).

from itertools import combinations, permutations, product
from sage.combinat.symmetric_group_algebra import a, b, e

@cached_function
def QS(field, n):
    # The n-th symmetric group algebra over the field field.
    return SymmetricGroupAlgebra(field, n)

@cached_function
def SBT(n):
    # List of standard bitableaux of size `n`,
    # in a specific order.
    # The order in which they are listed is NOT the
    # one in §6.7 of my notes. Rather, it runs through
    # the partitions in decreasing lexicographic order,
    # then runs through the standard tableaux of each
    # partition in whatever order
    # ``StandardTableaux.__iter__`` gives.
    return [(lam, U, V) for lam in Partitions(n)
            for U in reversed(StandardTableaux(lam))
            for V in reversed(StandardTableaux(lam))]

# Standard basis for n = 3 over QQ, ordered as in the notes:
stdbas3 = [QS(QQ, 3)(Permutation(x))
           for x in [[1,2,3], [2,1,3], [1,3,2], [3,2,1], [2,3,1], [3,1,2]]]

# Standard basis in general, different ordering:
def stdbas(field, n):
    return [QS(field, n)(Permutation(w))
            for w in Permutations(n)]

def wuv(U, V):
    # w_{U,V} for two n-tableaux U and V of equal shape.
    # This is the permutation of {1,2,...,n}
    # that yields U when applied to V.
    Uword = [a for row in U for a in row]
    Vword = [b for row in V for b in row]
    n = len(Uword)
    return Permutation([Uword[Vword.index(i)] for i in range(1, n+1)])

@cached_function
def EUV(field, U, V):
    # E_{U,V} for two n-tableaux U and V of equal shape.
    # See Definition 5.15.4 in my notes.
    CU = b(U)
    RV = a(V)
    SGA = parent(CU)
    w = wuv(U, V)
    w = SGA(w)
    return SGA.left_action_product(SGA.left_action_product(CU, w), RV)

@cached_function
def E(field, U):
    # The Young symmetrizer E_U for a given tableau U.
    return EUV(field, U, U)

# E_{U,V}-basis (= Young symmetrizer basis;
# see §6.4 in my notes):
@cached_function
def EUVbas(field, n):
    return [EUV(field, U, V) for (lam, U, V) in SBT(n)]

# E_{U,V}-basis for n = 3 over Q, labelled as in the notes:
EUVbas3 = EUVbas(QQ, 3)

def FwE(field, U, V, old=True):
    # F_U w_{U,V} E_V for two given tableaux U and V
    # of the same shape. See Theorem 5.21.2 and §6.5 in
    # my notes.
    RV = a(V)
    CV = b(V)
    absCV = sum(1 for _ in CV) # size of the column-symmetrizer group
    return absCV * RV * CV * RV * parent(CV)(wuv(U, V))

# F_U w_{U,V} E_V basis (§6.5 in my notes):
@cached_function
def FwEbas(field, n):
    return [FwE(field, U, V) for (lam, U, V) in SBT(n)]

def MurRow(field, U, V):
    # row-Murphy element \Nabla^\Row_{U,V} for two
    # given tableaux U and V of the same shape.
    # See Definition 6.8.2 in my notes.
    RV = a(V)
    SGA = parent(RV)
    return RV * SGA(wuv(U,V))

# row Murphy basis (§6.8 in my notes):
@cached_function
def MurRowbas(field, n):
    return [MurRow(field, U, V) for (lam, U, V) in SBT(n)]

def coeffs(basis, vector):
    # Coefficients of the given vector in the given basis.
    B = VectorTuple(basis)
    return B.coefficients(vector)

def coeff_matrix(basis):
    # Matrix containing the coefficients of the vectors
    # of ``basis`` in the standard basis.
    field = parent(basis[0]).base_ring()
    return Matrix(field,
                  [v.dense_coefficient_list() for v in basis])

def rep_matrix(basis, mapping, other_basis=None):
    # Matrix representing the linear map ``mapping`` in the
    # given basis ``basis``.
    # If ``other_basis`` is provided, then this is the
    # basis in which the images are represented.
    # If ``mapping`` is not provided, then the identity is
    # taken (so the result is a change-of-basis matrix).
    field = parent(basis[0]).base_ring()
    if other_basis is None:
        other_basis = basis
    if mapping is None:
        mapping = lambda v: v
    C = VectorTuple(other_basis)
    M = Matrix(field, [C.coefficients(mapping(v)) for v in basis])
    return M.transpose()

def rep_matrix_S(basis):
    # Matrix representing the antipode map `S` with
    # respect to the given basis ``basis``.
    return rep_matrix(basis, lambda v: v.antipode())

def Tsign(a):
    # The sign-twist of `a`, where `a` is in the
    # symmetric group algebra.
    # This is the result of multiplying every permutation
    # in `a` with its sign.
    SGA = parent(a)
    return SGA.sum_of_terms((m, m.sign() * c) for m, c in a)

def rep_matrix_Tsign(basis):
    # Matrix representing the sign-twist map with
    # respect to the given basis ``basis``.
    return rep_matrix(basis, Tsign)

def si(n, i):
    # The simple transposition `s_i` in `S_n`.
    return Permutation(list(range(1, i)) + [i+1,i] + list(range(i+2, n+1)))

def rep_matrix_L(a, basis):
    # The matrix representing left multiplication by
    # the given element `a` in the given basis ``basis``.
    return rep_matrix(basis, lambda v : v * a)

def rep_matrix_R(a, basis):
    # The matrix representing right multiplication by
    # the given element `a` in the given basis ``basis``.
    return rep_matrix(basis, lambda v : a * v)

def rep_matrix_Lsi(i, basis):
    # The matrix representing left multiplication by
    # `s_i` in the given basis ``basis``.
    SGA = parent(basis[0])
    n = SGA.n
    s = SGA(si(n, i))
    return rep_matrix(basis, lambda v : v * s)

def rep_matrix_Rsi(i, basis):
    # The matrix representing right multiplication by
    # `s_i` in the given basis ``basis``.
    SGA = parent(basis[0])
    n = SGA.n
    s = SGA(si(n, i))
    return rep_matrix(basis, lambda v : s * v)

@cached_function
def EparUn(field, n):
    # the list of the coefficients \gamma_{U,V}
    # such that
    # \sum_{U,V \in \SBT(n)} \gamma_{U,V} E_{U,V} = 1
    VT = VectorTuple(EUVbas(QQ, n))
    return VT.coefficients(QS(QQ, n).one())

def check_EparUn(n):
    # Just a unit test that EparUn does what it should.
    cs = EparUn(QQ, n)
    bas = EUVbas(QQ, n)
    return sum(c * bas[i]
               for (i, c) in enumerate(cs)) == 1

def W(field, U, V):
    # W_{U,V}.
    # this is the element sending the standard polytabloid e_V
    # to e_U and all other standard polytabloids to 0.
    # WARNING: In Sage's Specht module implementation, this is
    # rather the *antipode* of that element, since Sage's
    # Specht modules are right.
    lam = Partition([len(row) for row in U])
    n = sum(lam)
    SGA = QS(field, n)
    gamma = EparUn(field, n)
    SBTn = SBT(n)
    return SGA.sum(gamma[SBTn.index((lam, V, Q))] * EUV(field, U, Q)
                   for Q in StandardTableaux(lam))

# W_{U,V} (= Artin-Wedderburn-Specht) basis:
@cached_function
def Wbas(field, n):
    return [W(field, U, V) for (lam, U, V) in SBT(n)]

def check_Wbas(n):
    # Unit test that Wbas does what it should.
    # Verified up to n = 5.
    SBTn = SBT(n)
    Wbn = Wbas(QQ, n)
    Ebn = EUVbas(QQ, n)
    for (posUV, (lam, U, V)) in enumerate(SBTn):
        WUV = Wbn[posUV]
        for (posPQ, (mu, P, Q)) in enumerate(SBTn):
            WUVEPQ = Ebn[posPQ] * WUV
            if V != P and WUVEPQ != 0:
                print(U, V, P, Q)
                return False
            if V == P and WUVEPQ != Ebn[SBTn.index((lam, U, Q))]:
                print(U, V, P, Q)
                return False
    return True

r"""
Some old stuff specific to n = 3:

QS3 = QS(QQ, 3)
id3 = QS3(Permutation([1,2,3]))
s1 = QS3(Permutation([2,1,3]))
s2 = QS3(Permutation([1,3,2]))
t13 = QS3(Permutation([3,2,1]))
cyc1 = QS3(Permutation([2,3,1]))
cyc2 = QS3(Permutation([3,1,2]))
AWbas3 = [QQ(1)/6 * (id3 + s1 + s2 + t13 + cyc1 + cyc2),
          QQ(1)/3 * (id3 + s1 - t13 - cyc1),
          QQ(1)/3 * (-s1 + s2 + cyc1 - cyc2),
          QQ(1)/3 * (s2 - t13 - cyc1 + cyc2),
          QQ(1)/3 * (id3 - s1 + t13 - cyc2),
          QQ(1)/6 * (id3 - s1 - s2 - t13 + cyc1 + cyc2)]

# N,x,y,z,w,Nm = AWbas3

Murbas3 = [QS3.sum(QS3(w) for w in Permutations(3)),
           id3 + s1,
           s2 + cyc2,
           s2 + cyc1,
           id3 + t13,
           id3]
"""

@cached_function
def MurMat(field, n):
    # Matrix whose entries are the coordinates of the
    # elements of the row Murphy basis.
    MB = MurRowbas(field, n)
    return Matrix(field, [[v.coefficient(w) for w in Permutations(n)]
                          for v in MB])

# Munthe-Kaas-Stava basis.
# This is the basis of OSBB-products of the multilinear
# component of the tensor algebra, as in §5 of
# arXiv:2306.15582v3.
def MKS(field, w):
    # The Munthe-Kaas-Stava basis element corresponding
    # to the permutation w.
    n = len(w)
    SGA = QS(field, n)
    blocks = []
    block = []
    for (i, wi) in enumerate(w):
        if (not block) or wi < block[-1]:
            block.append(wi)
        else:
            last = block.pop()
            blocks.append((block, 0))
            blocks.append(([last, wi], 1))
            block = []
    if n > 0:
        blocks.append((block, 0))
    # print(blocks)
    # The list "blocks" now consists of pairs
    # (list, 0), which are to be symmetrized,
    # and (2-element-list, 1), which are to be
    # antisymmetrized.
    res = SGA.zero()
    for w in product(*[[perm for perm in permutations(block[0])] for block in blocks]):
        #print(w)
        perm = sum([block for block in w], tuple())
        num_of_flips = sum([1 for block in blocks
                            if block[1] == 1 and
                            tuple(block[0]) not in w])
        #print(perm, num_of_flips)
        res += (-1) ** num_of_flips * SGA(Permutation(list(perm)))
    return res

@cached_function
def MKSbas(field, n):
    return [MKS(field, w) for w in Permutations(n)]

# LRM-basis (Grinberg, Vassilieva, FPSAC 2026):

def lrmcomp(w):
    # The composition with D-set LRM' w.
    n = len(w)
    cur = [w[0]]
    curmin = w[0]
    for wi in w[1:]:
        if wi < curmin:
            cur.append(wi)
            curmin = wi
    I = sorted(cur) + [n+1]
    return Composition([Ii1 - Ii for (Ii1, Ii) in zip(I[1:], I)])

def LRMvec(field, w):
    # The element of the LRM-basis corresponding
    # to the permutation ``w``.
    n = len(w)
    B = DescentAlgebra(field, n).B()
    return w * B[lrmcomp(w)].to_symmetric_group_algebra()

@cached_function
def LRMbas(field, n):
    return [LRMvec(field, w) for w in Permutations(n)]

