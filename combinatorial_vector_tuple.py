# %attach /cygdrive/d/math/TEMPrepth/sage-subspace/combinatorial_vector_tuple.py

class CombinatorialVectorTuple():
    r"""
    List of vectors in some given module (not necessarily a concrete
    vector space, but something that implements ``monomial_coefficients``
    like a combinatorial free module does).
    Certain methods, however, will require the ground
    ring to be a field, or at least certain elements of it to be
    invertible.

    This probably implements some functionality from
    https://trac.sagemath.org/ticket/11111 , but it does so
    in a rather clumsy way.

    EXAMPLES:

    Let us consider some lists of vectors in the exterior algebra
    of a vector space::

        sage: E = ExteriorAlgebra(QQ, ["x","y","z","w"])
        sage: x,y,z,w = E.gens()
        sage: a = CombinatorialVectorTuple(QQ, E, [x, y, z, w])
        sage: a
        List of vectors [x, y, z, w] in Rational Field-module
        The exterior algebra of rank 4 over Rational Field
        sage: a.vectors
        [x, y, z, w]
        sage: a.dimension()
        4
        sage: a.echelon()
        List of vectors [w, z, y, x] in Rational Field-module
        The exterior algebra of rank 4 over Rational Field
        sage: a.coefficients(2*x + 3*y + 7*w)
        [2, 3, 0, 7]
        sage: a.coefficients(x*y) # None, because x*y is not in span(a).
        sage: a.echelon().coefficients(2*x + 3*y + 7*w)
        [7, 0, 3, 2]
        sage: aa = a.product(a); aa
        List of vectors
        [0, x^y, x^z, x^w, -x^y, 0, y^z, y^w, -x^z, -y^z, 0, z^w, -x^w, -y^w, -z^w, 0]
        in Rational Field-module The exterior algebra of rank 4 over Rational Field
        sage: aa.echelon()
        List of vectors [z^w, y^w, y^z, x^w, x^z, x^y] in Rational Field-module
        The exterior algebra of rank 4 over Rational Field
        sage: [a.power(i).dimension() for i in range(6)]
        [1, 4, 6, 4, 1, 0]
        sage: a.power(3).echelon()
        List of vectors [y^z^w, x^z^w, x^y^w, x^y^z] in Rational Field-module
        The exterior algebra of rank 4 over Rational Field
        sage: aa.contains(x*y - y*x)
        True
        sage: aa.contains(x*y - y*x + z)
        False
        sage: aca = a.commutator(a).echelon(); aca
        List of vectors [2*z^w, 2*y^w, 2*y^z, 2*x^w, 2*x^z, 2*x^y] in
        Rational Field-module The exterior algebra of rank 4 over Rational Field
        sage: aa.isbiggerthan(aca)
        True
        sage: aca.isbiggerthan(aa)
        True
        sage: aa.isequivalentto(aca)
        True

    Now, we do the same over `GF(2)`::

        sage: E = ExteriorAlgebra(GF(2), ["x","y","z","w"])
        sage: x,y,z,w = E.gens()
        sage: a = CombinatorialVectorTuple(GF(2), E, [x, y, z, w])
        sage: a
        List of vectors [x, y, z, w] in Finite Field of size 2-module
        The exterior algebra of rank 4 over Finite Field of size 2
        sage: a.vectors
        [x, y, z, w]
        sage: a.dimension()
        4
        sage: a.echelon()
        List of vectors [w, z, y, x] in Finite Field of size 2-module
        The exterior algebra of rank 4 over Finite Field of size 2
        sage: a.coefficients(2*x + 3*y + 7*w)
        [0, 1, 0, 1]
        sage: a.echelon().coefficients(2*x + 3*y + 7*w)
        [1, 0, 1, 0]
        sage: aa = a.product(a); aa
        List of vectors
        [0, x^y, x^z, x^w, x^y, 0, y^z, y^w, x^z, y^z, 0, z^w, x^w, y^w, z^w, 0]
        in Finite Field of size 2-module The exterior algebra of rank 4
        over Finite Field of size 2
        sage: aa.echelon()
        List of vectors [z^w, y^w, y^z, x^w, x^z, x^y]
        in Finite Field of size 2-module The exterior algebra of rank 4
        over Finite Field of size 2
        sage: [a.power(i).dimension() for i in range(6)]
        [1, 4, 6, 4, 1, 0]
        sage: a.power(3).echelon()
        List of vectors [y^z^w, x^z^w, x^y^w, x^y^z]
        in Finite Field of size 2-module The exterior algebra of rank 4
        over Finite Field of size 2
        sage: aa.contains(x*y - y*x)
        True
        sage: aa.contains(x*y - y*x + z)
        False
        sage: aca = a.commutator(a).echelon(); aca
        List of vectors []
        in Finite Field of size 2-module The exterior algebra of rank 4
        over Finite Field of size 2
        sage: aa.isbiggerthan(aca)
        True
        sage: aca.isbiggerthan(aa)
        False
        sage: aa.isequivalentto(aca)
        False

    Echelon form vs. reduced echelon form::

        sage: a = CombinatorialVectorTuple(GF(2), E, [x+y,x,y])
        sage: a.echelon_reduced()
        List of vectors [y, x]
        in Finite Field of size 2-module The exterior algebra of rank 4 over Finite Field of size 2
        sage: a.echelon()
        List of vectors [x + y, x]
        in Finite Field of size 2-module The exterior algebra of rank 4 over Finite Field of size 2

    TODO: Polynomial ring, through free abelian monoid.
    
    TODO: Symmetric group algebra.
    
    TODO: Free algebra, e.g., check Dynkin idempotent.
    """
    
    def __init__(self, basering, module, xs):
        # Create a list of vectors.
        # Syntax: ``CombinatorialVectorTuple(R, M, xs)``, where
        # ``R`` is the base ring, ``M`` is the module
        # containing the vectors, and ``xs`` is the list
        # of the vectors.
        # NOTE: Make sure ``xs`` is never mutated. The
        # safe way is to use
        # ``CombinatorialVectorTuple(R, M, xs.copy())``.
        # TODO: Should we use self.vectors = tuple(xs)?
        self.vectors = xs[:]
        self.basering = basering
        self.undermodule = module
        return object.__init__(self)
    
    def __repr__(self):
        return "List of vectors " + str(self.vectors) + " in " + str(self.basering) + "-module " + str(self.undermodule)

    def base_ring(self):
        return self.basering
    
    def underlying_module(self):
        return self.undermodule
    
    def list(self):
        return self.vectors[:]
    
    def show(self):
        print('List of vectors in module:')
        print self.undermodule
        print('over the base ring')
        print self.basering
        print('The vectors are:')
        for i in self.vectors: print i
    
    def reduction_blind(self, v):
        # Computes the reduction of a vector v modulo self, under the
        # assumption that self is in echelon form already.
        # It is assumed that the module knows division by elements of
        # the base ring and equality checking.
        R = self.basering
        w = v
        xs = self.vectors # the vectors in the list - in echelon form,
                          # with highest terms strictly decreasing.
        for t in xs:
            wcoeffs = w.monomial_coefficients()
            tleader = max(t.monomial_coefficients().items())
            if tleader[0] in wcoeffs:
                w = w - R(wcoeffs[tleader[0]]) / R(tleader[1]) * t
        return w
    
    @cached_method
    def echelon(self):
        # Returns an echelon form of self. This is a list of vectors
        # which spans the same submodule as self, but has its
        # sequence of highest terms strictly decreasing.
        # It is assumed that the module knows division by elements of
        # the base ring and equality checking.
        R = self.basering
        echeloned = [] # The echeloned variable contains a list of
                       # vectors already brought into echelon form.
                       # This list will grow step by step until it
                       # spans the same submodule as self.
        for v in self.vectors:
            w = v
            # Now reduce v modulo echeloned:
            for t in echeloned:
                wcoeffs = w.monomial_coefficients()
                tleader = max(t.monomial_coefficients().items())
                if tleader[0] in wcoeffs:
                    w = w - R(wcoeffs[tleader[0]]) / R(tleader[1]) * t
            # Now w is the reduction of v modulo echelon.
            # If w == 0, then v was linearly dependent on echelon,
            #      and we don't have to do anything.
            # If w != 0, then we now add w to echelon.
            if w != 0:
                echeloned.append(w)
                # w might have been added to the wrong place, so
                # let us sort. I know this is not the best way;
                # if the bisect module would allow for keys, then
                # this would be easy to improve.
                echeloned.sort(key=lambda x : max(x.monomial_coefficients().keys()), reverse=True)
        res = CombinatorialVectorTuple(self.basering, self.undermodule, echeloned)
        res.echelon.set_cache(res)
        return res

    def echelon2(self):
        # Another way to get an echelon form. Possibly sometimes
        # faster than echelon. It saves the result into the cache
        # of echelon.
        R = self.basering
        from bisect import insort
        echeloned = [] # The echeloned variable contains a list of
                       # pairs (i, v), with i being a monomial and v
                       # being a vector such that the highest monomial
                       # appearing in v with nonzero coefficient is i.
                       # This list will grow step by step. At every
                       # step k, the vectors v appearing in the list
                       # will be a basis of some submodule of the
                       # module spanned by self (namely, of the
                       # submodule spanned by the first k vectors of
                       # self). At the end, the vectors v appearing in
                       # the list will be a basis of the submodule
                       # spanned by self. At every step, the list will
                       # be sorted by increasing i, so the v's form a
                       # basis in row echelon form (whence the name
                       # "echelon").
        for v in self.vectors:
            w = v
            # Now reduce v modulo the span of the vectors in echeloned:
            for i, t in reversed(echeloned):
                wcoeffs = w.monomial_coefficients()
                if i in wcoeffs:
                    tleader = t.monomial_coefficients()[i]
                    w = w - R(wcoeffs[i]) / R(tleader) * t
            # Now w is the reduction of v modulo the span of the
            # v's in echeloned.
            # If w == 0, then v was linearly dependent on the v's in
            #      echeloned, and we don't have to do anything.
            # If w != 0, then we now add w to echeloned.
            if w != 0:
                insort(echeloned, (max(w.monomial_coefficients().keys()), w))
        # Now forget the i's in echeloned. This is probably not very
        # efficient.
        res = CombinatorialVectorTuple(self.basering, self.undermodule, [t for i, t in reversed(echeloned)])
        res.echelon.set_cache(res)
        return res

    @cached_method
    def echelon_reduced(self):
        # Returns the reduced echelon form of self. This is an
        # echelon form of self such that the highest term of
        # any of its vectors occurs in no other of its vectors.
        # It is assumed that the module knows division by elements of
        # the base ring and equality checking.
        R = self.basering
        ech = self.echelon().vectors[:]
        k = len(ech)
        for i in range(k-1, -1, -1):
            w = ech[i]
            # Now reduce v modulo ech[i+1:]:
            for t in ech[i+1:]:
                wcoeffs = w.monomial_coefficients()
                tleader = max(t.monomial_coefficients().items())
                if tleader[0] in wcoeffs:
                    w = w - R(wcoeffs[tleader[0]]) / R(tleader[1]) * t
            ech[i] = w
        res = CombinatorialVectorTuple(self.basering, self.undermodule, ech)
        res.echelon.set_cache(res)
        return res

    @cached_method
    def echelon_reducer(self):
        # Return the matrix whose i-th row contains the coefficients
        # of the i-th vector of self.echelon() with respect to the
        # spanning set self.
        # The matrix is returned as a list of lists (the inner lists
        # being its rows).
        from itertools import izip
        echeloned = [] # The echeloned variable contains a list of
                       # pairs (v, hs). The v components are vectors
                       # already brought into echelon form. The hs
                       # component of each such pair contains the
                       # coefficients of v with respect to self.
                       # This list will grow step by step until its
                       # v components span the same submodule as
                       # self.
        R = self.basering
        zero = R.zero()
        one = R.one()
        l = len(self.vectors)
        for (i, v) in enumerate(self.vectors):
            w = v
            whs = [zero] * l
            whs[i] = one
            # Now reduce v modulo echelon:
            for (t, hs) in echeloned:
                wcoeffs = w.monomial_coefficients()
                tleader = max(t.monomial_coefficients().items())
                if tleader[0] in wcoeffs:
                    coeff = R(wcoeffs[tleader[0]]) / R(tleader[1])
                    w = w - coeff * t
                    # Subtract coeff * hs from the vector whs:
                    whs = [whs_i - coeff * hs_i
                           for (whs_i, hs_i) in izip(whs, hs)]
            # Now w is the reduction of v modulo echelon.
            # If w == 0, then v was linearly dependent on echelon,
            #      and we don't have to do anything.
            # If w != 0, then we now add w to echelon.
            if w != 0:
                echeloned.append((w, whs))
                # w might have been added to the wrong place, so
                # let us sort. I know this is not the best way;
                # if the bisect module would allow for keys, then
                # this would be easy to improve.
                echeloned.sort(key=lambda x : max(x[0].monomial_coefficients().keys()), reverse=True)
        return [e[1] for e in echeloned]

    @cached_method
    def syzygies(self):
        # Returns a dictionary whose items `i: xs` stand for the
        # syzygies of self. More specifically, an item `i: xs` means that
        # the `i`-th vector in self equals the linear combination
        # of the first `i-1` vectors with coefficients taken from
        # `xs`. (Note that `xs` is a length-`i-1` list.)
        from itertools import izip
        echeloned = [] # The echeloned variable contains a list of
                       # pairs (v, hs). The v components are vectors
                       # already brought into echelon form. The hs
                       # component of each v contains the
                       # coefficients of v with respect to self.
                       # This list will grow step by step until its
                       # v components span the same submodule as
                       # self.
        syzzies = {}
        l = len(self.vectors)
        R = self.basering
        zero = R.zero()
        one = R.one()
        for (i, v) in enumerate(self.vectors):
            w = v
            whs = [zero] * l
            whs[i] = one
            # Now reduce v modulo echelon:
            for (t, hs) in echeloned:
                wcoeffs = w.monomial_coefficients()
                tleader = max(t.monomial_coefficients().items())
                if tleader[0] in wcoeffs:
                    coeff = R(wcoeffs[tleader[0]]) / R(tleader[1])
                    w = w - coeff * t
                    # Subtract coeff * hs from the vector whs:
                    whs = [whs_i - coeff * hs_i
                           for (whs_i, hs_i) in izip(whs, hs)]
            # Now w is the reduction of v modulo echelon.
            # If w == 0, then v was linearly dependent on echelon,
            #      and we record a syzygy.
            if w == 0:
                syzzies[i] = [-s for s in whs[:i]]
            # If w != 0, then we now add w to echelon.
            else:
                echeloned.append((w, whs))
                # w might have been added to the wrong place, so
                # let us sort. I know this is not the best way;
                # if the bisect module would allow for keys, then
                # this would be easy to improve.
                echeloned.sort(key=lambda x : max(x[0].monomial_coefficients().keys()), reverse=True)
        return syzzies

    def coefficients(self, v):
        # Computes a list [a_1, a_2, ..., a_k] of coefficients
        # such that v = a_1 v_1 + a_2 v_2 + ... + a_k v_k, where
        # self == [v_1, v_2, ..., v_k].
        # Return None (TODO: or raise ValueError)
        # if there is no such list.
        ech = self.echelon()
        l = len(ech.vectors)
        ring = self.basering
        blinds = ech.coefficients_blind(v)
        k = len(self.vectors)
        red = self.echelon_reducer()
        res = [ring.sum((blinds[i] * red[i][j]
                         for i in range(l)))
                for j in range(k)]
        if self.undermodule.sum(res[i] * self.vectors[i] for i in range(len(self.vectors))) == v:
            return res

    def coefficients_blind(self, v):
        # Computes a list [a_1, a_2, ..., a_k] of coefficients
        # such that v = a_1 v_1 + a_2 v_2 + ... + a_k v_k, where
        # self == [v_1, v_2, ..., v_k], under the assumption that
        # self is in echelon form, and that v is a linear
        # combination of the vectors of self.
        R = self.basering
        w = v
        coeffs = [R.zero()] * len(self.vectors)
        for (i, t) in enumerate(self.vectors):
            wcoeffs = w.monomial_coefficients()
            tleader = max(t.monomial_coefficients().items())
            if tleader[0] in wcoeffs:
                coeff = R(wcoeffs[tleader[0]]) / R(tleader[1])
                w = w - coeff * t
                coeffs[i] = coeff
        return coeffs

    def reduction(self, v):
        # Computes the reduction of a vector v modulo self.
        return self.echelon().reduction_blind(v)

    def contains_blind(self, v):
        # Finds out whether a vector v lies in the submodule generated
        # by self, under the assumption that self is in echelon form
        # already.
        return (self.reduction_blind(v) == 0)

    def contains(self, v):
        # Finds out whether a vector v lies in the submodule generated
        # by self.
        return (self.echelon().reduction_blind(v) == 0)
    
    def isbiggerthan(self, anotherlist, verbose=False):
        # Finds out whether the submodule spanned by self contains that
        # spanned by anotherlist (another list of vectors).
        xs = self.echelon()
        for y in anotherlist.list():
            if not xs.contains_blind(y):
                if verbose == True:
                    print "The offending vector is: "
                    print y
                return False
                break
        return True
    
    def issmallerthan(self, anotherlist, verbose=False):
        # Finds out whether the submodule spanned by self is contained
        # in that spanned by anotherlist (another list of vectors).
        return anotherlist.isbiggerthan(self, verbose=verbose)
    
    def isequivalentto(self, anotherlist, verbose=False):
        # Finds out whether the submodule spanned by self equals
        # that spanned by anotherlist (another list of vectors).
        # (For instance, self.isequivalentto(self.echelon()) should
        # always return True.)
        return (anotherlist.isbiggerthan(self, verbose=verbose)
                and self.isbiggerthan(anotherlist, verbose=verbose))
    
    def add(self, anotherlist):
        # Gives the disjoint union of self with anotherlist (another list
        # of vectors, which of course should be in the same module).
        # This union spans the sum of the respective submodules.
        us = self.vectors[:]
        us.extend(anotherlist.list())
        return CombinatorialVectorTuple(self.basering, self.undermodule, us)
    
    def intersection_blind(self, anotherlist):
        # Gives the intersection of the span of self with the span
        # of anotherlist (another list of vectors, which of course
        # should be in the same module), as an echelonized list of
        # vectors.
        # This assumes that self and anotherlist are in echelon form
        # already.
        vs = self.vectors
        ws = anotherlist.list()
        R = self.basering
        M = self.undermodule
        n = len(vs)
        m = len(ws)
        vsws = vs + ws
        syzzies = CombinatorialVectorTuple(R, M, vsws).syzygies()
        reslist = []
        for (k, syz) in syzzies.iteritems():
            # Throw the vs-part of syz away.
            k2 = k - n
            syz2 = syz[n:]
            reslist.append(ws[k2] - M.sum(syz2[i] * ws[i] for i in range(k2)))
        return CombinatorialVectorTuple(R, M, reslist)

    def intersection(self, anotherlist):
        # Gives the intersection of the span of self with the span
        # of anotherlist (another list of vectors, which of course
        # should be in the same module), as an echelonized list of
        # vectors.
        return self.echelon().intersection_blind(anotherlist.echelon())

    def product(self, anotherlist, op=operator.mul):
        # Gives the list formed by pairwise products of vectors in
        # ``self`` with vectors in ``anotherlist``.
        # Here, ``op`` is required to be a binary operation from
        # `M \times M` to `M`, where `M` is the underlying module
        # of ``self``. "Product" is understood to mean "image under
        # ``op``".
        # (By default, ``op`` is standard multiplication, which
        # assumes that `M` is an algebra.)
        # If ``op`` is bilinear, then the list returned by this
        # method spans the product of the respective submodules
        # (in the sense in which, e. g., the product of ideals is
        # defined).
        us = self.vectors
        vs = anotherlist.list()
        ws = [op(p, q) for p in us for q in vs]
        return CombinatorialVectorTuple(self.basering, self.undermodule, ws)
        
    def commutator(self, anotherlist, op=operator.mul):
        # Gives the list formed by pairwise commutators of vectors in
        # ``self`` with vectors in ``anotherlist``.
        # Here, ``op`` is required to be a binary operation from
        # `M \times M` to `M`, where `M` is the underlying module
        # of ``self``. "Commutator" is understood to mean "commutator
        # with respect to ``op``" (that is, the commutator of `a`
        # and `b` is defined as `op(a, b) - op(b, a)`).
        # (By default, ``op`` is standard multiplication, which
        # assumes that `M` is an algebra.)
        # If ``op`` is bilinear, then the list returned by this
        # method spans the commutator of the respective submodules.
        us = self.vectors
        vs = anotherlist.list()
        ws = [op(p, q) - op(q, p) for p in us for q in vs]
        return CombinatorialVectorTuple(self.basering, self.undermodule, ws)
        
    def power(self, n):
        # Returns the n-th power of the list with respect to the
        # above-defined product function.
        # This assumes that multiplication is actual multiplication.
        if n == 0:
            return CombinatorialVectorTuple(self.basering, self.undermodule, [self.undermodule.one()])
        elif n == 1:
            return self
        else:
            m = int(n) / int(2)
            M = n - m
            return self.power(m).product(self.power(M))
    
    def dimension(self):
        # Gives the dimension of the submodule generated by self.
        return len(self.echelon().list())
    
    def image(self, f):
        # Returns the list of the images of the vectors under a morphism f.
        # The module in which they lie is the codomain of f.
        ys = [f(x) for x in self.vectors]
        return CombinatorialVectorTuple(self.basering, f.codomain(), ys)

    def kernel(self, f):
        # Returns a basis of the kernel of a morphism f (restricted
        # to the span of self).
        img = self.image(f)
        syzzies = img.syzygies()
        vects = self.vectors
        R = self.basering
        M = self.undermodule
        ker = []
        for i, xs in syzzies.iteritems():
            ker.append(vects[i] - M.sum((c * vects[j] for (j, c) in enumerate(xs))))
        return CombinatorialVectorTuple(R, M, ker)

def subalgcomp(A, U, n):
    # INPUT:
    # A: a graded algebra over a base ring.
    # U: a list of listsofvectors, where each vector in U[i] lies in the
    # (i+1)-th graded component of A.
    # n: an integer.
    # OUTPUT:
    # a list of listsofvectors, the i-th of which spans the (i+1)-th
    # graded component of the subalgebra of A generated by U. The list
    # has length n.
    # (As should be clear from the above description, generators in the
    # 0-th graded component are not supported. This is because there is
    # no algorithm for finding the subalgebra generated by them in the
    # general case.)
    BR = A.base_ring()
    A = [0]*n
    l = len(U)
    if l > n: l = n    # forget about the generators that are too high for us to use
    for i in range(l):
        A[i] = U[i]
    for i in range(l, n):
        A[i] = CombinatorialVectorTuple(BR, A, [])
    for i in range(n):
        for j in range(1,i+1):
            if i - j < l:
                A[i] = A[i].add(A[j-1].product(U[i-j])).echelon()
    return A

def gradedideal(A, U, V, n):
    # INPUT:
    # A: a graded algebra over a base ring.
    # U: a list of listsofvectors, where each vector in U[i] lies in the
    # (i+1)-th graded component of A.
    # V: a list of listsofvectors, where each vector in V[i] lies in the
    # vector subspace of A generated by U[i].
    # n: an integer.
    # OUTPUT:
    # a list of listsofvectors, the i-th of which spans the (i+1)-th
    # graded component of the two-sided ideal generated by V inside the
    # subalgebra of A generated by U. The list has length n.
    # (As should be clear from the above description, generators in the
    # 0-th graded component are not supported. This is because there is
    # no algorithm for finding the subalgebra generated by them in the
    # general case.)
    BR = A.base_ring()
    B = subalgcomp(A, U, n)
    A = [0]*n
    l = len(V)
    if l > n: l = n    # forget about the generators that are too high for us to use
    for i in range(l):
        A[i] = V[i]
    for i in range(l, n):
        A[i] = CombinatorialVectorTuple(BR, A, [])
    for i in range(n):
        for j in range(1,i+1):
            A[i] = A[i].add(B[j-1].product(A[i-j]).add(A[i-j].product(B[j-1]))).echelon()
    return A

r"""
DefaultBase = QQ

n = 6
QSn = SymmetricGroupAlgebra(QQ, n)

def Epk(pi):
    return tuple([i for i in range(n) if (i == 0 or pi[i] > pi[i-1]) and (i == n-1 or pi[i] > pi[i+1])])

def epk(pi):
    return len(Epk(pi))

def E(xs):
    return QSn.sum(QSn.basis()[pi] for pi in Permutations(n) if Epk(pi) == xs)

def is_lacunar(xs):
    return all(x+1 not in xs for x in xs)

from itertools import combinations

lacsubs = [xs for k in range(n+1) for xs in combinations(range(n), k) if is_lacunar(xs)]

lacsubs

Epk_span = CombinatorialVectorTuple(QQ, QSn, [E(xs) for xs in lacsubs])
"""

