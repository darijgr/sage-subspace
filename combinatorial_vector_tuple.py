# %attach /cygdrive/d/math/TEMPrepth/sage-subspace/combinatorial_vector_tuple.py

from sage.misc.cachefunc import cached_method
import operator

class VectorTuple():
    r"""
    Vector tuple.

    This is a list of vectors in some given module
    (not necessarily a concrete
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
        sage: a = VectorTuple([x, y, z, w], ambient=E)
        sage: a
        Vector tuple [x, y, z, w] in Rational Field-module
        The exterior algebra of rank 4 over Rational Field
        sage: a.list()
        [x, y, z, w]
        sage: a.dimension()
        4
        sage: a.echelon()
        Vector tuple [w, z, y, x] in Rational Field-module
        The exterior algebra of rank 4 over Rational Field
        sage: a.coefficients(2*x + 3*y + 7*w)
        [2, 3, 0, 7]
        sage: a.coefficients(x*y) # None, because x*y is not in span(a).
        sage: a.echelon().coefficients(2*x + 3*y + 7*w)
        [7, 0, 3, 2]
        sage: aa = a.product(a); aa
        Vector tuple
        [0, x^y, x^z, x^w, -x^y, 0, y^z, y^w, -x^z, -y^z, 0, z^w, -x^w, -y^w, -z^w, 0]
        in Rational Field-module The exterior algebra of rank 4 over Rational Field
        sage: aa.echelon()
        Vector tuple [z^w, y^w, y^z, x^w, x^z, x^y] in Rational Field-module
        The exterior algebra of rank 4 over Rational Field
        sage: [a.power(i).dimension() for i in range(6)]
        [1, 4, 6, 4, 1, 0]
        sage: a.power(3).echelon()
        Vector tuple [y^z^w, x^z^w, x^y^w, x^y^z] in Rational Field-module
        The exterior algebra of rank 4 over Rational Field
        sage: aa.span_contains(x*y - y*x)
        True
        sage: aa.span_contains(x*y - y*x + z)
        False
        sage: aca = a.commutator(a).echelon(); aca
        Vector tuple [2*z^w, 2*y^w, 2*y^z, 2*x^w, 2*x^z, 2*x^y] in
        Rational Field-module The exterior algebra of rank 4 over Rational Field
        sage: aa.span_contains_as_subset(aca)
        True
        sage: aca.span_contains_as_subset(aa)
        True
        sage: aa.span_equals(aca)
        True
        sage: aca.monicized()
        Vector tuple [z^w, y^w, y^z, x^w, x^z, x^y] in Rational Field-module
        The exterior algebra of rank 4 over Rational Field

    Now, we do the same over `GF(2)`::

        sage: E = ExteriorAlgebra(GF(2), ["x","y","z","w"])
        sage: x,y,z,w = E.gens()
        sage: a = VectorTuple([x, y, z, w], ambient=E)
        sage: a
        Vector tuple [x, y, z, w] in Finite Field of size 2-module
        The exterior algebra of rank 4 over Finite Field of size 2
        sage: a.list()
        [x, y, z, w]
        sage: a.dimension()
        4
        sage: a.echelon()
        Vector tuple [w, z, y, x] in Finite Field of size 2-module
        The exterior algebra of rank 4 over Finite Field of size 2
        sage: a.coefficients(2*x + 3*y + 7*w)
        [0, 1, 0, 1]
        sage: a.echelon().coefficients(2*x + 3*y + 7*w)
        [1, 0, 1, 0]
        sage: aa = a.product(a); aa
        Vector tuple
        [0, x^y, x^z, x^w, x^y, 0, y^z, y^w, x^z, y^z, 0, z^w, x^w, y^w, z^w, 0]
        in Finite Field of size 2-module The exterior algebra of rank 4
        over Finite Field of size 2
        sage: aa.echelon()
        Vector tuple [z^w, y^w, y^z, x^w, x^z, x^y]
        in Finite Field of size 2-module The exterior algebra of rank 4
        over Finite Field of size 2
        sage: [a.power(i).dimension() for i in range(6)]
        [1, 4, 6, 4, 1, 0]
        sage: a.power(3).echelon()
        Vector tuple [y^z^w, x^z^w, x^y^w, x^y^z]
        in Finite Field of size 2-module The exterior algebra of rank 4
        over Finite Field of size 2
        sage: aa.span_contains(x*y - y*x)
        True
        sage: aa.span_contains(x*y - y*x + z)
        False
        sage: aca = a.commutator(a).echelon(); aca
        Vector tuple []
        in Finite Field of size 2-module The exterior algebra of rank 4
        over Finite Field of size 2
        sage: aa.span_contains_as_subset(aca)
        True
        sage: aca.span_contains_as_subset(aa)
        False
        sage: aa.span_equals(aca)
        False

    Echelon form vs. reduced echelon form::

        sage: a = VectorTuple([x+y,x,y])
        sage: a.echelon_reduced()
        Vector tuple [y, x]
        in Finite Field of size 2-module The exterior algebra of rank 4 over Finite Field of size 2
        sage: a.echelon()
        Vector tuple [x + y, x]
        in Finite Field of size 2-module The exterior algebra of rank 4 over Finite Field of size 2

    For our next example, we shall work with polynomials.
    Since ``PolynomialRing`` currently does not yield a
    ``CombinatorialFreeModule``, we cannot use ``VectorTuple``
    with elements of the former.
    Instead, we implement polynomials as elements of the
    monoid algebra of a free abelian monoid
    (we set the ``prefix`` and ``bracket`` arguments
    merely to make the output less clumsy)::

        sage: Mon = FreeAbelianMonoid(["a", "b", "c", "d"])
        sage: P = Mon.algebra(QQ, prefix="", bracket=False)
        sage: a, b, c, d = P.gens()
        sage: a*b + b*a
        2*F['a']*F['b']

    Apart from the awkward output and the lack of
    polynomial-specific functionality, this is a perfectly
    fine implementation of the polynomial ring in four
    variables `a, b, c, d` over `\QQ`.
    The order on monomials is the lexicographic order
    on their sorted representations (i.e., to compare
    two monomials, we write them as words of the form
    `aa\ldots a bb\ldots b cc\ldots c dd\ldots d`, and
    compare these words lexicographically, where
    `a < b < c < d`). Beware that this is not a monomial
    order in the sense of commutative algebra! ::

        sage: a < a*a
        True
        sage: a*b < a*a*b
        False

    (See
    :meth:`~sage.monoids.indexed_free_monoid.IndexedMonoidElement._richcmp_`
    for the definition of the order.)

    We now define vector tuples and reduce some
    polynomials modulo them::

        sage: xs = VectorTuple([a*b + b*a, a*c + c*a, a*d + d*a])
        sage: xs
        Vector tuple
        [2*F['a']*F['b'], 2*F['a']*F['c'], 2*F['a']*F['d']]
        in Rational Field-module Algebra of Free abelian monoid
        indexed by {'a', 'b', 'c', 'd'} over Rational Field
        sage: h = (a+b+c+d)**2
        sage: xs.reduction(h)
        F['a']^2 + F['b']^2 + 2*F['b']*F['c'] + 2*F['b']*F['d']
         + F['c']^2 + 2*F['c']*F['d'] + F['d']^2
        sage: xs.coefficients(a*(b-c))
        [1/2, -1/2, 0]
        sage: xs.coefficients(b*c)
        sage: xs.coefficients(b*c+a*d)
        sage: xs.coefficients(P.zero())
        [0, 0, 0]
        sage: xs.quo_rem(b*c+a*d)
        ([0, 0, 1/2], F['b']*F['c'])

        sage: xs = VectorTuple([a*b - b*c, b*c - b*d, a*d])
        sage: xs
        Vector tuple
        [F['a']*F['b'] - F['b']*F['c'],
         F['b']*F['c'] - F['b']*F['d'], F['a']*F['d']]
        in Rational Field-module Algebra of Free abelian monoid
        indexed by {'a', 'b', 'c', 'd'} over Rational Field
        sage: h = (a+b+c+d)**2
        sage: xs.reduction(h)
        F['a']^2 + 6*F['a']*F['b'] + 2*F['a']*F['c']
         + F['b']^2 + F['c']^2 + 2*F['c']*F['d'] + F['d']^2

        sage: xs = VectorTuple([d**2, (c+d)**2, (b+c+d)**2])
        sage: xs
        Vector tuple
        [F['d']^2, F['c']^2 + 2*F['c']*F['d'] + F['d']^2,
         F['b']^2 + 2*F['b']*F['c'] + 2*F['b']*F['d']
          + F['c']^2 + 2*F['c']*F['d'] + F['d']^2] 
        in Rational Field-module Algebra of Free abelian monoid
        indexed by {'a', 'b', 'c', 'd'} over Rational Field
        sage: h = (a+b+c+d)**2 + (c-d)**2
        sage: xs.reduction(h)
        F['a']^2 + 2*F['a']*F['b'] + 2*F['a']*F['c']
         + 2*F['a']*F['d'] + 2*F['c']^2
        sage: xs.quo_rem(b**2 + 2*b*c)
        ([0, 0, 0], F['b']^2 + 2*F['b']*F['c'])
        sage: xs.quo_rem(b * (b+2*c+2*d))
        ([0, -1, 1], 0)
        sage: xs.quo_rem(b * (b+c+d))
        ([0, -1/2, 1/2], 1/2*F['b']^2)

    A few other methods::

        sage: xs = VectorTuple([a, a + 2*b, b + 3*c])
        sage: xs.monicized()
        Vector tuple
        [F['a'], 1/2*F['a'] + F['b'], 1/3*F['b'] + F['c']]
        in Rational Field-module Algebra of Free abelian
        monoid indexed by {'a', 'b', 'c', 'd'} over Rational Field
        sage: xs.echelon().monicized().list()
        [F['c'], F['b'], F['a']]

    TODO: Polynomial ring, through free abelian monoid.
    
    TODO: Symmetric group algebra.
    
    TODO: Free algebra, e.g., check Dynkin idempotent.
    """
    
    def __init__(self, xs, ambient=None):
        r"""
        Create a vector tuple.

        Syntax: ``VectorTuple(xs, ambient=M)``, where
        ``R`` is the base ring, ``M`` is the module
        containing the vectors, and ``xs`` is the list
        of the vectors.
        If ``ambient`` is not provided, then
        ``parent(xs[0])`` is being used as default.
        """
        self._vectors = xs[:]
        if ambient is None:
            ambient = parent(xs[0])
            if not all(ambient.is_parent_of(x) for x in xs):
                raise ValueError("ambient does not contain all the xs elements")
        self._basering = ambient.base_ring()
        self._ambient = ambient
        return object.__init__(self)
    
    def __getitem__(self, i):
        """
        Return the `i`-th vector `v_i`.
        """
        return self._vectors[i]

    def __len__(self):
        """
        Return the number `n` of vectors in ``self``.
        """
        return len(self._vectors)

    def __iter__(self):
        """
        Iterate over the vectors in ``self``.
        """
        return iter(self._vectors[:])
        # Do I need [:] here?

    def ambient(self):
        """
        Return the ambient module `M` of ``self``.
        """
        return self._ambient

    def base_ring(self):
        r"""
        Return the base ring `R` of ``self``.
        """
        return self._basering
    
    def list(self, copy=True):
        r"""
        Return the list of vectors in ``self``.

        If ``copy`` is set to ``False``, the original
        list is returned, otherwise a copy.
        """
        if copy:
            return self._vectors[:]
        return self._vectors

    def __repr__(self):
        return "Vector tuple " + str(self._vectors) + " in " + str(self._basering) + "-module " + str(self._ambient)

    def show(self):
        print('Vector tuple in module:')
        print self._ambient
        print('over the base ring')
        print self._basering
        print('The vectors are:')
        for i in self._vectors: print i

    def span(self):
        """
        Return the span of ``self``, as a
        :class:`AbstractSubspace`.
        """
        raise NotImplementedError # TODO
    
    def reduction_blind(self, v):
        r"""
        Compute ``self.reduction(v)``, under the assumption
        that ``self`` is already in echelon form.
        It is assumed that the module knows division by elements of
        the base ring and equality checking.
        """
        R = self._basering
        w = v
        xs = self._vectors # the vectors in ``self`` in echelon form,
                           # with highest terms strictly decreasing.
        for t in xs:
            wcoeffs = w.monomial_coefficients()
            tleader = max(t.monomial_coefficients().items())
            if tleader[0] in wcoeffs:
                w = w - R(wcoeffs[tleader[0]]) / R(tleader[1]) * t
        return w
    
    @cached_method
    def echelon(self, reduced=True, monic=True):
        r"""
        Return an echelon form of ``self``. This is a vector tuple
        which spans the same submodule as ``self``, but has its
        sequence of highest terms strictly decreasing.
        It is assumed that the module knows division by elements of
        the base ring and equality checking.
        """
        R = self._basering
        echeloned = [] # The echeloned variable contains a list of
                       # vectors already brought into echelon form.
                       # This list will grow step by step until it
                       # spans the same submodule as self.
        for v in self._vectors:
            w = v
            # Now reduce v modulo the span of the vectors in echeloned:
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
        res = VectorTuple(echeloned, ambient=self._ambient)
        res.echelon.set_cache(res)
        return res

    def echelon2(self):
        r"""
        Different implementation of :meth:`echelon`.
        Returns the same result (or, rather, a vector tuple
        that represents the same list of vectors).
        """
        R = self._basering
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
        for v in self._vectors:
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
        res = VectorTuple([t for i, t in reversed(echeloned)], ambient=self._ambient)
        res.echelon.set_cache(res)
        return res

    @cached_method
    def echelon_reduced(self):
        r"""
        Return the reduced echelon form of ``self``. This is an
        echelon form of ``self`` such that the highest term of
        any of its vectors occurs in no other of its vectors.
        It is assumed that the module knows division by elements of
        the base ring and equality checking.
        """
        R = self._basering
        ech = self.echelon()._vectors[:]
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
        res = VectorTuple(ech, ambient=self._ambient)
        res.echelon.set_cache(res)
        res.echelon_reduced.set_cache(res)
        return res

    def monicized(self):
        r"""
        Return the vector list obtained from ``self`` by
        dividing each nonzero vector in ``self`` by its leading
        coefficient.
        """
        M = self._ambient
        R = self._basering
        mon = [(x if x == 0
                else ~(R(max(x.monomial_coefficients().items())[1])) * x)
               for x in self._vectors]
        return VectorTuple(mon, ambient=M)

    @cached_method
    def echelon_reducer(self):
        r"""
        Return the matrix whose `i`-th row contains the coefficients
        of the `i`-th vector of ``self.echelon()`` with respect to the
        spanning set ``self``.
        The matrix is returned as a list of lists (the inner lists
        being its rows).
        """
        from itertools import izip
        echeloned = [] # The echeloned variable contains a list of
                       # pairs (v, hs). The v components are vectors
                       # already brought into echelon form. The hs
                       # component of each such pair contains the
                       # coefficients of v with respect to self.
                       # This list will grow step by step until its
                       # v components span the same submodule as
                       # self.
        R = self._basering
        zero = R.zero()
        one = R.one()
        l = len(self._vectors)
        for (i, v) in enumerate(self._vectors):
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
        r"""
        Return a dictionary whose items `i: xs` stand for the
        syzygies of ``self``. More specifically, an item `i: xs` means that
        the `i`-th vector in ``self`` equals the linear combination
        of the first `i-1` vectors with coefficients taken from
        `xs`. (Note that `xs` is a length-`i-1` list.)
        """
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
        l = len(self._vectors)
        R = self._basering
        zero = R.zero()
        one = R.one()
        for (i, v) in enumerate(self._vectors):
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
        r"""
        Given a vector `v \in V`, return a list
        `[a_0, a_1, \ldots, a_{n-1}]` of scalars in `R`
        such that
        `v = a_0 v_0 + a_1 v_1 + \cdots + a_{n-1} v_{n-1}`
        if such a list exists; otherwise return ``None``.
        
        This list may not be unique; it is, however, computed
        deterministically from ``self``.
        """
        ech = self.echelon()
        l = len(ech._vectors)
        ring = self._basering
        blinds = ech.coefficients_blind(v)
        k = len(self._vectors)
        red = self.echelon_reducer()
        res = [ring.sum((blinds[i] * red[i][j]
                         for i in range(l)))
               for j in range(k)]
        if self._ambient.sum(res[i] * self._vectors[i] for i in range(len(self._vectors))) == v:
            return res

    def coefficients_blind(self, v):
        r"""
        Compute ``self.coefficients(v)``, under the assumption
        that ``self`` is in echelon form, and that ``v`` is a
        linear combination of the vectors of ``self``.
        """
        R = self._basering
        w = v
        coeffs = [R.zero()] * len(self._vectors)
        for (i, t) in enumerate(self._vectors):
            wcoeffs = w.monomial_coefficients()
            tleader = max(t.monomial_coefficients().items())
            if tleader[0] in wcoeffs:
                coeff = R(wcoeffs[tleader[0]]) / R(tleader[1])
                w = w - coeff * t
                coeffs[i] = coeff
        return coeffs

    def reduction(self, v):
        r"""
        Compute the reduction of a vector ``v`` modulo ``self``.

        This is the unique vector ``w`` satisfying
        `v \equiv w \mod V` that contains none of the
        leading monomials of the vectors in ``self.echelon()``.
        """
        return self.echelon().reduction_blind(v)

    def quo_rem(self, v):
        r"""
        Given a vector `v \in V`, return a pair
        `([a_0, a_1, \ldots, a_{n-1}], w)`, where
        `[a_0, a_1, \ldots, a_{n-1}]` is a list of scalars in
        `R`, and where `w \in V` is a vector such that
        `v = a_0 v_0 + a_1 v_1 + \cdots + a_{n-1} v_{n-1} + w`
        and such that `w` contains none of the leading
        monomials of the vectors in ``self.echelon()``.
        Note that `w = 0` if `v` lies in the span of ``self``,
        and that changing `v` by an element of `V` must leave
        `w` unchanged.
        """
        red = self.reduction(v)
        return (self.coefficients(v - red), red)
        # TODO: Optimize! This is doing a lot of unnecessary work.

    def span_contains_blind(self, v):
        r"""
        Compute ``self.span_contains(v)``, under the assumption that
        ``self`` is already in echelon form.
        """
        return (self.reduction_blind(v) == 0)

    def span_contains(self, v):
        r"""
        Check whether a vector `v \in V` lies in the span
        of ``self``.
        """
        return (self.echelon().reduction_blind(v) == 0)
    
    def span_contains_as_subset(self, anotherlist, verbose=False):
        r"""
        Check whether the span of ``self`` contains the
        span of a further vector tuple ``anotherlist`` as
        a subspace.
        """
        xs = self.echelon()
        for y in anotherlist.list(copy=False):
            if not xs.span_contains_blind(y):
                if verbose == True:
                    print "The offending vector is: "
                    print y
                return False
                break
        return True
    
    def span_is_contained_in(self, anotherlist, verbose=False):
        r"""
        Check whether the span of ``self`` is contained
        in the span of a further vector tuple ``anotherlist``.
        """
        return anotherlist.span_contains_as_subset(self, verbose=verbose)
    
    def span_equals(self, anotherlist, verbose=False):
        r"""
        Check whether the span of ``self`` equals
        the span of a further vector tuple ``anotherlist``.
        """
        return (anotherlist.span_contains_as_subset(self, verbose=verbose)
                and self.span_contains_as_subset(anotherlist, verbose=verbose))
    
    def concatenate(self, anotherlist):
        r"""
        Return the list of vectors obtained by
        concatenating the vector lists ``self`` and ``anotherlist``.

        This assumes that the vector lists ``self`` and ``anotherlist``
        have the same base ring and the same ambient
        module.
        """
        us = self._vectors[:]
        us.extend(anotherlist.list(copy=False))
        return VectorTuple(us, ambient=self._ambient)
    
    def intersection_blind(self, anotherlist):
        r"""
        Compute ``self.intersection(anotherlist)``, under
        the assumption that ``self`` and ``anotherlist``
        are already in echelon form.
        """
        vs = self._vectors
        ws = anotherlist.list(copy=False)
        M = self._ambient
        n = len(vs)
        vsws = vs + ws
        syzzies = VectorTuple(vsws, ambient=M).syzygies()
        reslist = []
        for (k, syz) in syzzies.iteritems():
            # Throw the vs-part of syz away.
            k2 = k - n
            syz2 = syz[n:]
            reslist.append(ws[k2] - M.sum(syz2[i] * ws[i] for i in range(k2)))
        return VectorTuple(reslist, ambient=M)

    def intersection(self, anotherlist):
        r"""
        Given a further vector tuple ``anotherlist``
        (whose ambient module is ``self.ambient``),
        return a vector tuple that spans the intersection
        of the span of ``self`` with the span of
        ``anotherlist``.

        This implementation returns an echelonized vector
        tuple.
        """
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
        us = self._vectors
        vs = anotherlist.list(copy=False)
        ws = [op(p, q) for p in us for q in vs]
        return VectorTuple(ws, ambient=self._ambient)
        
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
        us = self._vectors
        vs = anotherlist.list(copy=False)
        ws = [op(p, q) - op(q, p) for p in us for q in vs]
        return VectorTuple(ws, ambient=self._ambient)
        
    def power(self, n):
        # Returns the n-th power of the list with respect to the
        # above-defined product function.
        # This assumes that multiplication is actual multiplication.
        if n == 0:
            return VectorTuple([self._ambient.one()], ambient=self._ambient)
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
        ys = [f(x) for x in self._vectors]
        return VectorTuple(ys, ambient=f.codomain())

    def kernel(self, f):
        # Returns a basis of the kernel of a morphism f (restricted
        # to the span of self).
        img = self.image(f)
        syzzies = img.syzygies()
        vects = self._vectors
        M = self._ambient
        ker = []
        for i, xs in syzzies.iteritems():
            ker.append(vects[i] - M.sum((c * vects[j] for (j, c) in enumerate(xs))))
        return VectorTuple(ker, ambient=M)

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
    A = [0]*n
    l = len(U)
    if l > n: l = n    # forget about the generators that are too high for us to use
    for i in range(l):
        A[i] = U[i]
    for i in range(l, n):
        A[i] = VectorTuple([], ambient=A)
    for i in range(n):
        for j in range(1,i+1):
            if i - j < l:
                A[i] = A[i].concatenate(A[j-1].product(U[i-j])).echelon()
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
    B = subalgcomp(A, U, n)
    A = [0]*n
    l = len(V)
    if l > n: l = n    # forget about the generators that are too high for us to use
    for i in range(l):
        A[i] = V[i]
    for i in range(l, n):
        A[i] = VectorTuple([], ambient=A)
    for i in range(n):
        for j in range(1,i+1):
            A[i] = A[i].concatenate(B[j-1].product(A[i-j]).concatenate(A[i-j].product(B[j-1]))).echelon()
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

Epk_span = VectorTuple([E(xs) for xs in lacsubs], ambient=QSn)
"""

