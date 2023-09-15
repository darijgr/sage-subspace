"""
Generic accessors for vectors in SageMath

Eventually this API should be available as a collection of methods
implemented by all ModulesWithBasis in Sage.

The implementation below works for CombinatorialFreeModule, and should
be generalized progressively
"""

from sage.structure.parent cimport Parent
from sage.structure.element cimport Element
from sage.rings.polynomial.polydict cimport ETuple
from sage.groups.perm_gps.permgroup_element cimport PermutationGroupElement

import sage.combinat.tableau
from sage.combinat.words.word import Word

from sage.combinat.free_module import CombinatorialFreeModule

cpdef items_of_vector(Element v):
    """
    Return an iterator over the pairs ``(index, coefficient)`` for `v`.

    INPUT:

    - ``v`` -- an element of some vector space or free module

    EXAMPLES:

    This method handles indexed free module elements::

        sage: E = CombinatorialFreeModule(QQ, [1,2,4,8,16])
        sage: v = E.an_element(); v
        2*B[1] + 2*B[2] + 3*B[4]
        sage: list(items_of_vector(v))
        [(1, 2), (2, 2), (4, 3)]

    free module elements::

        sage: v = vector([4,0,1,2])
        sage: list(items_of_vector(v))
        [(0, 4), (2, 1), (3, 2)]

        sage: v = vector([4,0,1,2], sparse=True)
        sage: list(items_of_vector(v))
        [(0, 4), (2, 1), (3, 2)]

    multivariate polynomials::

        sage: P = QQ['x,y,z']
        sage: x,y,z = P.gens()
        sage: p = (x+y+1)^2; p
        x^2 + 2*x*y + y^2 + 2*x + 2*y + 1
        sage: list(items_of_vector(p))
        [((1, 0, 0), 2),
         ((1, 1, 0), 2),
         ((0, 0, 0), 1),
         ((2, 0, 0), 1),
         ((0, 1, 0), 2),
         ((0, 2, 0), 1)]

    univariate polynomials::

        sage: P = ZZ['x']
        sage: x = P.gen()
        sage: (x+2)^3
        x^3 + 6*x^2 + 12*x + 8
        sage: list(items_of_vector(_))
        [(0, 8), (1, 12), (2, 6), (3, 1)]

    elements of quotients::

        sage: C = CyclotomicField(5)
        sage: z = C.gen()
        sage: p = (z+2)^2; p
        zeta5^2 + 4*zeta5 + 4
        sage: list(items_of_vector(p))
        [(0, 4), (1, 4), (2, 1)]
    """
    if isinstance(v, CombinatorialFreeModule.Element):
        return v
    else:
        try:
            return v.dict().items()
        except AttributeError:
            return items_of_vector(v.lift())

def leading_item(v):
    return v.leading_item()

def coeff(v, i):
    return v[i]

def is_zero(v):
    """
    Return whether `v` is zero
    """
    return not v