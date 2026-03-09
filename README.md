# sage-subspace

SageMath code to do linear algebra with spans of vectors in arbitrary vector spaces (not just row or column vectors).

The code isn't particularly efficient, but it has served me well over the years.

Overview:

* ``combinatorial_vector_tuple.py``: The actual code.
  The class ``VectorTuple`` is the main workhorse.
  Other classes are experimental.

* ``sga-bases.py``: Some bases of the symmetric group
  algebra. Scratch work, but might be useful.

* ``irdelta_test.py``, ``vectors.pyx``: Old scratch
  work whose purpose I don't remember.

-- Darij