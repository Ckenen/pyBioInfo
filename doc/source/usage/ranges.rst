.. _usage-ranges:

Ranges and Genomic Ranges (the Range class hierarchy)
======================================================

pyBioInfo provides a small family of "Range" classes, headed by
``GRange`` — a genomic range which also carries all the plain
one-dimensional interval features. They all share the same coordinate
conventions:

* **Half-open interval** ``[start, end)``: ``start`` is the inclusive,
  0-based coordinate and ``end`` is the exclusive coordinate, hence
  ``len(range) == end - start``.
* **Span vs. effective length**: in the block-structured classes, the
  ``start``/``end`` attributes denote the outer span from the first
  block's start to the last block's end, while ``len()`` counts only
  the within-block length (gaps between blocks are excluded).

Class hierarchy
---------------

.. code-block:: text

   IRange          minimal 1-D half-open interval: read-only
                   start / end coordinates only
     └── CRange    chromosome range: + the chrom attribute
           └── GRange  genomic range: + name / strand and multi-block
                       (blocks) structure
               ├── TRange  thick range: + thick (e.g. BED12 thick columns)
               └── MRange  multiple block groups: + blocks_array
                           (e.g. paired-end fragments)

   BlockTools      static toolkit for lists of blocks (the foundation
                   of the classes above; see :ref:`usage-blocktools`)

.. _usage-irange:

IRange: minimal interval
------------------------

``IRange(start, end)`` is the plain-coordinate base class of the
Range family: a half-open interval on the number line exposing the
two read-only attributes ``start`` and ``end`` plus the
index/position conversions between absolute coordinates and 0-based
offsets within the interval.

Instances also implement a fixed 2-element sequence protocol:
``len()`` returns the interval length, iteration and ``[0]``/``[1]``
indexing yield the start and end coordinates (item assignment is
forbidden), and ``str()`` renders as ``"<ClassName>: <start>-<end>"``.

.. code-block:: python

   >>> from pyBioInfo.Range import IRange

   >>> obj = IRange(10, 20)
   >>> obj.start, obj.end       # read-only coordinates
   (10, 20)

**Construction constraints**: ``start >= 0`` and ``start < end``,
otherwise a ``ValueError`` is raised.

**Index/position conversion**:

.. code-block:: python

   >>> obj.get_index(12)        # absolute coordinate -> 0-based index
   2
   >>> obj.get_position(-1)     # index -> absolute coordinate, negative ok
   19

**Sequence protocol**:

.. code-block:: python

   >>> len(obj)                 # interval length
   10
   >>> start, end = obj         # iterable, equivalent to obj[0], obj[1]
   (10, 20)
   >>> str(obj)
   'IRange: 10-20'

**Ordering**: instances are ordered first by ``start``, then by
``end`` — a total order, so they can be used directly with
``sorted`` / ``min`` and friends. The ordering is implemented by the
``_compare`` method, which subclasses extend with higher-level keys
(e.g. the chromosome name in ``CRange``).

**Distance relations** (each compares against another ``IRange``;
all are decided by the ``get_distance`` method on the coordinates
alone — negative means overlap, zero adjacency, positive a gap):

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Method
     - Meaning
   * - ``is_overlap(other)``
     - the two intervals share at least one position
   * - ``is_adjacent(other)``
     - the two intervals touch end-to-start (no overlap, no gap)
   * - ``is_away(other)``
     - the two intervals are separated by a non-zero gap
   * - ``is_contact(other)``
     - overlapping or adjacent (no gap in between)

.. code-block:: python

   >>> obj.is_overlap(IRange(15, 25))
   True
   >>> obj.is_away(IRange(25, 30))
   True

.. _usage-crange:

CRange: chromosome range
------------------------

``CRange(chrom, start, end)`` extends ``IRange`` with the chromosome
name ``chrom`` (a plain read/write attribute) and the optional strand
``strand``.

.. code-block:: python

   >>> from pyBioInfo.Range import CRange

   >>> obj = CRange("chr1", 10, 20)
   >>> obj.chrom, obj.start, obj.end
   ('chr1', 10, 20)

**Strand**: ``CRange(chrom, start, end, strand=...)`` accepts one of
the four values ``None``, ``"."``, ``"+"`` and ``"-"`` (validated on
assignment; invalid values raise a ``ValueError``), and offers the
convenience properties ``is_forward`` / ``is_reverse``.

.. code-block:: python

   >>> obj = CRange("chr1", 10, 20, strand="+")
   >>> obj.strand, obj.is_forward, obj.is_reverse
   ('+', True, False)

**Ordering**: from CRange onwards the ordering becomes
chromosome-aware — instances are grouped by ``chrom`` in
lexicographic order first, then ordered by the inherited ``start``
and ``end`` coordinates within the same chromosome (plain ``IRange``
operands fall back to the coordinate-only comparison).

**Strand-aware index/position conversion**: the inherited
``get_index(position, strandness=True)`` / ``get_position(index,
strandness=True)`` become strand-aware here — on the reverse strand
the internal indexing runs from the end towards the start, like a
sequence read on the minus strand; negative indexes count from the
strand-aware end. Pass ``strandness=False`` to keep the plain
forward counting.

.. code-block:: python

   >>> obj = CRange("chr1", 10, 20, strand="-")
   >>> obj.get_index(19)        # end - 1 is the first base on "-"
   0
   >>> obj.get_position(0)
   19
   >>> obj.get_index(19, strandness=False)
   9

**Distance relations**: the distance becomes chromosome-aware too —
ranges on different chromosomes can never overlap or touch, so they
are always reported as ``is_away`` (never ``is_overlap`` /
``is_adjacent`` / ``is_contact``); only ranges on the same chromosome
(compared against plain ``IRange`` operands: the coordinate-only
fallback) may overlap or be adjacent.

**Text output**: ``format(fmt, ...)`` dispatches to the three
writers ``to_bed_string(ncol)``, ``to_gff_string(attris)`` and
``to_tbl_string(header)``, all implemented here. The optional columns
they need (``name``, ``score``, ``thick``, ``color``, ``blocks`` for
BED; ``source``, ``feature``, ``score``, ``frame``, ``name`` for GFF;
``"blocks"`` and ``"thick"`` keys for TBL) are resolved via attribute
lookup — attributes not provided by the class are treated as ``None``
(rendered as ``"."`` or an empty cell where the format requires it),
and a range without a block structure is rendered as a single block
covering the whole span.

.. code-block:: python

   >>> obj.to_bed_string(ncol=4)
   'chr1\t10\t20\tNone'

.. _usage-grange:

GRange: genomic range
---------------------

``GRange`` extends ``CRange`` with feature naming and a **blocks**
structure: a genomic range lives on a chromosome (``chrom``), carries
an optional feature name (``name``) and strand (``strand``, both
inherited from ``CRange``), and its interval may be split
into several non-overlapping blocks (e.g. the exons of a transcript);
the overall ``start``/``end`` span coordinates are inherited from
``IRange``.

.. code-block:: python

   >>> from pyBioInfo.Range import GRange

   >>> obj = GRange("chr1", blocks=[(10, 15), (20, 25)])
   >>> obj.start, obj.end       # outer span: first block start to last block end
   (10, 25)
   >>> len(obj)                 # effective length: within-block positions only
   10
   >>> obj.chrom, obj.name, obj.strand
   ('chr1', None, None)
   >>> start, end = obj         # iterable as a 2-element (start, end) sequence
   (10, 25)

**Construction** (either way):

* pass ``blocks`` (``start``/``end`` may be omitted; if given they
  must agree with the first/last block);
* pass only ``start``/``end``, which degenerates to a single-block
  range.

The block list must satisfy: non-empty, ``start >= 0`` for every
block, strictly increasing order between blocks (zero-length gaps
are allowed), and ``start < end`` for every block (zero-length
blocks are not), otherwise a ``ValueError`` is raised (validation is
performed by :ref:`BlockTools <usage-blocktools>`).

**Chromosome attributes**: ``chrom`` and ``name`` are plain read/
write attributes. ``strand`` is inherited from ``CRange`` (validated
on assignment with the ``is_forward`` / ``is_reverse`` conveniences).

**Annotation attributes**: ``score``, ``source``, ``feature``,
``frame`` and ``color`` are all optional (default ``None``, rendered
as ``"."`` during formatting); they are mainly populated by
subclasses or the IO modules (BED/GFF readers).

**Ordering**: the chromosome-aware ordering is inherited from
``CRange``: instances are grouped by ``chrom`` in lexicographic
order, then sorted by ``start`` and ``end`` within a chromosome —
i.e. natural genomic coordinate order. This is a total order, so
ranges can be used directly with ``sorted`` / ``min`` and friends.
The 2-element sequence protocol (iteration, ``obj[0]`` / ``obj[1]``)
is inherited from ``IRange``; only ``len()`` is overridden to return
the total block length.

**Distance relations**: the span-based, chromosome-aware relations
(``is_overlap`` / ``is_adjacent`` / ``is_away`` / ``is_contact``)
are inherited from the base classes (IRange / CRange). GRange adds
the block-level containment tests (each compares against another
``GRange``):

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Method
     - Meaning
   * - ``is_contain(other)``
     - this range's blocks fully cover ``other``'s blocks
   * - ``is_coincide(other)``
     - ``other``'s block boundaries align with this range's blocks
       (directional; see below)

**Block-aware indexing**: the indices of ``get_index`` /
``get_position`` run over the *concatenated blocks*, automatically
skipping inter-block gaps, and honor ``strandness`` for the reverse
strand:

.. code-block:: python

   >>> obj = GRange("chr1", blocks=[(10, 15), (20, 25)])
   >>> obj.get_index(21, strandness=False)   # skips the 15-20 gap
   6
   >>> obj.get_position(6, strandness=False)
   21

**Block-level relations**: when both operands are ``GRange``
instances, ``is_contain`` / ``is_coincide`` are decided on block
boundaries (not the outer span):

* ``is_contain(other)``: the blocks of this range **fully cover**
  every block (or block fragment) of ``other``;
* ``is_coincide(other)``: the block boundaries of ``other`` are
  aligned with blocks of this range — ``other`` may start later in
  its first block and end earlier in its last block (i.e. ``other``
  is a "trimmed" subset of this range's block layout). The test is
  **directional**: ``a.is_coincide(b)`` being true does not imply
  ``b.is_coincide(a)``.

**Text output**: ``format()`` supports three formats — ``"BED"``
(pass ``ncol`` to select 1-12 columns; default is the 12-column
BED12), ``"GFF"`` (pass ``attris`` to append attribute columns) and
``"TBL"`` (requires ``header`` naming the column keys) — each
producing a single tab-delimited line. All three writers are
inherited from ``CRange`` (which resolves the optional columns via
attribute lookup), with the block and thick columns filled from the
``blocks`` structure and ``thick`` interval.

.. _usage-trange:

TRange: thick range
-------------------

``TRange`` extends ``GRange`` with an optional **thick** range, i.e.
the coding portion marked by columns 7-8 of BED12 (such as the CDS
interval of a transcript).

.. code-block:: python

   >>> from pyBioInfo.Range import TRange

   >>> obj = TRange("chr1", blocks=[(10, 20), (40, 50)], thick=(40, 50))
   >>> obj.thick
   (40, 50)

**Construction validation**: ``thick=(start, end)`` must fall inside
the block structure — ``thick_start >= start``, ``thick_end <= end``,
and both positions ``thick_start`` and ``thick_end - 1`` must be
covered by some block; otherwise a ``ValueError`` is raised.

**``indexes()``**: returns the indices of the thick range on the
concatenated blocks as the tuple ``(0, thick_start_index,
thick_end_index, length)`` (the two middle entries are ``None``
when no thick range is set). With ``strandness`` (default ``True``)
the counting is reversed on the minus strand; results are cached
per strand. These indices are commonly used to map genomic
coordinates onto the concatenated "transcript coordinate system".

.. code-block:: python

   >>> obj.indexes()
   (0, 10, 20, 20)

**Text output**: the writers are inherited from ``CRange``, whose
attribute lookup picks up the ``thick`` interval to fill the BED
columns 7-8 (falling back to ``start`` when absent) and the
``"thick"`` key of ``format("TBL")``.

.. _usage-mrange:

MRange: multiple block groups
-----------------------------

``MRange`` is a genomic range composed of **multiple groups of
blocks**; each group is kept separately in ``blocks_array`` (a
typical case being the alignment blocks of the two mates of a
paired-end fragment).

.. code-block:: python

   >>> from pyBioInfo.Range import MRange

   >>> obj = MRange("chr1", blocks_array=[[(10, 20), (50, 60)],
   ...                                    [(30, 40)]])
   >>> obj.blocks_array        # the two original groups
   (((10, 20), (50, 60)), ((30, 40),))
   >>> obj.blocks              # overall blocks = sorted union (gap<=0 sutured)
   ((10, 20), (30, 40), (50, 60))

**Construction** (either way):

* pass ``blocks_array`` (an iterable of block groups): the overall
  ``blocks`` are the union of all groups, produced via
  ``BlockTools.sorted`` + ``BlockTools.suture(gap=0)`` (adjacent or
  overlapping blocks are merged); each group is validated
  independently;
* pass only ``start``/``end``: a single-group, single-block range
  is constructed.

**``coincide_from(other)``**: tests whether **every group** of
``self`` coincides with ``other.blocks`` (a directional test with
the same rules as ``BlockTools.is_coincide``; see
:ref:`usage-blocktools`) — all groups must coincide for the result
to be ``True``. This class is the base of ``Fragment`` (the fam-file
fragment) in ``pyBioInfo.IO.File``.

.. _usage-blocktools:

BlockTools: block toolkit
=========================

``BlockTools`` (in ``pyBioInfo.Utils``) is the static toolkit behind
every block-related operation. A "block" is simply a pair
``(start, end)`` denoting a half-open interval, and lists of blocks
are usually required to be sorted and non-overlapping (see
``check_blocks``). The construction validation and index arithmetic
of ``GRange`` / ``TRange`` / ``MRange`` are all backed by it.

.. code-block:: python

   >>> from pyBioInfo.Utils import BlockTools

   >>> BlockTools.check_blocks([(10, 20), (30, 40)])   # silent when valid
   >>> BlockTools.get_length([(10, 20), (30, 40)])     # total block length (no gaps)
   20
   >>> BlockTools.gaps([(10, 20), (30, 40)])           # inter-block gaps
   [(20, 30)]

Method reference
----------------

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Method
     - Description
   * - ``check_blocks(blocks)``
     - validate a block list: non-empty, start>=0, strictly
       increasing, start<end; raises on violation
   * - ``is_valid_blocks(blocks)``
     - same rules, returns ``bool`` instead of raising
   * - ``check_sorted(blocks)``
     - validate coordinate sorting (by start first, then end)
   * - ``is_sorted(blocks)``
     - same as above, returns ``bool``
   * - ``sorted(blocks)``
     - return a new, sorted list of blocks
   * - ``get_length(blocks)``
     - sum of block lengths (gaps excluded)
   * - ``get_index(blocks, pos)``
     - absolute coordinate -> 0-based index on the concatenated
       blocks (``reverse`` supported)
   * - ``get_position(blocks, i)``
     - index -> absolute coordinate (negative indices and
       ``reverse`` supported)
   * - ``is_include(blocks, pos)``
     - test whether a coordinate is covered by the blocks
   * - ``move(blocks, step)``
     - shift all blocks by ``step`` (start+step must stay >= 0)
   * - ``gaps(blocks)``
     - list the gaps between adjacent blocks
   * - ``suture(blocks, gap=0)``
     - suture: merge adjacent blocks whose gap is <= ``gap``
   * - ``fusion(blocks1, blocks2)``
     - union of two block lists (sorted + gap=0 sutured)
   * - ``is_contain(b1, b2)``
     - whether ``b1`` fully covers every block of ``b2``
   * - ``is_coincide(b1, b2)``
     - whether the block boundaries of ``b2`` align with ``b1``
       (directional; see below)
   * - ``clip(blocks, start, end)``
     - clip blocks to the window ``[start, end)``, with
       extend/supply options
   * - ``extend(...)``
     - extend blocks to both sides (not implemented yet; raises
       ``NotImplementedError``)

The directional semantics of ``is_coincide``
--------------------------------------------

``is_coincide(blocks1, blocks2)`` tests whether **blocks2 coincides
with blocks1**: every block boundary of blocks2 must be aligned with
a block of blocks1, but the following are allowed:

* the first block of blocks2 may **start later** (trimmed at the 5'
  end);
* the last block of blocks2 may **end earlier** (trimmed at the 3'
  end);
* blocks of blocks1 in the middle may be **skipped** (not covered by
  blocks2).

The test is therefore directional: ``is_coincide(a, b)`` being true
does not imply ``is_coincide(b, a)``. A typical use case is checking
whether the alignment blocks of a fragment "fall inside" the block
layout of a gene or transcript.

Window clipping with ``clip``
-----------------------------

``clip(blocks, start, end)`` keeps only the blocks (or block
fragments) within the window ``[start, end)``. When the window
extends beyond the original blocks, three padding strategies are
available:

* ``extend`` / ``extend_left`` / ``extend_right``: stretch towards
  ``template`` (a template block list, e.g. a gene structure); the
  stretched regions still come from the template blocks;
* ``supply`` / ``supply_left`` / ``supply_right``: pad the parts
  beyond the template range with plain filler blocks (bounded by the
  ``limited`` range);
* ``suture=True``: apply a final gap=0 suture to the result.

.. code-block:: python

   >>> BlockTools.clip([(10, 20), (30, 40)], start=15, end=35)
   [(15, 20), (30, 35)]
