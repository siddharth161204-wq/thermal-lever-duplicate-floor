"""Reading MATLAB v7.3 (HDF5) structs of the fast-charging data with h5py, for Stage 6 (DECISIONS.md D83, D84).

The structs are those of data.matr.io/1, one file per batch, each holding a 1 x n struct array `batch` whose fields are
stored as n x 1 datasets of object references. A MATLAB char array is a uint16 dataset with MATLAB_class 'char'; a
cell array is a dataset of references with MATLAB_class 'cell'; an empty array carries MATLAB_empty = 1. A MATLAB
column vector of k entries is stored with HDF5 shape (1, k), a row vector with (k, 1).

Runs on Python 3.10 with numpy and h5py only (the laptop's shell, where the files are).
"""
import numpy as np
import h5py


def mclass(obj):
    c = obj.attrs.get("MATLAB_class", b"")
    return c.decode() if isinstance(c, bytes) else str(c)


def deref(f, x):
    return f[x] if isinstance(x, h5py.Reference) else x


_HEAP = {}


def string_heap(f):
    """The MATLAB string objects of a file, decoded from '#subsystem#/MCOS' (MATLAB_class FileWrapper__).

    A string object is stored in the struct as a uint32 reference [0xDD000000, 2, 1, 1, object id, class id]; its
    characters are in the MCOS cell array, whose first cell is the metadata and second an empty placeholder, the
    object of id k holding its value in cell k + 1 (counting from 0) as a uint64 array: a version (1), the number of
    dimensions, the dimensions, the number of characters of each string, then the UTF-16 characters four to a word.
    The mapping is checked against the site's notes (DECISIONS.md D84): the 2017-06-30 cells it names by channel and
    barcode decode to that channel and barcode."""
    key = f.filename
    if key in _HEAP:
        return _HEAP[key]
    heap = {}
    if "#subsystem#" in f and "MCOS" in f["#subsystem#"]:
        refs = np.asarray(f["#subsystem#"]["MCOS"][()]).ravel()
        for j in range(2, len(refs)):
            o = f[refs[j]]
            if not isinstance(o, h5py.Dataset) or o.dtype != np.uint64:
                continue
            w = np.asarray(o[()]).ravel()
            if len(w) < 4 or int(w[0]) != 1:
                continue
            nd = int(w[1])
            ns = int(np.prod(w[2:2 + nd]))
            lens = [int(v) for v in w[2 + nd:2 + nd + ns]]
            chars = w[2 + nd + ns:].astype("<u8").view("<u2")
            out, pos = [], 0
            for n in lens:
                out.append("".join(chr(int(c)) for c in chars[pos:pos + n]))
                pos += n
            heap[j - 1] = out[0] if ns == 1 else "|".join(out)
    _HEAP[key] = heap
    return heap


def text(f, x):
    """A char array, a MATLAB string object, or a cell holding either, as (str, class); (None, class) otherwise."""
    obj = deref(f, x)
    if isinstance(obj, h5py.Group):
        return None, "group"
    cls = mclass(obj)
    if int(obj.attrs.get("MATLAB_empty", 0)):
        return "", cls + ":empty"
    if cls == "char":
        a = np.asarray(obj[()]).ravel(order="C")
        return "".join(chr(int(c)) for c in a), cls
    if cls == "string" and int(obj.attrs.get("MATLAB_object_decode", 0)) == 3:
        w = np.asarray(obj[()]).ravel()
        if len(w) == 6 and int(w[0]) == 0xDD000000 and int(w[1]) == 2 and int(w[2]) == 1 and int(w[3]) == 1:
            v = string_heap(f).get(int(w[4]))
            return v, "string" if v is not None else "string:undecoded"
        return None, "string:array"
    if cls == "cell":
        parts = [text(f, r) for r in np.asarray(obj[()]).ravel(order="F")]
        if all(p[0] is not None for p in parts):
            return "|".join(p[0] for p in parts), "cell:" + ",".join(sorted(set(p[1] for p in parts)))
        return None, "cell:" + ",".join(sorted(set(p[1] for p in parts)))
    return None, cls


def orientation(ds):
    """'column' for HDF5 shape (1, k), 'row' for (k, 1); the length itself is not returned."""
    sh = ds.shape
    if len(sh) != 2 or (sh[0] != 1 and sh[1] != 1):
        raise ValueError(f"not a MATLAB vector: {ds.name}")
    return "column" if sh[0] == 1 else "row"


def prefix(ds, m):
    """The first m entries of a MATLAB vector (fewer if it is shorter), read without reading the rest."""
    if orientation(ds) == "column":
        return np.asarray(ds[0, :m], dtype=float)
    return np.asarray(ds[:m, 0], dtype=float)
