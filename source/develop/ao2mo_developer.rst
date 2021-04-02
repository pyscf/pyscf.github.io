.. _developer_ao2mo:

***********************
AO to MO transformation
***********************

*Modules*: :mod:`ao2mo`

Overview
========
The :mod:`ao2mo` module exposes important functionality for transforming 2-electron integrals, i.e. the 4-index transformation from atomic orbitals (AOs) to molecular orbitals (MOs).
The most important functions are :func:`kernel`, :class:`load`, :func:`restore`. 
The function :func:`kernel` takes in a :class:`Mole` object or the AO integrals stored in a :class:`numpy.ndarray` and a set of MO coefficients.
If ``erifile`` is specified, the integrals will be stored in a HDF5 file with the given name, otherwise the integrals are returned as a numpy array.
::
    from pyscf import gto, scf, ao2mo
    import tempfile
    mol = gto.Mole()
    mol.build(
        atom = 'H 0 0 0; F 0 0 1.1',  # in Angstrom
        basis = 'ccpvdz',
        symmetry = True,
    )

    mf = scf.RHF(mol)
    mf.kernel()

    orb = mf.mo_coeff
    # get the two-electron integrals as a numpy array
    eri = ao2mo.kernel(mol, orb)

    ftmp = tempfile.NamedTemporaryFile()
    # saves the two-electron integrals in the file ftmp.name
    ao2mo.kernel(mol, orb, ftmp.name)
    

The :func:`kernel` function can also be invoked as the :meth:`ao2mo` method of the :class:`Mole` object
::
    mol.ao2mo(orb) 

The class :class:`load` allows users to load two-electron integrals from a hdf5 file.
::
    ftmp = tempfile.NamedTmpFile()
    ao2mo.kernel(mol, orb, ftmp)
    with ao2mo.load(ftmp) as eri:
        print(eri.shape)

The function :func:`restore` restores a set of integrals to the desired permutation symmetry.
::
    import numpy
    with ao2mo.load(ftmp) as eri:
        eri1 = ao2mo.restore(1, numpy.asarray(eri), orb.shape[1])
        eri4 = ao2mo.restore('4', numpy.asarray(eri), orb.shape[1)
        eri8 = ao2mo.restore('s8', numpy.asarray(eri), orb.shape[1])
        print(eri1.shape)
        print(eri4.shape)
        print(eri8.shape)

The kernel function
===================
.. py:function:: kernel(eri_or_mol, mo_coeffs, erifile=None, dataname='eri_mo', intor='int2e', *args, **kwargs)

    :arg eri_or_mol: This is either a four-dimensional array that stores the AO integrals explicitly or a :class:`Mole` object. If it is an array, everything will be kept in memory, and an incore algorithm will be used. If it is a :class:`Mole` object, AO integrals will be computed on the fly, and the outcore algorithm will be used.
    :arg mo_coeffs: This can be either a single set of MO coefficients in a numpy array, or a list of four sets of MO coefficients. Each of the four sets of MO coefficients correspond to an index in (ij|kl). If only one is provided, the four indices will correspond to the same MO coefficients.
    :keyword erifile: This is the name of the HDF5 file to store the integrals in. If the ``eri_or_mol`` argument is a numpy array, :func:`kernel` will call the incore algorithm to perform the transformation, and this argument will not be used. If specified, the integrals will be stored in the HDF5 file or the related group. If not specified, PySCF will use an anonymous temp file and return a ``numpy.ndarray`` in the end.
    :type erifile: str or :class:`h5py.Group` object or :class:`h5py.File` object
    :keyword str dataname: ``dataname`` labels the integrals stored in the erifile. The integrals can be reused by assigning different datanames. If the erifile already contains the given dataname, the old integrals will be overwritten. 
    :keyword str intor: the name of the integral you want to evaluate. More details can be found in :mod:`gto`.
    :keyword int comp: the component of the integral to be evaluated. It is closely related to ``intor``, more details can also be found in :mod:`gto`.

    :keyword aosym: It specifies to what level will the :meth:`ao2mo` utilize symmetry in ao integrals. Supported symmetry labels are the same as those in :func:`restore`. For more details see :ref:`Transform integrals between symmetries` section below. Default aosym is 's4'.
    :type aosym: int or str

    :keyword bool compact: When this is ``True``, the returned MO integrals have (up to 4-fold) permutation symmetry. When this is ``False``, the function will abandon any permutation symmetry, and return the "plain" MO integrals without any permutation symmetry.

Load the integrals
==================
Since integrals are stored in a HDF5 file when the outcore algorithm is used,
Pyscf provides the :class:`load` to help users access the integrals stored in the HDF5 file.
It takes a flexible ranges of objects including a ``str`` which is the name of the HDF5 file, a :class:`h5py.File` object, a :class:`h5py.Group` object and a numpy array.
The dataname can also be taken as a second argument.
In this way, it helps users access integrals stored in memory and in a file on the same footing.
The :class:`load` class should only be used within a ``with`` statement.
::
    eri = ao2mo.kernel(mol, orb)
    print(type(eri))
    with ao2mo.load(erifile) as eri:
        print(eri.shape)
    
    ao2mo.kernel(mol, orb, erifile = 'hf.h5', dataname = 'test')
    # load 2e integrals by filename and dataname
    with ao2mo.load('hf.h5', 'test') as eri:
        print(eri.shape)
    
    import numpy
    erirand = numpy.random.random((5,5,5,5))
    # load 2e integrals from numpy array
    with ao2mo.load(erirand) as eri:
        print(eri.shape)

Transform integrals between symmetries
======================================
The function :func:`restore` can transform 2e integrals or any 4-index quantity (e.g. 2rdm) between different permutation symmetries.
It takes a symmetry label, a ``numpy.ndarray`` as the 4-index quantity and the dimension ``norb`` of this quantity (the four indices
must have the same dimension).
The symmetry label specifies the output symmetry the user desires, while
the symmetry of the input array is determined by its shape.
The relation between different symmetries and shape is described below.

"Plain" MO integrals or 's1' symmetry have a shape of (norb, norb, norb, norb).
If a pair of indices have permutation symmetry, then only the triangular index is used and the shape is (npair, npair) with npair = norb*(norb+1)/2.
Thus integrals with 's2ij' and 's2kl' symmetries have a shape of (npair, norb, norb) or (norb, norb, npair) respectively.
If the integrals have 's4' symmetry, then there is permutation symmetry between both ij indices and kl indices and the shape will be (npair, npair).
If the integral have 's8' symmetry, then permutation symmetry also exists between the ij pair and kl pair.
This will yield a one-dimensional array with length npair*(npair+1)/2.
:func:`restore` determines the symmetry of the input array based on the above rules, and the shape of the output is also determined in this way.
If the input does not correspond to any of the above shapes, PySCF will throw an error.

Listed are the symmetry labels that can be used, whichcan be either be a ``str`` or ``int``.

============== ====
's8', '8', 8   8-fold symmetry
's4, '4', 4    4-fold symmetry
's2kl', '2kl'  2-fold symmetry between ij indices.
's2ij', '2ij'  2-fold symmetry between kl indices.
's1', '1', 1   1-fold symmetry or no symmetry.
============== ====

Note
====
The examples in this document can be found as a single python script in :source:`examples/ao2mo/02-ao2mo_doc.py`,
more examples can also be found at :source:`examples/ao2mo/02-ao2mo_doc.py` directory.
