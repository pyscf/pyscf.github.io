.. _developer_ao2mo:

***********************
AO to MO transformation
***********************

*Modules*: :mod:`ao2mo`

Overview
========
The :mod:`ao2mo` module exposes users with :func:`kernel`, :class:`load`, :func:`restore` to handle 4-index AO to MO transformation.
The function :func:`kernel` takes in :class:`Mole` object or the ao integrals stored in a :class:`numpy.ndarray` and a set of mo coefficients
If ``erifile`` is specified, the integrals will be stored in a HDF5 file, otherwise it will be stored in a numpy array.
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
    # eri is a numpy array
    eri = ao2mo.kernel(mol, orb)

    ftmp = tempfile.NamedTemporaryFile()
    # saves the two-electron integrals in the file ftmp.name
    ao2mo.kernel(mol, orb, ftmp.name)
    

It can also be called as a method :meth:`ao2mo` of the :class:`Mole` object
::
    mol.ao2mo(orb) 

The class :class:`load` helps users load two-electron integrals from a hdf5 file.
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

    :arg eri_or_mol: It can either be a four-dimensional array that stores the AO integrals explicitly or a :class:`Mole` object. If it is an array, everything will be kept in memory, the incore algorithm will be used. If it is an :class:`Mole` object, the AO integrals will be computed on the fly, the outcore algorithm will be used.
    :arg mo_coeffs: It can be either a single set of mo coefficients in numpy array or a list of four sets of mo coefficients. Each of the four sets of mo coefficients correspond to a index in (ij|kl). If only one is provided, the four indices will correspond to the same mo coefficients.
    :keyword erifile: It is the name of the hdf5 file in which the integrals are stored.If the ``eri_or_mol`` argument is an numpy array, :func:`kernel` will call the incore algorithm to perform the transformation,this argument will then be of no use. If specified, the integrals will be stored in the HDF5 file or the related group.If not specified, pyscf will use an anonymous temp file and returns a ``numpy.ndarray`` in the end.
    :type erifile: str or :class:`h5py.Group` object or :class:`h5py.File` object
    :keyword str dataname: ``dataname`` labels the integrals stored in the erifile. The integrals can be reused by assigning different dataname. If the erifile already contains the given dataname, the old integrals will be overwritten. 
    :keyword str intor: the name of the integral you want to evaluate. More details can be found at :mod:`gto`.
    :keyword int comp: the component of the integral to be evaluated. It is closely related to ``intor``, more details can also be found at :mod:`gto`.

    :keyword aosym: It specifies to what level can the :meth:`ao2mo` utilizes symmetry in ao integrals. Supported symmetry labels are the same as those in :func:`restore`. For more details see :ref:`Transform integrals between symmetries` section below. Default aosym is 's4'. 
    :type aosym: int or str

    :keyword bool compact: When it is ``True``, the returned MO integrals has (up to 4-fold) permutation symmetry. When it is ``False``, the function will abandon any permutation symmetry, and return the "plain" MO integrals without any permutation symmetry.

Load the integrals
==================
Since integrals are stored in a HDF5 file when outcore algorithm is used.
Pyscf provides the :class:`load` to help users access the integrals stored in the HDF5 file.
It takes a flexible ranges of objects including a ``str`` which is the name of the HDF5 file, a :class:`h5py.File` object, a :class:`h5py.Group` object and a numpy array.
The dataname can also be taken as a second argument.
In this way, it helps users to treat integrals stored in memory and in file on the same footing.
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
It takes a symmetry label, a ``numpy.ndarray`` as the 4-index quantity and the ``norb`` of this quantity.
The symmetry label specifies the symmetry the user desires.
The symmetry of the input array is determined by its shape.
The four indices of this quantity must have the same dimension.
The relation between different symmetry and shape is explianed below.

The "plain" MO integrals or 's1' symmetry have the shape of (norb, norb, norb, norb).
If a pair of indices have permutation symmetry, then this pair of indices reduces to one index, the dimension becomes npair = norb*(norb+1)/2.
So integrals having 's2ij' and 's2kl' symmetries have the shape of (npair, norb, norb) or (norb, norb, npair) respectively.
If a integral has 's4' symmetry, then there are permutation symmetries between both ij indices and kl indices and it will have the shape of (npair, npair).
If a integral has 's8' symmetry, then the permutation symmetry also exists between ij pair and kl pair.
It will be a one-dimension array with length npair*(npair+1)/2.
The function determines the symmetry of the input array by its shape based on the above rule, and the shape of the output is also determined in this way.
If the input doesn't have the shape of any of the symmetries, pyscf will throw an error.

Listed are the symmetry label that this function takes. It can be either ``str`` or ``int``.

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