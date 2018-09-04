Chrom_tools
-----------

To use, simply do::

    >>> import chrom_tools
    >>> c = chrom_tools.ChromatogramsContainer()
    >>> ok, msg = c.load_chromatogram(folder + 'LAG525 AEC VCS MuLV Run1 001.csv')
    >>> chrom = c[-1]
    >>> HETPs = chrom.HETP_calc('Cond', [0,20], 19.8, plot = True)