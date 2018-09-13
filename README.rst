Chrom_tools
-----------

#Usage

To plot chromatograms, simply do::

    >>> import chrom_tools
    >>> c = chrom_tools.ChromatogramsContainer()
    >>> ok, msg = c.load_chromatogram(folder + 'LAG525 AEC VCS MuLV Run1 001.csv')
    >>> chrom = c[-1]
    >>> HETPs = chrom.HETP_calc('Cond', [0,20], 19.8, plot = True)
	
To plot analytical data:

- plot eluates and load of the selected runs:

    >>> import chrom_tools
    >>> c = chrom_tools.ChromatogramsContainer()
    >>> ok, msg = c.load_chromatogram(folder + 'LAG525 AEC VCS MuLV Run1 001.csv')
    >>> chrom = c[-1]
    >>> HETPs = chrom.HETP_calc('Cond', [0,20], 19.8, plot = True)