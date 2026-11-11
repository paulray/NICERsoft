"""Reading NICER filter (.mkf) files."""

import numpy as np
import astropy.io.fits as pyfits
from astropy.table import Table


def read_mkf(mkfname, hdu=1):
    """Read an .mkf file into an astropy Table.

    Some .mkf files declare TNULL values that cannot be represented in the
    column's unsigned integer type (e.g. TNULL = -1 on a uint8 column).
    With numpy >= 2, astropy's Table.read raises a TypeError on these when
    it sets the masked column's fill_value. Since such a null value can never
    occur in the data, it is dropped before converting to a Table.
    """
    with pyfits.open(mkfname) as hdul:
        tabhdu = hdul[hdu]
        for col in tabhdu.columns:
            if col.null is None:
                continue
            # Use the scaled dtype (TZERO applied), which is what Table.read uses
            dt = tabhdu.data[col.name].dtype
            if dt.kind in "iu":
                info = np.iinfo(dt)
                if not (info.min <= col.null <= info.max):
                    col.null = None
        return Table.read(tabhdu)
