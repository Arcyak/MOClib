"""
Instructions:


1. Extract table in VOTable format (see the examples in section 4, p.3)).
Choose an API todo that (eg: astropy, astroquery or pyvo).
Note: choose a table having position and query a reasonable number of records.


2. Build an astropy table with an added HEALPix column in the order that you choose (eg: order 10) You
can use the astropy-healpix library to compute the HEALPix number from ra,dec.


3. Create a simple HEALPix pixelation from the HEALPix column (so a kind of MOC but in a unique
order). Your code must contains a text serialization.
Example of HEALPix pixelation serialization: 6/17407 18090 18773 19456


4. Extend the code to build a real MOC (with HEALPix in different order) and apply it to the HEALPix
column. Your code must contains a text serialization (like example 2.2).
Todo that, you can (but you can follow an other algorithm) have a recursive approach (see 2.2) that
manages cells in a structure (made of dictionaries, lists, objects) and functions/methods that manages
node operations (like to add or delete a node). Then you iterate the HEALPix values and feed one by one
the MOC tree.

5. Illustrate the result in a Jupyter notebook with plots (matplotlib) or in Aladin (module ipyaladin).
In the notebook, use your serialization and mocpy (especially for functions like plots or moc-operation
only available in mocpy).


Note: to test, you can use the mocpy library and check if the results are similar. You can also plot the result
using astropy (see the mocpy doc.). The mocpy library is able to check some inconsistencies like overlapping
elements.

Note: for this project a debugger could help (for instance pdb)
"""

from astropy.table import Table


def get_votable(catalogue: str, out_max: int = 10000):
    """get VOTable from a VizieR catalogue for further studies.
    :param catalogue : VizieR catalogue name
    :param out_max: max number of objects, default: 10000
    """

    url = f'https://vizier.cds.unistra.fr/viz-bin/votable?-source={catalogue}&-out.max={out_max}'
    return Table.read(url)


# self.nside = 2**self.order


if __name__ == '__main__':
    print(get_votable('II/7A/catalog', 100))

