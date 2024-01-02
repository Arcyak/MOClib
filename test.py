"""
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
import time
from typing import List, Dict

from astropy import log
from astropy.table import Table, Column
from astropy_healpix import HEALPix
from astropy.coordinates import Galactic, SkyCoord
import astropy.units as u

import re

log.setLevel('ERROR')


class Ratio:
    def __init__(self, __order: int = 10) -> None:
        """
        :param __order: Order of HEALPix calculation, default:10
        """
        self.n_side = 2**__order
        self.hp_object = HEALPix(nside=self.n_side, order='nested', frame=Galactic())
        self.moc_tree = Moc_tree(order=order)

    def build_healpix_table(self, __table: Table) -> Table:
        """Build an astropy table with an added HEALPix column in the order that you choose (eg: order 10) You
        can use the astropy-healpix library to compute the HEALPix number from ra,dec.

        :param __table: Table
        """

        # Create a new column name HEALPix from the RA and DE columns
        __RA: str = next((col for col in __table.columns if re.match(r'.*RA.*', col)), None)
        __DE: str = next((col for col in __table.columns if re.match(r'.*DE.*', col)), None)

        new_column = (
            Column([self.hp_object.skycoord_to_healpix(SkyCoord(__table[__RA][i], __table[__DE][i], unit="deg"))
                    for i in range(len(__table))], name="HEALPix", dtype=int))

        # Add the column to the table
        __table.add_column(new_column)

        # Build the MOC tree
        self.moc_tree.build_moc_tree(__table)

        return __table


class Moc_tree:
    """TODO MOC TREE"""
    def __init__(self, order: int) -> None:
        self.orders: Dict[int, List[int]] = {order: []}

    def build_moc_tree(self, table: Table) -> None:
        """Build a MOC tree from a VOTable having a HEALPix column"""

        # Check if the HEALPix column exist
        assert 'HEALPix' in table.columns, \
            'No HEALPix column in the given table. Please build it before using this method'

        # Add HEALPix indices to the MOC tree
        order = self.orders.copy()
        order_keys = sorted(order.keys())

        for i, healpix in enumerate(table['HEALPix']):
            order[order_keys[0]].append(healpix)

        # Add HEALPix indices to the MOC tree in a hierarchical way
        for current_order in range(order_keys[0], order_keys[-1]):
            current_level = order[current_order]
            next_order = current_order + 1
            next_level = []
            for ipix in current_level:
                # Split each HEALPix cell into its sub-cells of the next order
                next_level.extend(list(HEALPix(nside=2**(next_order), order='nested').subpix2pix(ipix)))

            order[next_order] = next_level

        self.orders = order

    def serialize_moc(self) -> str:
        """Serialize the MOC tree"""

        moc_serialization = ""
        for order, ipix_list in self.orders.items():
            moc_serialization += f"{order}/{'-'.join(map(str, ipix_list))} "

        return moc_serialization.rstrip()

    def add_node(self, order: int, ipix: int):
        """Add a node to the tree"""
        if order not in self.orders:
            self.orders[order] = [ipix]
        else:
            self.orders[order].append(ipix)

    def del_node(self, order: int, ipix: int):
        """Del a node from the tree"""
        if order in self.orders and ipix in self.orders[order]:
            self.orders[order].remove(ipix)

    pass


def get_votable(__catalogue: str, __out_max: int = 10000) -> Table:
    """Get VOTable from a VizieR catalogue for further studies.
    TODO: CONVERT COORDINATES TO J2000 (if necessary) !

    :param __catalogue: Name of the catalogue
    :param __out_max: Max number of entries in the table, default:10000
    """

    __url = f'https://vizier.cds.unistra.fr/viz-bin/votable?-source={__catalogue}&-out.max={__out_max}'

    # Query the VOTable from VizieR
    __table = Table.read(__url, format='votable')
    assert __table is not None and len(__table) > 0, f"Table {__catalogue} is empty or does not exist."

    # Check the coordinate system and if needed converts the RA and DE to J2000
    __table = update_coordinate_system(__table)

    return __table


def serialize_healpix(__table: Table, __order: int) -> str:
    """Serialize the HEALPix from a HEALPix column for a given order."""

    __serialization = f'{__order}/'
    for i, healpix in enumerate(__table['HEALPix']):
        __serialization += f'{healpix} '

    return __serialization


def update_coordinate_system(__table: Table):
    # Check if RA and DEC columns exist
    __RA: str = next((col for col in __table.columns if re.match(r'.*RA1950.*|_RA$', col)), None)
    __DE: str = next((col for col in __table.columns if re.match(r'.*DE1950.*|_DE$', col)), None)

    if __RA and __DE:
        # Convert B1950 coordinates to J2000 if __RA and __DEC columns exist
        b1950_coords = SkyCoord(ra=__table[__RA], dec=__table[__DE], unit=(u.deg, u.deg), frame='fk4')
        j2000_coords = b1950_coords.transform_to('fk5')

        # Update the table with J2000 coordinates
        __table['RAJ2000'] = j2000_coords.ra.deg
        __table['DEJ2000'] = j2000_coords.dec.deg

        # Remove B1950 columns if needed
        __table.remove_columns([__RA, __DE])

        print("Coordinates converted to J2000.\n")
    return __table


# Function to serialize the MOC tree
def serialize_moc_tree(moc_tree: Moc_tree) -> str:
    return moc_tree.serialize_moc()

if __name__ == '__main__':
    # 'II/7A/catalog' : UBVRIJKLMNH Photoelectric Catalogue (Morel+ 1978)
    # 'J/ApJ/831/67/table1' : Galactic CHaMP. III. 12CO dense clump properties (Barnes+, 2016)
    # 'J/ApJ/804/L15' : SDSS-DR7 broad-line QSOs (Sun+, 2015)
    # 'VII/26D/catalog' : Uppsala General Catalogue of Galaxies (UGC) (1973)
    # 'I/298/table3' :  Catalog of Northern stars with annual proper motions larger than 0.15" (2005)
    catalogue = 'II/7A/catalog'
    out_max = 10000
    order = 10

    hp = Ratio(order)
    t = get_votable(__catalogue=catalogue, __out_max=out_max)
    t1 = hp.build_healpix_table(t)

    healpix_table = hp.build_healpix_table(t1)

    # Serialize the MOC tree
    moc_tree_serialization = serialize_moc_tree(hp.moc_tree)

    # Display the serialized MOC tree
    print("Serialized MOC Tree:", moc_tree_serialization)
    print('-------------')
    print(serialize_healpix(t, order))
