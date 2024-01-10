"""
Python custom MOC library that creates MOC objects from HEALPix.

Author: Yannick HÉNIN
University of Strasbourg
"""

import re

import astropy.units as u
import numpy as np
from astropy import log
from astropy.coordinates import SkyCoord
from astropy.io import votable
from astropy.table import Table, Column
from astropy_healpix import HEALPix

log.setLevel("ERROR")


class Loader:
    """Loader class, allow one to query a Table from Vizier and do preliminary work on the Table object"""
    def __init__(self, _catalogue: str, _out_max: int = 10000, _order: int = 10, _source: str = 'CDS') -> None:
        """
        :param _catalogue: Name of the catalogue
        :param _out_max: Max number of entries in the table, default:10000
        :param _order: Order of HEALPix calculation, default:10
        """

        self.__catalogue = _catalogue

        match _source:
            case 'CDS':
                # CDS Vizier (Strasbourg)
                self.__url = f"https://vizier.cds.unistra.fr/viz-bin/votable?-source={_catalogue}&-out.max={_out_max}"

            case "CFA":
                # CFA Vizier (Harvard)
                self.__url = f"https://vizier.cfa.harvard.edu/viz-bin/votable?-source={_catalogue}&-out.max={_out_max}"

        self.__n_side = 2**_order
        self.__healpix = HEALPix(nside=self.__n_side, order="nested", frame="icrs")

    def get_votable(self) -> Table:
        """Get VOTable from a Vizier catalogue for further studies"""

        # Query the VOTable from Vizier
        print(f"Querying: '{self.__catalogue}' from Vizier")
        table = Table.read(self.__url, format="votable")
        assert table is not None and len(table) > 0, f"Table {self.__catalogue} is empty or does not exist."
        vtable = votable.parse(self.__url)

        # Check the coordinate system and if needed converts the RA and DE to J2000
        for coosys in vtable.iter_coosys():
            print(f"Coordinate system of the table: {coosys.system}\n")
            table = update_coordinate_system(table, coosys.system)

        return table

    def build_healpix_table(self, _table: Table, _rm_duplicates: bool = True) -> Table:
        """Build an astropy table with an added HEALPix column in an order.
        Use the astropy-healpix library to compute the HEALPix number from ra,dec

        :param _table: Table from which to build the HEALPix column
        :param _rm_duplicates: bool to decide whether to remove the duplicated HEALPix numbers or not
        """

        tmp = _table.copy()

        # Create a new column name HEALPix from the RA and DE columns
        ra: str = next((col for col in _table.columns if re.match(r".*RA.*", col)), None)
        de: str = next((col for col in _table.columns if re.match(r".*DE.*", col)), None)

        healpix_col = self.__healpix.lonlat_to_healpix(lon=np.array(_table[ra])*u.deg, lat=np.array(_table[de])*u.deg)
        new_column = Column(healpix_col, name="HEALPix", dtype=int)

        # Add the column to the table
        tmp.add_column(new_column)

        # Remove the duplicated Healpix lines
        if not _rm_duplicates:
            print(f"Duplicates were not removed.\n")
        else:
            tmp = remove_duplicates(tmp)

        return tmp


class Moc_tree:
    """Moc_tree class able to separate the HEALPix from a Table object into HEALPix of different orders"""
    def __init__(self, _order: int) -> None:
        self.nodes = {_order: []}

    def build_moc_tree(self, _table: Table) -> None:
        """Build a MOC tree from a VOTable having a HEALPix column

        :param _table: Table used to build the MOC tree"""

        # Check if the HEALPix column exist
        assert "HEALPix" in _table.columns, \
            "No HEALPix column in the given table. Please build it before using this method"

        # Get the initial order and update the dict
        _order: int = list(self.nodes.keys())[0]
        self.nodes[_order] = list(_table["HEALPix"])

        while _order in list(self.nodes.keys()):
            tmp = self.nodes[_order].copy()

            for ipix in tmp:
                n_ipix = ipix >> 2

                # Check the neighbouring cells
                child_nodes = [(n_ipix << 2), (n_ipix << 2) + 1, (n_ipix << 2) + 2, (n_ipix << 2) + 3]
                if all(item in self.nodes[_order] for item in child_nodes):

                    # We found 4 consecutive healpix, we can merge the neighbouring cells
                    self.add_node(_order-1, n_ipix)
                    self.del_nodes(_order, child_nodes)

            _order -= 1

    def add_node(self, _order: int, _ipix: int) -> None:
        """Add a node to the tree

        :param _order: order of the HEALPix node
        :param _ipix: value of the HEALPix node
        """

        if _order not in self.nodes:
            self.nodes[_order] = [_ipix]
        else:
            self.nodes[_order].append(_ipix)

    def del_nodes(self, _order: int, _node_list: list) -> None:
        """Del child nodes from the tree

        :param _order: order of the HEALPix node
        :param _node_list: list of neighbouring HEALPix nodes
        """

        for item in _node_list:
            self.nodes[_order].remove(item)

    def serialize_moc(self) -> str:
        """Serialize the MOC tree"""

        moc_serialization = ""

        # Sort the items of the dict to retrieve the different HEALPix orders in the good order
        for _order, ipix_list in sorted(self.nodes.items()):
            moc_serialization += serialize_healpix(ipix_list, _order)+"\n"

        return moc_serialization


def update_coordinate_system(_table: Table, _coosys: str) -> Table:
    """Update the coordinate system of the table to ICRS if needed.

    :param _table: Table from which to update the coordinates
    :param _coosys: coordinate system"""

    # Detect B1950 coordinates
    coord_system = re.match(r".*FK4.*", _coosys)
    if coord_system is not None:
        ra: str = next((col for col in _table.columns if re.match(r".*RA.*", col)), None)
        de: str = next((col for col in _table.columns if re.match(r".*DE.*", col)), None)

        # Convert B1950 coordinates to J2000
        old_coords = SkyCoord(ra=_table[ra], dec=_table[de], unit=(u.hourangle, u.deg), frame="fk4")
        new_coords = old_coords.transform_to("icrs")

        # Update the table with J2000 coordinates
        _table["RAJ2000"] = new_coords.ra
        _table["DEJ2000"] = new_coords.dec

        # Remove the columns
        _table.remove_columns([ra, de])

        print("Coordinates converted to icrs.")

    return _table


def remove_duplicates(_table: Table) -> Table:
    """Remove ipix duplicates in a table and sorts the table

    :param _table: Table from which to remove duplicates"""

    _, indices = np.unique(_table["HEALPix"], return_index=True)
    print(f"Removed {len(_table) - len(indices)} duplicates.\n")

    _table = _table[indices]
    return _table[np.where(_table["HEALPix"] != -1)]


def serialize_healpix(_healpix_table: list, _order: int) -> str:
    """Serialize the HEALPix from a HEALPix list for a given order

    :param _healpix_table: list of HEALPix values
    :param _order: order of the HEALPix values"""

    serialization = f"{_order}/"
    for healpix in _healpix_table:
        serialization += f"{healpix} "

    return serialization


def serialize_moc_tree(moc_tree: Moc_tree) -> str:
    """Serialize the MOC tree

    :param moc_tree: Moc_tree from which is to be serialized"""

    return moc_tree.serialize_moc()
