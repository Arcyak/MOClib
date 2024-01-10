# MOClib - A custom MOC library
__Author:__ Yannick Hénin<br>
__University of Strasbourg__

## Table of Content
- [Description](#description)
- [Examples](#examples)
- [Installation](#installation)

## Description

MOClib is a custom Python MOC library that uses the Astropy library to build MOC trees and serialization from HEALPix cells.
The library contains actually two classes: `python Loader and Moc_tree`. <br>
The Loader class is able to query a VOTable from Vizier, update the coordinate system if needed and create the HEALPix number associated to each object.
The Moc_tree class takes care of processing the cells from the HEALPix column of a astropy.Table object. A MOC serialization can be obtained from the Moc_tree class for further analysis.


## Example

A Jupyter Notebook can be found in the package.<br>
It provides more details on how to use the library and what can be done with it.<br>
First, make sure you have installed jupyter on your python environnement: `pip install jupyter` before running the notebook.


## Installation

To install the package, please use pip or conda :

```shell
> pip install MOClib/
OR
> conda install MOClib/
```

You can already try to load ipyaladin in a notebook.

```python
from moc_library import Loader, Moc_tree

loader = Loader(_catalogue='II/7A/catalog', _out_max=10000, _order=7, _source="CDS")
votable = loader.get_votable()
...
moc_tree = Moc_tree(_order=7)
moc_tree.build_moc_tree(votable)
...
```
