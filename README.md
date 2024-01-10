# MOClib - A custom MOC library
__Author:__ Yannick Hénin<br>
__University of Strasbourg__

## Table of Content
- [Description](#description)
- [Examples](#examples)
- [Installation](#installation)

## Description

Building a Python MOC library that creates a MOC from list of HEALPix cells.  Note: the project asks to create a library and to illustrate it into a notebook. The best way to provide a library, consists to  make a module that can be reused. As a bonus, I also package the module

## Example

A Jupyter Notebook can be found in the package.<br>
It provides more details on how to use the library and what can be done with it.<br>
First, make sure you have installed jupyter on your python environnement: `pip install jupyter` before running the notebook.


## Installation

To install the package, please use pip or conda :

```shell
> pip install MOClib/
OR
> conda install moc_library/
```

You can already try to load ipyaladin in a notebook.

```python
from moc_library import Loader, Moc_tree

loader = Loader(_catalogue= ,_out_max=10000, _order=7, _source="CDS")
votable = loader.get_votable()
...
moc_tree = Moc_tree(_order=7)
moc_tree.build_moc_tree(votable)
...
```
