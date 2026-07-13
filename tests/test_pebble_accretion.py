import pytest
import numpy as np
from pa3py.pebble_accretion import PebbleAccretionModule3
from pa3py import constants as c

def test_constants():
    """Verifica que las constantes estén bien definidas en CGS."""
    assert c.G == 6.6743e-8
    assert c.M_SUN > 1.9e33
    assert c.M_EARTH > 5.9e27
    assert c.AU > 1.4e13

def test_module_initialization():
    """Prueba básica de existencia de la clase."""
    assert PebbleAccretionModule3 is not None
    assert PebbleAccretionModule3.G_CGS == c.G
