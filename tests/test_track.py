import unittest

import numpy
from pyproj import CRS

from adcircpy.forcing.winds.best_track import ellipsoidal_distance


class TestTrack(unittest.TestCase):
    def test_distance(self):
        point_a = (-77.035, 38.889)
        point_b = (325750.9, 4306472.6)
        crs_a = CRS.from_epsg(4326)
        crs_b = CRS.from_epsg(26918)

        distance = ellipsoidal_distance(point_a, point_b, crs_a, crs_b)

        assert numpy.allclose(distance, 14823241.794380495)


if __name__ == '__main__':
    unittest.main()
