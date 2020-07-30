from datetime import datetime, timedelta
import pathlib

from haversine import haversine
from matplotlib import pyplot
from matplotlib.transforms import Bbox
import numpy
from pandas import DataFrame
from pyproj import CRS, Proj
from shapely.geometry import Point, Polygon
from tropycal.realtime import Realtime

from adcircpy.forcing.winds import atcf_id
from adcircpy.forcing.winds.base import WindForcing


class NHCAdvisory(WindForcing):
    storms = Realtime()

    def __init__(self, storm_id: str, start_date: datetime = None, end_date: datetime = None,
                 crs: CRS = None):
        self._storm_id = storm_id
        super().__init__(start_date, end_date, crs)

    def clip_to_bbox(self, bbox: Bbox):
        """
        Important: bbox must be expressed in Mercator projection (EPSG:3395)
        """
        assert isinstance(bbox, Bbox), f"bbox must be a {Bbox} instance."
        bbox_pol = Polygon([
            [bbox.xmin, bbox.ymin],
            [bbox.xmax, bbox.ymin],
            [bbox.xmax, bbox.ymax],
            [bbox.xmin, bbox.ymax],
            [bbox.xmin, bbox.ymin]
        ])
        _switch = True
        unique_dates = numpy.unique(self._df['datetime'])
        for _datetime in unique_dates:
            records = self._df[self._df['datetime'] == _datetime]
            radii = records['radius_of_last_closed_isobar'].iloc[0]
            radii = 1852. * radii  # convert to meters
            merc = Proj("EPSG:3395")
            x, y = merc(
                records['longitude'].iloc[0],
                records['latitude'].iloc[0])
            p = Point(x, y)
            pol = p.buffer(radii)
            if _switch:
                if not pol.intersects(bbox_pol):
                    continue
                else:
                    self.start_date = records['datetime'].iloc[0]
                    _switch = False
                    continue
                    # self.start_date = 
            else:
                if pol.intersects(bbox_pol):
                    continue
                else:
                    self.end_date = records['datetime'].iloc[0]
                    break

    def plot_trajectory(self, ax: pyplot.Axes = None, show: bool = False, color='k', **kwargs):
        kwargs.update({'color': color})
        if ax is None:
            fig = pyplot.figure()
            ax = fig.add_subplot(111)
        for i in range(len(self.speed)):
            # when dealing with nautical degrees, U is sine and V is cosine.
            U = self.speed.iloc[i] * numpy.sin(numpy.deg2rad(self.direction.iloc[i]))
            V = self.speed.iloc[i] * numpy.cos(numpy.deg2rad(self.direction.iloc[i]))
            ax.quiver(
                self.longitude.iloc[i], self.latitude.iloc[i], U, V, **kwargs)
            ax.annotate(
                self.df['datetime'].iloc[i],
                (self.longitude.iloc[i], self.latitude.iloc[i])
            )
        if show:
            ax.axis('scaled')
            pyplot.show()

    def write(self, path: str, overwrite: bool = False):
        path = pathlib.Path(path)
        if path.is_file() and not overwrite:
            raise Exception('Files exist, set overwrite=True to allow overwrite.')
        with open(path, 'w') as f:
            f.write(self.fort22)

    @property
    def storm_id(self):
        return self._storm_id

    @property
    def _storm_id(self):
        return f"{self.basin}{self.storm_number}{self.year}"

    @_storm_id.setter
    def _storm_id(self, storm_id):
        chars = 0
        for char in storm_id:
            if char.isdigit():
                chars += 1

        if chars == 4:
            _atcf_id = atcf_id(storm_id)
            if _atcf_id is None:
                raise Exception(f'No storm with id: {storm_id}')
            storm_id = _atcf_id

        self.storm = self.storms.get_storm(storm_id)

    @property
    def _start_date(self):
        return self.__start_date

    @_start_date.setter
    def _start_date(self, start_date):
        if start_date is not None:
            assert isinstance(start_date, datetime)
        else:
            start_date = self._df['datetime'].iloc[0]
        assert self._df['datetime'].iloc[0] <= start_date < self._df['datetime'].iloc[-1], \
            f"start_date must be {self._df['datetime'].iloc[0]} <= start_date ({start_date}) < " \
            f"{self._df['datetime'].iloc[-1]}"
        self.__start_date = start_date

    @property
    def _end_date(self):
        return self.__end_date

    @_end_date.setter
    def _end_date(self, end_date):
        if end_date is not None:
            assert isinstance(end_date, datetime)
        else:
            end_date = self._df['datetime'].iloc[-1]
        assert self._df['datetime'].iloc[0] < end_date <= self._df['datetime'].iloc[-1], \
            f"end_date must be {self._df['datetime'].iloc[0]} <= end_date ({end_date}) <= " \
            f"{self._df['datetime'].iloc[-1]}"
        assert end_date > self.start_date, \
            f"end_date ({end_date}) must be after start_date ({self.start_date})"
        self.__end_date = end_date

    @property
    def name(self):
        return self.df['name'].value_counts()[:].index.tolist()[0]

    @property
    def basin(self):
        return self.df['basin'].iloc[0]

    @property
    def storm_number(self):
        return self.df['storm_number'].iloc[0]

    @property
    def year(self):
        return self.df['datetime'].iloc[0].year

    @property
    def datetime(self):
        return self.df['datetime']

    @property
    def speed(self):
        return self.df['speed']

    @property
    def direction(self):
        return self.df['direction']

    @property
    def longitude(self):
        return self.df['longitude']

    @property
    def latitude(self):
        return self.df['latitude']

    @property
    def df(self):
        return self._df[(self._df['datetime'] >= self.start_date) &
                        (self._df['datetime'] <= self._file_end_date)]

    @property
    def _df(self):
        # https://www.nrlmry.navy.mil/atcf_web/docs/database/new/abdeck.txt
        try:
            return self.__df
        except AttributeError:
            data = {
                "basin"                       : self.storm['wmo_basin'],
                "storm_number"                : int(self.storm['id'][2:4]),
                "datetime"                    : self.storm['date'],
                "record_type"                 : self.storm['special'],
                "latitude"                    : self.storm['lat'],
                "longitude"                   : self.storm['lon'],
                "max_sustained_wind_speed"    : self.storm['vmax'],
                "central_pressure"            : self.storm['mslp'],
                "development_level"           : self.storm['type'],
                "isotach"                     : None,
                "quadrant"                    : None,
                "radius_for_NEQ"              : None,
                "radius_for_SEQ"              : None,
                "radius_for_SWQ"              : None,
                "radius_for_NWQ"              : None,
                "background_pressure"         : None,
                "radius_of_last_closed_isobar": None,
                "radius_of_maximum_winds"     : None,
                "name"                        : self.storm['name'],
                "direction"                   : [],
                "speed"                       : []
            }
            data = self._compute_velocity(data)
            # data = self._transform_coordinates(data)
            self.__df = DataFrame(data=data)
            return self.__df

    @property
    def fort22(self):
        record_number = self._generate_record_numbers()
        fort22 = ''
        for i, (_, row) in enumerate(self.df.iterrows()):
            latitude = row['latitude']
            if latitude >= 0:
                latitude = f'{int(latitude / 0.1):>4}N'
            else:
                latitude = f'{int(latitude / -0.1):>4}S'
            longitude = row['longitude']
            if longitude >= 0:
                longitude = f'{int(longitude / 0.1):>5}E'
            else:
                longitude = f'{int(longitude / -0.1):>5}W'

            background_pressure = row["background_pressure"]
            if background_pressure is None:
                background_pressure = self.df["background_pressure"].iloc[i - 1]
            if background_pressure is not None:
                if background_pressure <= row["central_pressure"] < 1013:
                    background_pressure = 1013
                elif background_pressure <= row["central_pressure"] >= 1013:
                    background_pressure = int(row["central_pressure"] + 1)
                else:
                    background_pressure = int(row["background_pressure"])
            else:
                background_pressure = ""

            row.extend([
                f'{row["basin"]:<2},',
                f'{row["storm_number"]:>3},',
                f'{format(row["datetime"], "%Y%m%d%H"):>11},',
                f'{"":3},',
                f'{row["record_type"]:>5},',
                f'{int((row["datetime"] - self.start_date) / timedelta(hours=1)):>4},',
                f'{latitude:>5},',
                f'{longitude:>5},',
                f'{int(row["max_sustained_wind_speed"]):>4},',
                f'{int(row["central_pressure"]):>5},',
                f'{row["development_level"] if row["development_level"] is not None else "":>3},',
                f'{int(row["isotach"]) if row["isotach"] is not None else "":>4},',
                f'{row["quadrant"] if row["quadrant"] is not None else "":>4},',
                f'{int(row["radius_for_NEQ"]) if row["radius_for_NEQ"] is not None else "":>5},',
                f'{int(row["radius_for_SEQ"]) if row["radius_for_SEQ"] is not None else "":>5},',
                f'{int(row["radius_for_SWQ"]) if row["radius_for_SWQ"] is not None else "":>5},',
                f'{int(row["radius_for_NWQ"]) if row["radius_for_NWQ"] is not None else "":>5},',
                f'{background_pressure:>5},',
                f'{int(row["radius_of_last_closed_isobar"]) if row["radius_of_last_closed_isobar"] is not None else "":>5},',
                f'{int(row["radius_of_maximum_winds"]) if row["radius_of_maximum_winds"] is not None else "":>4},'
                f'{"":>5},',  # gust
                f'{"":>4},',  # eye
                f'{"":>4},',  # subregion
                f'{"":>4},',  # maxseas
                f'{"":>4},',  # initials
                f'{row["direction"] if row["direction"] is not None else "":>3},',
                f'{row["speed"]:>4},',
                f'{row["name"]:^12},',
                # from this point forwards it's all aswip
                f'{record_number[i]:>4}'
            ])
            row = ', '.join(row)
            fort22 += f'{row}\n'
        return fort22

    @property
    def WTIMINC(self):
        return f'{self.start_date:%Y %m %d %H} ' \
               f'{self.df["storm_number"].iloc[0]} {self.BLADj} {self.geofactor}'

    @property
    def BLADj(self):
        try:
            return self.__BLADj
        except AttributeError:
            return 0.9

    @BLADj.setter
    def BLADj(self, BLADj: float):
        BLADj = float(BLADj)
        assert 0 <= BLADj <= 1
        self.__BLADj = BLADj

    @property
    def geofactor(self):
        try:
            return self.__geofactor
        except AttributeError:
            return 1

    @geofactor.setter
    def geofactor(self, geofactor: float):
        geofactor = float(geofactor)
        assert 0 <= geofactor <= 1
        self.__geofactor = geofactor

    def _generate_record_numbers(self):
        record_number = [1]
        for i in range(1, len(self.datetime)):
            if self.datetime.iloc[i] == self.datetime.iloc[i - 1]:
                record_number.append(record_number[-1])
            else:
                record_number.append(record_number[-1] + 1)
        return record_number

    @property
    def _file_end_date(self):
        unique_dates = numpy.unique(self._df['datetime'])
        for date in unique_dates:
            if date >= self.end_date:
                return date

    @staticmethod
    def _compute_velocity(data: {}):
        """
        Output has units of meters per second.
        """

        merc = Proj("EPSG:3395")
        x, y = merc(data['longitude'], data['latitude'])
        t = data['datetime']
        unique_datetimes = numpy.unique(t)
        for i, _datetime in enumerate(unique_datetimes):
            indexes, = numpy.where(numpy.asarray(t) == _datetime)
            for idx in indexes:
                if indexes[-1] + 1 < len(t):
                    dx = haversine((y[idx], x[indexes[-1] + 1]), (y[idx], x[idx]), unit='nmi')
                    dy = haversine((y[indexes[-1] + 1], x[idx]), (y[idx], x[idx]), unit='nmi')
                    dt = ((t[indexes[-1] + 1] - t[idx]) / timedelta(hours=1))
                    vx = numpy.copysign(dx / dt, x[indexes[-1] + 1] - x[idx])
                    vy = numpy.copysign(dy / dt, y[indexes[-1] + 1] - y[idx])
                else:
                    dx = haversine((y[idx], x[indexes[0] - 1]), (y[idx], x[idx]), unit='nmi')
                    dy = haversine((y[indexes[0] - 1], x[idx]), (y[idx], x[idx]), unit='nmi')
                    dt = ((t[idx] - t[indexes[0] - 1]) / timedelta(hours=1))
                    vx = numpy.copysign(dx / dt, x[idx] - x[indexes[0] - 1])
                    vy = numpy.copysign(dy / dt, y[idx] - y[indexes[0] - 1])
                speed = numpy.sqrt(dx ** 2 + dy ** 2) / dt
                bearing = (360. + numpy.rad2deg(numpy.arctan2(vx, vy))) % 360
                data['speed'].append(int(numpy.around(speed, 0)))
                data['direction'].append(int(numpy.around(bearing, 0)))
        return data
