from io import StringIO
import urllib.request

from pandas import read_csv


def atcf_id(storm_id: str):
    url = 'ftp://ftp.nhc.noaa.gov/atcf/archive/storm.table'
    res = urllib.request.urlopen(url)
    df = read_csv(
        StringIO("".join([_.decode('utf-8') for _ in res])),
        header=None,
        # usecols=[]
    )
    name = f"{storm_id[:-4].upper():>10}"
    year = f"{storm_id[-4:]:>5}"
    entry = df[(df[0].isin([name]) & df[8].isin([year]))]
    if len(entry) == 0:
        return None
    else:
        return entry[20].tolist()[0].strip()
