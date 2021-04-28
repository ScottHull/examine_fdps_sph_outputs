import pandas as pd

path = "/Users/scotthull/Desktop/2000.csv"

df = pd.read_csv(path)

EARTH_MASS = 5.972 * 10 ** 24
LUNAR_MASS = 7.34767309 * 10 ** 22

data = {
    "DISK": {
        "mass": 0,
        "particles": 0
    },
    "PLANET": {
        "mass": 0,
        "particles": 0
    },
    "ESCAPE": {
        "mass": 0,
        "particles": 0
    },
}

for row in df.index:
    label = df['label'][row]
    mass = df['mass'][row]

    data[label]['particles'] += 1
    data[label]['mass'] += mass

print(
    "DISK: {} PARTICLES, {} MASS ({} LUNAR MASSES)\n"
    "PLANET: {} PARTICLES, {} MASS ({} EARTH MASSES)\n"
    "ESCAPE: {} PARTICLES, {} MASS ({} LUNAR MASSES)\n".format(
        data['DISK']['particles'], data['DISK']['mass'], data['DISK']['mass'] / LUNAR_MASS,
        data['PLANET']['particles'], data['PLANET']['mass'], data['PLANET']['mass'] / EARTH_MASS,
        data['ESCAPE']['particles'], data['ESCAPE']['mass'], data['ESCAPE']['mass'] / LUNAR_MASS,
    )
)
