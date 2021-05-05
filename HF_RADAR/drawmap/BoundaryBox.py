class BoundaryBox:
    lat_min: float
    lat_max: float
    lon_min: float
    lon_max: float

    def __init__(self):
        self.lat_min = -90.
        self.lat_max = 90.
        self.lon_min = -180.
        self.lon_max = 180.

    def middle_lat(self):
        x = self.lat_min + 0.5 * (self.lat_max - self.lat_min)
        print(x)
        return x

    def middle_lon(self):
        return self.lon_min + 0.5 * (self.lon_max - self.lon_min)


