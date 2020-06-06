import math


# all linear units in meters, angles in radians
class FoucaultTest:
    # zones is radius of each zone on couder mask, assumed to be ordered inner to outer
    # if rms_zones is set to true, the rms zone diameter will be used instead of the mean
    def __init__(self, apperture, f_ratio, zones, wavelength=5.62e-6,rms_zones=False):
        self.apperture = apperture
        self.focal_length = apperture * f_ratio
        self.radius = 2 * self.focal_length
        self.f_ratio = f_ratio
        self.wavelength = wavelength
        self.zones = zones
        self.zones_inner = [0.0] + zones[:-1]
        self.zone_means = [sqrt(i**2 + j**2) if rms_zones else (i + j) / 2 for
                         i, j in zip(self.zones_inner, zones)]
        self.disk_theoretical = f_ratio * 1.22 * wavelength
        self.transverse_coefficients = [i / (4 * self.focal_length) for i in self.zone_means]
        self.longitudinal_theoretical = [i**2 / self.radius for i in self.zone_means]  # FIXME: should be radius of curvature?

    def balanced_transverse(self, data):
        longitudinal_residual = [i - j for i, j in zip(data, self.longitudinal_theoretical)]
        transverse = [i * j for i, j in zip(longitudinal_residual, self.transverse_coefficients)]
        for i in range(len(self.zones)):
            for j in filter(lambda x: x != i, range(len(self.zones))):
                a1, a2 = longitudinal_residual[i], longitudinal_residual[j]
                b1, b2 = self.transverse_coefficients[i], self.transverse_coefficients[j]
                c = -(a1 * b1 + a2 * b2) / (b1 + b2)
                t = lambda x, y: (x + c) * y
                transverse = [t(k, l) for k, l in zip(longitudinal_residual, self.transverse_coefficients)]
                if ((max(transverse) == transverse[i] or min(transverse) == transverse[i]) and
                    (max(transverse) == transverse[j] or min(transverse) == transverse[j])):
                    return transverse

    def relative_transverse(self, data):
        return [i / self.disk_theoretical for i in self.balanced_transverse(data)]

    def slope_error(self, data):
        return [-i / self.focal_length for i in self.balanced_transverse(data)]

    def surface_error(self, data):
        res = []
        acc = 0
        se = self.slope_error(data)
        for i in range(len(self.zones)):
            acc += se[i] * (self.zones[i] - self.zones_inner[i])
            res.append(acc)
        return res


# generate zone radii for couder mask such that all zones are of equal area
def auto_zones(apperture, quantity):
    r = apperture / 2
    a = math.pi * (r ** 2)
    az = a / quantity
    return [math.sqrt((az * i) / math.pi) for i in range(1, quantity + 1)]
