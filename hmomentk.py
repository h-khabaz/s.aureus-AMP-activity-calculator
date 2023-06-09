import array
import collections
import functools
import math
import operator
import random
import statistics
import typing
import numpy

kd={'A': 1.8,'R': -4.5,'N': -3.5,'D': -3.5,'C': 2.5,'Q': -3.5,'E': -3.5,'G': -0.4,'H': -3.2,'I': 4.5,'L': 3.8,'K': -3.9,'M': 1.9,'F': 2.8,'P': -1.6,'S': -0.8,'T': -0.7,'W': -0.9,'Y': -1.3,'V': 4.2}
def hydrop_moment(self, window: int = 11, angle: int = 100) -> float:
        window = min(11, len(self))
        scale = kd
        lut = [scale.get(aa, 0.0) for aa in self._CODE1]
        angles = [(100 * i) % 360 for i in range(window)]

        if numpy is None:
            angsin = [math.sin(math.radians(theta)) for theta in angles]
            angcos = [math.cos(math.radians(theta)) for theta in angles]
        else:
            angsin = numpy.sin(numpy.radians(angles))
            angcos = numpy.cos(numpy.radians(angles))

        maxnorm = 0.0
        for i in range(len(self.sequence) - window + 1):
            # compute sin and cos of angles
            if numpy is None:
                sumsin = sumcos = 0
                for aa, s, c in zip(self.encoded[i:i+window], angsin, angcos):
                    sumsin += lut[aa]*s
                    sumcos += lut[aa]*c
            else:
                hvec = numpy.take(lut, self.encoded[i:i+window])
                sumsin = numpy.sum(hvec * angsin)
                sumcos = numpy.sum(hvec * angcos)
            # compute only the distance component (this way we can avoid
            # computing the square root in each iteration)
            norm = sumsin**2 + sumcos**2
            if norm > maxnorm:
                maxnorm = norm

        # compute the angular moment from the norm
        return math.sqrt(maxnorm) / window